from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "mechanism_drug_work/fda_repurposing_work/cache"
WORK.mkdir(parents=True, exist_ok=True)
DATA = (
    ROOT
    / "mechanism_drug_work/reference_data/oncoPredict/full/verified_extracted/"
    "DataFiles/DataFiles"
)
REG = ROOT / "mechanism_drug_work/reference_data/regulatory_validation"

SEEDS = (17, 43, 97)
MIN_N = 100
MIN_RESPONSE_SD = 0.5
MIN_MEDIAN_SPEARMAN = 0.25
MIN_POSITIVE_FRACTION = 0.80
RIDGE_ALPHAS = np.logspace(-2, 4, 13)
N_BOOTSTRAP = 100


@dataclass
class HN6Input:
    tpm: pd.DataFrame
    expression: pd.DataFrame
    groups: pd.Series
    sample_audit: pd.DataFrame


def normalize_name(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).lower().strip()
    text = text.replace("5-fu", "5 fluorouracil")
    text = re.sub(r"[._/\\-]+", " ", text)
    text = re.sub(r"[^a-z0-9+ ]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


SALT_SUFFIXES = {
    "hydrochloride",
    "hydrobromide",
    "mesylate",
    "maleate",
    "fumarate",
    "tartrate",
    "succinate",
    "acetate",
    "citrate",
    "phosphate",
    "sulfate",
    "sulphate",
    "nitrate",
    "besylate",
    "tosylate",
    "malate",
    "hydrate",
    "monohydrate",
    "dihydrate",
}


def strip_salt(value: object) -> str:
    name = normalize_name(value)
    bits = name.split()
    while len(bits) > 1 and bits[-1] in SALT_SUFFIXES:
        bits.pop()
    return " ".join(bits)


def load_hn6(path: Path) -> HN6Input:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    tpm_cols = [
        c
        for c in raw.columns
        if c.startswith("tpm_HN6_") and not c.startswith("average_")
    ]
    if len(tpm_cols) != 12:
        raise RuntimeError(f"Expected 12 TPM columns, found {len(tpm_cols)}")
    if raw["gene_symbol"].isna().any() or raw["gene_symbol"].duplicated().any():
        raise RuntimeError("Gene symbols are missing or duplicated")
    tpm = raw.set_index("gene_symbol")[tpm_cols].astype(float)
    sample_names = [c.removeprefix("tpm_").replace("Vetor", "Vector") for c in tpm_cols]
    tpm.columns = sample_names
    groups = []
    for sample in sample_names:
        if "_KD" in sample:
            groups.append("KD")
        elif "_NC" in sample:
            groups.append("NC")
        elif "_OE" in sample:
            groups.append("OE")
        elif "_Vector" in sample:
            groups.append("Vector")
        else:
            raise RuntimeError(f"Unrecognized HN6 sample: {sample}")
    group_series = pd.Series(groups, index=sample_names, name="group")
    if group_series.value_counts().to_dict() != {
        "KD": 3,
        "NC": 3,
        "OE": 3,
        "Vector": 3,
    }:
        raise RuntimeError(f"Unexpected group counts: {group_series.value_counts().to_dict()}")
    if "PDE2A" not in tpm.index:
        raise RuntimeError("PDE2A is absent from TPM matrix")
    pde2a = tpm.loc["PDE2A"]
    means = pde2a.groupby(group_series).mean()
    if not (means["KD"] < means["NC"] and means["OE"] > means["Vector"]):
        raise RuntimeError(f"PDE2A perturbation direction failed: {means.to_dict()}")
    expression = np.log2(tpm + 1.0)
    sample_audit = pd.DataFrame(
        {
            "sample": sample_names,
            "raw_column": tpm_cols,
            "group": group_series.values,
            "tpm_sum": tpm.sum(axis=0).values,
            "PDE2A_TPM": pde2a.values,
            "PDE2A_log2TPM1": expression.loc["PDE2A"].values,
        }
    )
    return HN6Input(tpm, expression, group_series, sample_audit)


def read_rds_matrix(path: Path) -> pd.DataFrame:
    result = pyreadr.read_r(str(path))
    if len(result) != 1:
        raise RuntimeError(f"Expected one R object in {path}")
    return next(iter(result.values()))


def load_fda_names() -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    drugs = pd.read_csv(REG / "drugsatfda/Products.txt", sep="\t", encoding="latin1")
    orange = pd.read_csv(REG / "orangebook/products.txt", sep="~", encoding="latin1")
    rows = []
    for source, frame, ingredient_col, trade_col in [
        ("Drugs@FDA", drugs, "ActiveIngredient", "DrugName"),
        ("OrangeBook", orange, "Ingredient", "Trade_Name"),
    ]:
        for _, r in frame[[ingredient_col, trade_col]].dropna(subset=[ingredient_col]).iterrows():
            ingredient = str(r[ingredient_col]).strip()
            if ";" in ingredient:
                continue
            rows.append(
                {
                    "source": source,
                    "active_ingredient": ingredient,
                    "trade_name": str(r.get(trade_col, "")),
                    "normalized": normalize_name(ingredient),
                    "base_normalized": strip_salt(ingredient),
                }
            )
    table = pd.DataFrame(rows).drop_duplicates()
    exact_map = (
        table.sort_values(["source", "active_ingredient"])
        .drop_duplicates("normalized")
        .set_index("normalized")["active_ingredient"]
        .to_dict()
    )
    base_map = (
        table.sort_values(["source", "active_ingredient"])
        .drop_duplicates("base_normalized")
        .set_index("base_normalized")["active_ingredient"]
        .to_dict()
    )
    return table, exact_map, base_map


def split_response_name(column: str) -> tuple[str, str]:
    match = re.match(r"^(.*)_([0-9]+)$", column)
    if not match:
        return column, ""
    return match.group(1), match.group(2)


def build_database_drug_table() -> pd.DataFrame:
    _, fda_exact, fda_base = load_fda_names()
    ann = pd.read_csv(DATA / "GLDS/GDSCv1/screened_compunds_rel_8.2.csv")
    ann["DRUG_ID"] = ann["DRUG_ID"].astype(str)
    ann_by_id = ann.set_index("DRUG_ID")
    prism = pd.read_csv(REG / "prism_secondary_treatment_info.csv", low_memory=False)
    prism = prism.sort_values("dose").drop_duplicates("broad_id")
    prism_names = set(prism["name"].map(normalize_name))
    prism_ids = prism.groupby(prism["name"].map(normalize_name))["broad_id"].apply(
        lambda x: ";".join(sorted(set(map(str, x))))
    )
    ctrp = read_rds_matrix(DATA / "Training Data/CTRP2_Res.rds")
    ctrp_names = {normalize_name(c): c for c in ctrp.columns}
    rows = []
    for database, response_file in [
        ("GDSC1", "GDSC1_Res.rds"),
        ("GDSC2", "GDSC2_Res.rds"),
    ]:
        response = read_rds_matrix(DATA / f"Training Data/{response_file}")
        for col in response.columns:
            raw_name, drug_id = split_response_name(str(col))
            annotation = ann_by_id.loc[drug_id] if drug_id in ann_by_id.index else None
            canonical = str(annotation["DRUG_NAME"]) if annotation is not None else raw_name
            synonyms = "" if annotation is None or pd.isna(annotation["SYNONYMS"]) else str(annotation["SYNONYMS"])
            aliases = [canonical, raw_name] + [x.strip() for x in synonyms.split(",") if x.strip()]
            exact_hits = [a for a in aliases if normalize_name(a) in fda_exact]
            base_hits = [a for a in aliases if strip_salt(a) in fda_base]
            if exact_hits:
                fda_status, method = True, "FDA active ingredient exact alias"
                active_ingredient = fda_exact[normalize_name(exact_hits[0])]
            elif base_hits:
                fda_status, method = True, "FDA active ingredient salt/base alias"
                active_ingredient = fda_base[strip_salt(base_hits[0])]
            else:
                fda_status, method = False, "no conservative FDA active-ingredient match"
                active_ingredient = ""
            is_small_molecule = not bool(
                re.search(r"(?:mab|cept)$", normalize_name(active_ingredient or canonical))
                or normalize_name(canonical) in {"cetuximab"}
            )
            if fda_status and not is_small_molecule:
                method = "excluded FDA biologic"
            norm = normalize_name(canonical)
            rows.append(
                {
                    "database": database,
                    "response_column": col,
                    "drug_id": drug_id,
                    "drug": canonical,
                    "normalized_drug": norm,
                    "synonyms": synonyms,
                    "known_target": "" if annotation is None else annotation["TARGET"],
                    "target_pathway": "" if annotation is None else annotation["TARGET_PATHWAY"],
                    "FDA_status": fda_status,
                    "FDA_match_method": method,
                    "active_ingredient": active_ingredient,
                    "small_molecule_status": is_small_molecule,
                    "PRISM_available": norm in prism_names,
                    "PRISM_ID": prism_ids.get(norm, ""),
                    "CTRP_available": norm in ctrp_names,
                    "CTRP_response_column": ctrp_names.get(norm, ""),
                }
            )
    for col in ctrp.columns:
        norm = normalize_name(col)
        base = strip_salt(col)
        if norm in fda_exact:
            fda_status, method = True, "FDA active ingredient exact name"
            active_ingredient = fda_exact[norm]
        elif base in fda_base:
            fda_status, method = True, "FDA active ingredient salt/base name"
            active_ingredient = fda_base[base]
        else:
            fda_status, method, active_ingredient = (
                False,
                "no conservative FDA active-ingredient match",
                "",
            )
        is_small_molecule = not bool(
            re.search(r"(?:mab|cept)$", normalize_name(active_ingredient or col))
        )
        rows.append(
            {
                "database": "CTRP2",
                "response_column": col,
                "drug_id": col,
                "drug": col,
                "normalized_drug": norm,
                "synonyms": "",
                "known_target": "",
                "target_pathway": "",
                "FDA_status": fda_status,
                "FDA_match_method": method,
                "active_ingredient": active_ingredient,
                "small_molecule_status": is_small_molecule,
                "PRISM_available": norm in prism_names,
                "PRISM_ID": prism_ids.get(norm, ""),
                "CTRP_available": True,
                "CTRP_response_column": col,
            }
        )
    return pd.DataFrame(rows)


def prepare_expression(
    training: pd.DataFrame, hn6: HN6Input, max_pcs: int = 100
) -> tuple[np.ndarray, np.ndarray, list[str], PCA, dict]:
    common = training.index.intersection(hn6.expression.index)
    train = training.loc[common].astype(np.float32)
    test = hn6.expression.loc[common].astype(np.float32)
    means = train.mean(axis=1)
    sds = train.std(axis=1, ddof=1)
    keep = (means > 3.0) & (sds > 0.25)
    train = train.loc[keep]
    test = test.loc[keep]
    train_z = train.sub(train.mean(axis=1), axis=0).div(train.std(axis=1, ddof=1), axis=0)
    test_sd = test.std(axis=1, ddof=1).replace(0, np.nan)
    test_z = test.sub(test.mean(axis=1), axis=0).div(test_sd, axis=0).fillna(0.0)
    n_components = min(max_pcs, train_z.shape[1] - 1, train_z.shape[0])
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=17)
    train_pc = pca.fit_transform(train_z.T)
    test_pc = pca.transform(test_z.T)
    meta = {
        "n_common_genes": len(common),
        "n_filtered_genes": train_z.shape[0],
        "n_pcs": n_components,
        "explained_variance": float(pca.explained_variance_ratio_.sum()),
    }
    return train_pc, test_pc, list(train.columns), pca, meta


def ood_metrics(train_pc: np.ndarray, test_pc: np.ndarray, k: int = 5) -> dict:
    dims = min(30, train_pc.shape[1])
    tr = train_pc[:, :dims]
    te = test_pc[:, :dims]
    within = cdist(tr, tr)
    np.fill_diagonal(within, np.inf)
    train_knn = np.partition(within, k - 1, axis=1)[:, :k].mean(axis=1)
    test_dist = cdist(te, tr)
    test_knn = np.partition(test_dist, k - 1, axis=1)[:, :k].mean(axis=1)
    lw = LedoitWolf().fit(tr)
    train_maha = np.sqrt(lw.mahalanobis(tr))
    test_maha = np.sqrt(lw.mahalanobis(te))
    return {
        "train_knn_p95": float(np.quantile(train_knn, 0.95)),
        "test_knn_median": float(np.median(test_knn)),
        "test_knn_max": float(np.max(test_knn)),
        "train_mahalanobis_p975": float(np.quantile(train_maha, 0.975)),
        "test_mahalanobis_median": float(np.median(test_maha)),
        "test_mahalanobis_max": float(np.max(test_maha)),
        "ood_flag": bool(
            np.median(test_knn) > np.quantile(train_knn, 0.95)
            or np.median(test_maha) > np.quantile(train_maha, 0.975)
        ),
    }


def group_delta(pred: np.ndarray, groups: pd.Series) -> tuple[float, float]:
    s = pd.Series(pred, index=groups.index)
    return (
        float(s[groups == "KD"].mean() - s[groups == "NC"].mean()),
        float(s[groups == "OE"].mean() - s[groups == "Vector"].mean()),
    )


def fit_drug_model(
    X: np.ndarray,
    y: np.ndarray,
    test_pc: np.ndarray,
    groups: pd.Series,
    seed_offset: int,
) -> dict:
    valid = np.isfinite(y)
    Xv, yv = X[valid], y[valid].astype(float)
    n = len(yv)
    response_sd = float(np.std(yv, ddof=1))
    if n < MIN_N or response_sd < MIN_RESPONSE_SD:
        return {
            "n_train": n,
            "response_sd": response_sd,
            "model_pass": False,
            "exclusion_reason": "n_train<100" if n < MIN_N else "response_sd<0.5",
        }
    fold_metrics = []
    for seed in SEEDS:
        kf = KFold(n_splits=5, shuffle=True, random_state=seed + seed_offset)
        for fold, (tr, va) in enumerate(kf.split(Xv), 1):
            model = RidgeCV(alphas=RIDGE_ALPHAS, cv=None)
            model.fit(Xv[tr], yv[tr])
            pred = model.predict(Xv[va])
            rho = spearmanr(yv[va], pred).statistic
            fold_metrics.append(
                {
                    "seed": seed,
                    "fold": fold,
                    "spearman": float(rho) if np.isfinite(rho) else 0.0,
                    "rmse": float(mean_squared_error(yv[va], pred) ** 0.5),
                    "alpha": float(model.alpha_),
                }
            )
    folds = pd.DataFrame(fold_metrics)
    median_rho = float(folds["spearman"].median())
    positive_fraction = float((folds["spearman"] > 0).mean())
    model_pass = median_rho >= MIN_MEDIAN_SPEARMAN and positive_fraction >= MIN_POSITIVE_FRACTION
    result = {
        "n_train": n,
        "response_sd": response_sd,
        "cv_median_spearman": median_rho,
        "cv_mean_rmse": float(folds["rmse"].mean()),
        "cv_positive_fraction": positive_fraction,
        "model_pass": model_pass,
        "exclusion_reason": "" if model_pass else "CV quality gate failed",
        "cv_fold_records": fold_metrics,
    }
    if not model_pass:
        return result
    full = RidgeCV(alphas=RIDGE_ALPHAS, cv=None).fit(Xv, yv)
    point = full.predict(test_pc)
    delta_kd, delta_oe = group_delta(point, groups)
    rng = np.random.default_rng(1000 + seed_offset)
    boot_kd, boot_oe = [], []
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, n)
        boot = Ridge(alpha=float(full.alpha_)).fit(Xv[idx], yv[idx])
        pred = boot.predict(test_pc)
        a, b = group_delta(pred, groups)
        boot_kd.append(a)
        boot_oe.append(b)
    boot_kd = np.asarray(boot_kd)
    boot_oe = np.asarray(boot_oe)
    result.update(
        {
            "alpha": float(full.alpha_),
            "delta_KD": delta_kd,
            "delta_OE": delta_oe,
            "std_delta_KD": delta_kd / response_sd,
            "std_delta_OE": delta_oe / response_sd,
            "bootstrap_KD_low": float(np.quantile(boot_kd, 0.025)),
            "bootstrap_KD_high": float(np.quantile(boot_kd, 0.975)),
            "bootstrap_OE_low": float(np.quantile(boot_oe, 0.025)),
            "bootstrap_OE_high": float(np.quantile(boot_oe, 0.975)),
            "bootstrap_KD_negative_fraction": float((boot_kd < 0).mean()),
            "bootstrap_KD_positive_fraction": float((boot_kd > 0).mean()),
            "bootstrap_OE_negative_fraction": float((boot_oe < 0).mean()),
            "bootstrap_OE_positive_fraction": float((boot_oe > 0).mean()),
            "predictions": dict(zip(groups.index, map(float, point))),
        }
    )
    return result


def run_database(
    database: str,
    expression_file: str,
    response_file: str,
    drug_table: pd.DataFrame,
    hn6: HN6Input,
    log2_training: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, np.ndarray, np.ndarray]:
    expression = read_rds_matrix(DATA / f"Training Data/{expression_file}")
    if log2_training:
        expression = np.log2(expression.astype(np.float32) + 1.0)
    response = read_rds_matrix(DATA / f"Training Data/{response_file}")
    response = response.reindex(expression.columns)
    train_pc, test_pc, samples, _, prep = prepare_expression(expression, hn6)
    ood = ood_metrics(train_pc, test_pc)
    selected = drug_table[
        (drug_table["database"] == database)
        & drug_table["FDA_status"]
        & drug_table["small_molecule_status"]
    ]
    rows, sample_rows = [], []
    for offset, r in enumerate(selected.itertuples(index=False)):
        fit = fit_drug_model(
            train_pc,
            response[r.response_column].to_numpy(float),
            test_pc,
            hn6.groups,
            offset,
        )
        base = r._asdict()
        records = fit.pop("cv_fold_records", [])
        predictions = fit.pop("predictions", {})
        base.update(fit)
        rows.append(base)
        for sample, value in predictions.items():
            sample_rows.append(
                {
                    "database": database,
                    "drug": r.drug,
                    "drug_id": r.drug_id,
                    "sample": sample,
                    "group": hn6.groups[sample],
                    "predicted_response": value,
                    "lower_is_more_sensitive": True,
                }
            )
    model_table = pd.DataFrame(rows)
    sample_table = pd.DataFrame(sample_rows)
    metadata = {"database": database, **prep, **ood, "n_fda_models": len(selected)}
    return model_table, sample_table, metadata, train_pc, test_pc


def run() -> None:
    hn6 = load_hn6(ROOT / "mechanism_drug_work/vendor_report_extract/gene_expression.xls")
    drug_table = build_database_drug_table()
    drug_table.to_csv(WORK / "database_drug_mapping.tsv", sep="\t", index=False)
    hn6.sample_audit.to_csv(WORK / "hn6_sample_audit.tsv", sep="\t", index=False)
    all_models, all_samples, metadata = [], [], []
    pca_arrays = {}
    for db, expr, res in [
        (
            "GDSC1",
            "GDSC1_Expr (RMA Normalized and Log Transformed).rds",
            "GDSC1_Res.rds",
        ),
        (
            "GDSC2",
            "GDSC2_Expr (RMA Normalized and Log Transformed).rds",
            "GDSC2_Res.rds",
        ),
        (
            "CTRP2",
            "CTRP2_Expr (TPM, not log transformed).rds",
            "CTRP2_Res.rds",
        ),
    ]:
        models, samples, meta, train_pc, test_pc = run_database(
            db, expr, res, drug_table, hn6, log2_training=(db == "CTRP2")
        )
        all_models.append(models)
        all_samples.append(samples)
        metadata.append(meta)
        pca_arrays[f"{db}_train"] = train_pc[:, :10]
        pca_arrays[f"{db}_test"] = test_pc[:, :10]
    pd.concat(all_models, ignore_index=True).to_csv(
        WORK / "all_model_results.tsv", sep="\t", index=False
    )
    pd.concat(all_samples, ignore_index=True).to_csv(
        WORK / "all_sample_predictions.tsv", sep="\t", index=False
    )
    (WORK / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    np.savez_compressed(WORK / "pca_coordinates.npz", **pca_arrays)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    run()
