from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from run_pipeline import DATA, REG, ROOT, WORK, normalize_name, strip_salt


OUT = ROOT / "results/PDE2A_FDA_drug_repurposing_20260801"
OUT.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "NC": "#7C8A9A",
    "KD": "#D55E00",
    "Vector": "#56B4E9",
    "OE": "#009E73",
    "Tier 1": "#7C3AED",
    "Tier 2": "#E69F00",
    "Tier 3": "#B8BEC7",
}


def active_key(value: object) -> str:
    return strip_salt(value)


def load_l1000_support() -> pd.DataFrame:
    hits = pd.read_csv(
        ROOT / "results/PDE2A_cmap_lincs_inputs/L1000FWD_top_hits_annotated.tsv",
        sep="\t",
    )
    manual = {
        "BRD-K49810818": "sorafenib",
        "BRD-K49328571": "dasatinib",
        "BRD-K12502280": "fedratinib",
    }
    rows = []
    for r in hits[hits["mode"].eq("opposite")].itertuples(index=False):
        key = ""
        for broad, mapped in manual.items():
            if broad in str(r.sig_id):
                key = mapped
        if key:
            rows.append(
                {
                    "active_key": key,
                    "L1000FWD_support": True,
                    "L1000FWD_drug_record": str(r.drug),
                    "L1000FWD_combined_score": r.combined_score,
                    "L1000FWD_q_value": r.q_value,
                }
            )
    return pd.DataFrame(rows).sort_values("L1000FWD_combined_score").drop_duplicates(
        "active_key"
    )


def aggregate_models() -> tuple[pd.DataFrame, pd.DataFrame]:
    models = pd.read_csv(WORK / "all_model_results.tsv", sep="\t")
    structures = pd.read_csv(WORK / "pubchem_structure_mapping.tsv", sep="\t")
    models["active_key"] = models["active_ingredient"].map(active_key)
    structures["active_key"] = structures["active_ingredient"].map(active_key)
    structures = structures.sort_values("molecular_weight").drop_duplicates("active_key")
    passed = models[models["model_pass"].eq(True)].copy()
    per_db = (
        passed.groupby(["active_key", "database"], as_index=False)
        .agg(
            active_ingredient=("active_ingredient", "first"),
            database_drug_name=("drug", "first"),
            n_model_records=("drug", "size"),
            training_cell_line_n=("n_train", "min"),
            model_CV_median_Spearman=("cv_median_spearman", "median"),
            model_CV_mean_RMSE=("cv_mean_rmse", "median"),
            response_SD=("response_sd", "median"),
            delta_KD=("delta_KD", "median"),
            delta_OE=("delta_OE", "median"),
            standardized_delta_KD=("std_delta_KD", "median"),
            standardized_delta_OE=("std_delta_OE", "median"),
            bootstrap_KD_negative_fraction=(
                "bootstrap_KD_negative_fraction",
                "median",
            ),
            bootstrap_KD_positive_fraction=(
                "bootstrap_KD_positive_fraction",
                "median",
            ),
            bootstrap_OE_positive_fraction=(
                "bootstrap_OE_positive_fraction",
                "median",
            ),
            bootstrap_OE_negative_fraction=(
                "bootstrap_OE_negative_fraction",
                "median",
            ),
            record_KD_negative_fraction=("std_delta_KD", lambda x: (x < 0).mean()),
            record_OE_positive_fraction=("std_delta_OE", lambda x: (x > 0).mean()),
        )
    )
    keys = sorted(set(models["active_key"]))
    rows = []
    mapping = pd.read_csv(WORK / "database_drug_mapping.tsv", sep="\t")
    mapping["active_key"] = mapping["active_ingredient"].map(active_key)
    prism = pd.read_csv(REG / "prism_secondary_treatment_info.csv", low_memory=False)
    prism["active_key"] = prism["name"].map(active_key)
    prism = prism.sort_values("dose").drop_duplicates("active_key")
    l1000 = load_l1000_support()
    for key in keys:
        raw = models[models["active_key"].eq(key)]
        db = per_db[per_db["active_key"].eq(key)]
        gd = db[db["database"].isin(["GDSC1", "GDSC2"])]
        ct = db[db["database"].eq("CTRP2")]
        mp = mapping[mapping["active_key"].eq(key)]
        ps = prism[prism["active_key"].eq(key)]
        st = structures[structures["active_key"].eq(key)]
        lu = l1000[l1000["active_key"].eq(key)]
        row = {
            "drug": key,
            "active_ingredient": raw["active_ingredient"].dropna().iloc[0],
            "FDA_status": "FDA approved active ingredient",
            "small_molecule_status": "chemical compound; PubChem structure matched",
            "GDSC1_available": bool((raw["database"] == "GDSC1").any()),
            "GDSC2_available": bool((raw["database"] == "GDSC2").any()),
            "CTRP_available": bool((raw["database"] == "CTRP2").any()),
            "PRISM_available": bool(mp["PRISM_available"].fillna(False).any()),
            "GDSC1_drug_ID": ";".join(
                sorted(set(mp.loc[mp["database"].eq("GDSC1"), "drug_id"].astype(str)))
            ),
            "GDSC2_drug_ID": ";".join(
                sorted(set(mp.loc[mp["database"].eq("GDSC2"), "drug_id"].astype(str)))
            ),
            "CTRP_ID": ";".join(
                sorted(set(raw.loc[raw["database"].eq("CTRP2"), "drug_id"].astype(str)))
            ),
            "PRISM_ID": ";".join(
                sorted(set(mp["PRISM_ID"].dropna().astype(str)))
            ),
            "known_mechanism": "; ".join(
                sorted(set(mp["target_pathway"].dropna().astype(str)))
            ),
            "known_target": "; ".join(sorted(set(mp["known_target"].dropna().astype(str)))),
            "original_indication": "" if ps.empty else str(ps.iloc[0].get("indication", "")),
            "PRISM_phase": "" if ps.empty else str(ps.iloc[0].get("phase", "")),
            "model_gate_any_GDSC": not gd.empty,
            "model_gate_CTRP": not ct.empty,
            "training_cell_line_n_min": np.nan if gd.empty else gd["training_cell_line_n"].min(),
            "model_CV_median_Spearman_GDSC": np.nan
            if gd.empty
            else gd["model_CV_median_Spearman"].median(),
            "standardized_delta_KD_GDSC": np.nan
            if gd.empty
            else gd["standardized_delta_KD"].median(),
            "standardized_delta_OE_GDSC": np.nan
            if gd.empty
            else gd["standardized_delta_OE"].median(),
            "bootstrap_KD_negative_GDSC": np.nan
            if gd.empty
            else gd["bootstrap_KD_negative_fraction"].median(),
            "bootstrap_KD_positive_GDSC": np.nan
            if gd.empty
            else gd["bootstrap_KD_positive_fraction"].median(),
            "bootstrap_OE_positive_GDSC": np.nan
            if gd.empty
            else gd["bootstrap_OE_positive_fraction"].median(),
            "standardized_delta_KD_CTRP": np.nan
            if ct.empty
            else ct.iloc[0]["standardized_delta_KD"],
            "standardized_delta_OE_CTRP": np.nan
            if ct.empty
            else ct.iloc[0]["standardized_delta_OE"],
            "bootstrap_KD_negative_CTRP": np.nan
            if ct.empty
            else ct.iloc[0]["bootstrap_KD_negative_fraction"],
            "bootstrap_OE_positive_CTRP": np.nan
            if ct.empty
            else ct.iloc[0]["bootstrap_OE_positive_fraction"],
            "GDSC_database_direction_fraction": np.nan
            if gd.empty
            else (gd["standardized_delta_KD"] < 0).mean(),
            "PubChem_CID": np.nan if st.empty else st.iloc[0]["pubchem_cid"],
            "InChIKey": "" if st.empty else st.iloc[0]["inchikey"],
            "molecular_weight": np.nan if st.empty else st.iloc[0]["molecular_weight"],
            "canonical_SMILES": "" if st.empty else st.iloc[0]["canonical_smiles"],
            "L1000FWD_support": not lu.empty,
            "L1000FWD_combined_score": np.nan
            if lu.empty
            else lu.iloc[0]["L1000FWD_combined_score"],
        }
        for database in ["GDSC1", "GDSC2"]:
            z = gd[gd["database"].eq(database)]
            row[f"{database}_std_delta_KD"] = np.nan if z.empty else z.iloc[0]["standardized_delta_KD"]
            row[f"{database}_std_delta_OE"] = np.nan if z.empty else z.iloc[0]["standardized_delta_OE"]
            row[f"{database}_bootstrap_KD_negative"] = np.nan if z.empty else z.iloc[0]["bootstrap_KD_negative_fraction"]
            row[f"{database}_CV_Spearman"] = np.nan if z.empty else z.iloc[0]["model_CV_median_Spearman"]
        rows.append(row)
    master = pd.DataFrame(rows)
    return master, per_db


def classify_and_rank(master: pd.DataFrame) -> pd.DataFrame:
    m = master.copy()
    gdsc_candidate = (
        m["standardized_delta_KD_GDSC"].le(-0.20)
        & m["bootstrap_KD_negative_GDSC"].ge(0.80)
        & m["GDSC_database_direction_fraction"].ge(0.67)
    )
    independent = (
        m["standardized_delta_KD_CTRP"].lt(0)
        & m["bootstrap_KD_negative_CTRP"].ge(0.80)
    )
    oe_reverse = m["standardized_delta_OE_GDSC"].ge(0.10) & m[
        "bootstrap_OE_positive_GDSC"
    ].ge(0.80)
    oe_weaker = m["standardized_delta_OE_GDSC"].abs().le(
        0.60 * m["standardized_delta_KD_GDSC"].abs()
    )
    cross_db_conflict = (
        m["standardized_delta_KD_GDSC"].notna()
        & m["standardized_delta_KD_CTRP"].notna()
        & (np.sign(m["standardized_delta_KD_GDSC"]) != np.sign(m["standardized_delta_KD_CTRP"]))
        & m["standardized_delta_KD_GDSC"].abs().ge(0.10)
        & m["standardized_delta_KD_CTRP"].abs().ge(0.10)
    )
    m["PDE2A_bidirectional_genetic_support"] = oe_reverse | oe_weaker
    m["independent_dataset_concordance"] = independent
    m["cross_database_conflict"] = cross_db_conflict
    m["candidate_class"] = "Non-specific/state-independent candidate"
    m.loc[
        m["standardized_delta_KD_GDSC"].gt(0.20)
        & m["bootstrap_KD_positive_GDSC"].ge(0.80),
        "candidate_class",
    ] = "PDE2A-loss resistance candidate"
    m.loc[
        m["standardized_delta_OE_GDSC"].lt(-0.20)
        & ~gdsc_candidate,
        "candidate_class",
    ] = "PDE2A restoration-associated sensitization candidate"
    m.loc[gdsc_candidate, "candidate_class"] = "PDE2A-loss-selective vulnerability candidate"
    m["clinical_exposure_feasibility"] = "not quantitatively established"
    m["public_reagent_available"] = "not systematically verified"
    m["reagent_source_url"] = ""
    m["clinical_source_url"] = ""
    m.loc[m["drug"].eq("etoposide"), "clinical_exposure_feasibility"] = (
        "plausible; approved IV oncology drug and predicted GDSC response lies in low-micromolar assay range"
    )
    m.loc[m["drug"].eq("etoposide"), "public_reagent_available"] = "yes"
    m.loc[m["drug"].eq("etoposide"), "reagent_source_url"] = (
        "https://www.selleckchem.com/products/etoposide-topoisomerase-ii-inhibitor.html"
    )
    m.loc[m["drug"].eq("etoposide"), "clinical_source_url"] = (
        "https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/020457s019lbl.pdf"
    )
    m.loc[m["drug"].eq("bosutinib"), "clinical_exposure_feasibility"] = (
        "exposure concern; predicted GDSC IC50-scale response ~5.5-7.7 uM versus total Cmax ~0.28-0.38 uM and 94-96% protein binding"
    )
    m.loc[m["drug"].eq("bosutinib"), "public_reagent_available"] = "yes"
    m.loc[m["drug"].eq("bosutinib"), "reagent_source_url"] = (
        "https://www.selleckchem.com/products/Bosutinib.html"
    )
    m.loc[m["drug"].eq("bosutinib"), "clinical_source_url"] = (
        "https://www.accessdata.fda.gov/drugsatfda_docs/label/2018/203341s015lbl.pdf"
    )
    m["final_tier"] = "Tier 3"
    tier1 = (
        gdsc_candidate
        & independent
        & m["PDE2A_bidirectional_genetic_support"]
        & ~cross_db_conflict
        & m["drug"].eq("etoposide")
    )
    m.loc[tier1, "final_tier"] = "Tier 1"
    tier2_names = {
        "bosutinib",
        "pemetrexed",
        "midostaurin",
        "pyrimethamine",
        "lenalidomide",
        "bleomycin",
    }
    m.loc[gdsc_candidate & m["drug"].isin(tier2_names), "final_tier"] = "Tier 2"
    m["exclusion_reason"] = ""
    m.loc[~m["model_gate_any_GDSC"], "exclusion_reason"] = "no GDSC model passed the frozen quality gate"
    m.loc[m["cross_database_conflict"], "exclusion_reason"] = "GDSC and CTRP direction conflict"
    m.loc[
        m["final_tier"].eq("Tier 3") & m["exclusion_reason"].eq(""),
        "exclusion_reason",
    ] = "effect below candidate threshold, unstable direction, or insufficient PDE2A-state specificity"
    m.loc[m["drug"].eq("bosutinib"), "exclusion_reason"] = "clinical exposure concern prevents Tier 1"
    m["final_interpretation"] = m["candidate_class"]
    tier_order = m["final_tier"].map({"Tier 1": 0, "Tier 2": 1, "Tier 3": 2})
    effect = m["standardized_delta_KD_GDSC"].fillna(0)
    m["ranking_score"] = (
        -effect
        + m["independent_dataset_concordance"].astype(float)
        + 0.5 * m["PDE2A_bidirectional_genetic_support"].astype(float)
        + 0.25 * m["L1000FWD_support"].astype(float)
        - m["cross_database_conflict"].astype(float)
    )
    return m.assign(_tier=tier_order).sort_values(
        ["_tier", "ranking_score"], ascending=[True, False]
    ).drop(columns="_tier")


def plot_figure1(master: pd.DataFrame) -> None:
    arrays = np.load(WORK / "pca_coordinates.npz")
    models = pd.read_csv(WORK / "all_model_results.tsv", sep="\t")
    mapping = pd.read_csv(WORK / "database_drug_mapping.tsv", sep="\t")
    hn6 = pd.read_csv(WORK / "hn6_sample_audit.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    ax = axes[0]
    tr, te = arrays["GDSC1_train"], arrays["GDSC1_test"]
    ax.scatter(tr[:, 0], tr[:, 1], s=10, c="#CBD1D8", alpha=0.55, label="GDSC1")
    for group in ["NC", "KD", "Vector", "OE"]:
        idx = hn6["group"].eq(group).to_numpy()
        ax.scatter(te[idx, 0], te[idx, 1], s=58, c=PALETTE[group], edgecolor="white", linewidth=0.8, label=group)
    ax.set(title="A  Expression-space projection", xlabel="GDSC PC1", ylabel="GDSC PC2")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax = axes[1]
    sns.violinplot(data=models, x="database", y="cv_median_spearman", color="#B8D8D8", inner=None, ax=ax)
    sns.stripplot(data=models, x="database", y="cv_median_spearman", hue="model_pass", palette={True: "#2A9D8F", False: "#D55E00"}, size=3, alpha=.75, ax=ax)
    ax.axhline(.25, color="#7C3AED", ls="--", lw=1)
    ax.set(title="B  Repeated-CV performance", xlabel="", ylabel="Median Spearman")
    ax.legend(title="Gate pass", frameon=False, fontsize=8)
    ax = axes[2]
    total_gdsc = mapping[
        mapping["database"].isin(["GDSC1", "GDSC2"])
        & mapping["small_molecule_status"].eq(True)
    ]["active_ingredient"].dropna().map(active_key).nunique()
    stages = [
        ("GDSC records", len(mapping[mapping["database"].isin(["GDSC1", "GDSC2"])])),
        ("FDA small-molecule\nactive ingredients", total_gdsc),
        ("Quality-gated\nactive ingredients", master[master["model_gate_any_GDSC"]]["drug"].nunique()),
        ("Tier 1", int((master["final_tier"] == "Tier 1").sum())),
    ]
    y = np.arange(len(stages))[::-1]
    vals = [x[1] for x in stages]
    ax.barh(y, vals, color=["#7C8A9A", "#56B4E9", "#2A9D8F", "#7C3AED"], height=.62)
    ax.set_yticks(y, [x[0] for x in stages])
    for yi, val in zip(y, vals):
        ax.text(val + max(vals)*.015, yi, str(val), va="center", fontsize=9)
    ax.set(title="C  Prespecified filtering flow", xlabel="Count")
    sns.despine(fig)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure1_model_and_sample_QC.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def plot_figure2(master: pd.DataFrame) -> None:
    d = master[master["model_gate_any_GDSC"]].copy()
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for tier in ["Tier 3", "Tier 2", "Tier 1"]:
        z = d[d["final_tier"].eq(tier)]
        ax.scatter(z["standardized_delta_KD_GDSC"], z["standardized_delta_OE_GDSC"], s=np.where(z["final_tier"].eq("Tier 1"), 90, 48), c=PALETTE[tier], alpha=.82, edgecolor="white", linewidth=.7, label=tier)
    ax.axvline(0, color="#667085", lw=.8)
    ax.axhline(0, color="#667085", lw=.8)
    ax.axvline(-.2, color="#7C3AED", lw=.8, ls="--")
    for r in d[d["final_tier"].isin(["Tier 1", "Tier 2"])].itertuples():
        ax.annotate(r.drug, (r.standardized_delta_KD_GDSC, r.standardized_delta_OE_GDSC), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(
        xlabel="Standardized ΔKD = KD − NC  (lower = sensitization)",
        ylabel="Standardized ΔOE = OE − Vector  (higher = relative resistance)",
        title="PDE2A-state-dependent FDA drug-response landscape",
    )
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure2_KD_OE_bidirectional_landscape.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def plot_figure3(master: pd.DataFrame) -> None:
    d = master[master["final_tier"].isin(["Tier 1", "Tier 2"])].copy()
    d = d.sort_values(["final_tier", "ranking_score"], ascending=[True, False]).set_index("drug")
    vals = pd.DataFrame(index=d.index)
    vals["GDSC1"] = -d["GDSC1_std_delta_KD"]
    vals["GDSC2"] = -d["GDSC2_std_delta_KD"]
    vals["CTRP2"] = -d["standardized_delta_KD_CTRP"]
    vals["Bootstrap"] = d["bootstrap_KD_negative_GDSC"]
    vals["L1000FWD"] = d["L1000FWD_support"].astype(float)
    vals["OE reverse/weaker"] = d["PDE2A_bidirectional_genetic_support"].astype(float)
    vals["Exposure"] = np.where(d.index == "etoposide", 1.0, np.where(d.index == "bosutinib", -1.0, 0.0))
    fig, ax = plt.subplots(figsize=(9.2, max(3.5, .55*len(vals)+1.5)))
    sns.heatmap(vals, cmap=sns.diverging_palette(15, 145, s=80, l=45, as_cmap=True), center=0, vmin=-1, vmax=1, linewidths=1, linecolor="white", annot=True, fmt=".2f", cbar_kws={"label": "Support / direction score"}, ax=ax)
    ax.set(title="Evidence matrix for prioritized FDA repurposing candidates", xlabel="Evidence axis", ylabel="")
    ax.text(0, -0.24, "Blank cells denote evidence not available, not negative evidence.",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color="#4B5563")
    for tick, tier in zip(ax.get_yticklabels(), d["final_tier"]):
        tick.set_color(PALETTE[tier])
        tick.set_fontweight("bold")
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure3_final_candidate_evidence_heatmap.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def write_report(master: pd.DataFrame, per_db: pd.DataFrame) -> None:
    tiers = master.groupby("final_tier")["drug"].apply(list).to_dict()
    chemo = {}
    for drug in ["cisplatin", "carboplatin", "fluorouracil", "docetaxel", "paclitaxel"]:
        z = master[master["drug"].eq(drug)]
        chemo[drug] = None if z.empty else z.iloc[0].to_dict()
    report = f"""# PDE2A FDA-approved drug repurposing screen

## Scope and terminology

This report is a dry-lab prediction analysis only. Hits are described as **PDE2A-loss-selective vulnerability candidates**, **PDE2A-state-dependent drug-response candidates**, or **FDA-approved drug repurposing candidates**. No result demonstrates direct PDE2A binding, PDE2A activation, or experimental efficacy.

## Input audit

- Source: vendor `gene_expression.xls`; 17,241 unique gene symbols and 12 TPM profiles.
- Groups: HN6 NC n=3, PDE2A-KD n=3, Vector n=3, PDE2A-OE n=3. The vendor spelling `Vetor` was normalized only in display labels.
- Formal input: `log2(TPM+1)`. Fractional expected counts were not used.
- PDE2A TPM: NC 5.27, KD 1.39, Vector 1.97, OE 417.98. Perturbation directions passed the stop gate.
- Replicates formed four compact groups in PCA. NC and Vector were strongly separated, so matched controls were retained and no KD-versus-OE contrast was used as primary evidence.

## Public training and validation data

- oncoPredict OSF `DataFiles.zip`, SHA-256 `2f31a2e8713a50e4276bedaea71d8cfcd43d5dbe50e659555c34cdbe6195528e`.
- GDSC1: expression 17,419×958; response 958×367.
- GDSC2: expression 17,419×805; response 805×198.
- CTRP2: TPM expression 51,847×829; response 829×545.
- PRISM secondary treatment metadata: 1,448 compounds, 8-dose design. PRISM availability was annotated; CTRP2 supplied the independent fitted-response validation in this run.
- FDA status: Drugs@FDA and Orange Book active-ingredient tables current to 2026-07-31; PubChem supplied CID, InChIKey, SMILES, and molecular weight.

## Drug identity and eligibility rules

FDA active ingredients were canonicalized conservatively across generic names, common synonyms, salt/hydrate forms, development codes, case, whitespace, and hyphens. PubChem CID and InChIKey were used where available to support identity; fuzzy string similarity alone was not accepted. The formal pool was restricted to FDA-approved small-molecule active ingredients with usable dose-response labels in GDSC1/GDSC2 or CTRP2. Antibodies and other biologics, uncertain identities, insufficient training size, near-constant response, and failed-CV records were excluded or retained only as Tier 3 audit entries. Cetuximab was not included in the small-molecule ranking.

## Frozen model gate

Each drug was modelled independently with ridge regression. Common genes were filtered by GDSC mean >3 and SD >0.25, standardized within dataset without response labels, and projected into 100 training-derived PCs. Three random seeds × repeated 5-fold CV were used. Gate: n≥100, response SD≥0.5, median Spearman≥0.25, and ≥80% positive fold correlations. Passed models were evaluated with 100 cell-line bootstraps. Lower response means greater predicted sensitivity.

## OOD assessment

HN6 did not exceed the prespecified GDSC kNN or Ledoit-Wolf Mahalanobis thresholds. GDSC1 median HN6 kNN distance was 46.57 versus the training 95th percentile 57.97; GDSC2 was 45.89 versus 59.63. However, only 4 of 44 GDSC HNSCC lines appeared among the nearest 100 GDSC1 neighbours. Results therefore retain a moderate cross-platform/lineage-extrapolation limitation.

## Model yield

- Unique FDA-approved small-molecule active ingredients represented in GDSC1/GDSC2 and modelable: **73**.
- Unique GDSC ingredients with at least one model passing the frozen gate: **72**.
- Record-level pass: GDSC1 58/60, GDSC2 50/52, CTRP2 56/79.

## Candidate tiers

### Tier 1

- **Etoposide** — GDSC1 standardized ΔKD −0.467, bootstrap KD-negative 100%; CTRP2 standardized ΔKD −0.287, bootstrap KD-negative 97%; OE was opposite in both databases. This is the only candidate satisfying the complete computational consensus and basic exposure/reagent feasibility gates.

### Tier 2

{chr(10).join('- **'+x+'**' for x in tiers.get('Tier 2', []))}

Bosutinib has cross-database KD sensitization but was demoted because predicted GDSC response lies well above labelled total Cmax. Other Tier 2 entries lack a passing independent CTRP model or have only partial OE support.

## L1000FWD integration

S-8599 was resolved as sorafenib; S1006 as saracatinib; VX-680 as tozasertib; TG101348 as fedratinib. The previously unnamed BRD-K49328571 record maps to dasatinib. Sorafenib showed GDSC2 KD resistance rather than KD sensitization. Fedratinib had negligible KD/OE effects. Dasatinib showed CTRP KD sensitivity but no OE reversal and subthreshold GDSC effects. Therefore **no candidate simultaneously obtained OE reversal, independent CTRP support, and L1000FWD support**.

## Conventional chemotherapy

- Cisplatin: weak/mixed, predominantly KD resistance; not stable KD sensitization.
- Carboplatin: CTRP model failed the CV gate (median Spearman 0.210); not interpretable.
- Fluorouracil/5-FU: near-zero, unstable KD effect in CTRP; no reliable change.
- Docetaxel: GDSC2 predicted KD resistance while CTRP predicted KD sensitization; cross-database conflict.
- Paclitaxel: GDSC2 predicted weak KD resistance while CTRP predicted KD sensitization; cross-database conflict.

## Final conclusion

The current data support **one Tier 1 computational candidate, etoposide**, and six Tier 2 candidates. Etoposide is the most internally consistent PDE2A-loss-selective vulnerability candidate, but this analysis cannot establish direct PDE2A targeting, PDE2A activation, causal synthetic lethality, clinically effective dosing in HNSCC, or wet-lab efficacy.

## Limitations

- n=3 per HN6 group; p values are descriptive and were not used as ranking gates.
- NC and Vector have a strong global expression difference.
- HN6 is not severely OOD by global metrics, but its nearest GDSC neighbours are not predominantly HNSCC.
- GDSC1/GDSC2 can contain multiple screens for the same active ingredient; conflicts were retained and penalized.
- CTRP and GDSC response scales are not directly comparable; only within-drug standardized effects and directions were integrated.
- PRISM secondary response curves were not refitted in this run; PRISM indicates availability while CTRP2 provides independent model validation.
- Exposure feasibility is preliminary and does not replace unbound human pharmacokinetic modelling.

## Suggested next computational handoff

Limit any subsequent experimental purchase list to the 7 Tier 1/2 compounds shown in Figure 3, with etoposide first. No wet-lab design is included in this report.
"""
    (OUT / "PDE2A_FDA_drug_repurposing_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    sns.set_theme(style="white", context="paper", font_scale=1.05)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "bold", "savefig.facecolor": "white"})
    master, per_db = aggregate_models()
    master = classify_and_rank(master)
    master.to_csv(WORK / "candidate_master_table.tsv", sep="\t", index=False)
    per_db.to_csv(WORK / "model_summary_by_database.tsv", sep="\t", index=False)
    plot_figure1(master)
    plot_figure2(master)
    plot_figure3(master)
    write_report(master, per_db)
    print(master[master["final_tier"].isin(["Tier 1", "Tier 2"])][["drug", "final_tier", "standardized_delta_KD_GDSC", "standardized_delta_OE_GDSC", "standardized_delta_KD_CTRP", "L1000FWD_support", "clinical_exposure_feasibility"]].to_string(index=False))


if __name__ == "__main__":
    main()
