#!/usr/bin/env python3
"""Checkpoint 8: STRING-only PPI context for the frozen 84 FDR-significant DEGs.

No external nodes are added. The primary network uses score >= 0.400 and the
score >= 0.700 network is sensitivity-only. Centrality values are descriptive;
this script never computes a composite candidate score.
"""

from __future__ import annotations

import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests

SEED = 20260829
SPECIES = 9606
BASE = "https://string-db.org/api/json"
USER_AGENT = "PDE2A-bulkRNAseq-checkpoint8/1.0"


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty frozen input: {path}")


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="")


def api_post(method: str, data: dict, raw_path: Path) -> list[dict]:
    """Call STRING once, or reuse the frozen raw JSON response."""
    if raw_path.is_file() and raw_path.stat().st_size > 0:
        return json.loads(raw_path.read_text(encoding="utf-8"))
    response = requests.post(
        f"{BASE}/{method}", data=data, timeout=180,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()
    raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


def api_get(method: str, raw_path: Path) -> list[dict]:
    if raw_path.is_file() and raw_path.stat().st_size > 0:
        return json.loads(raw_path.read_text(encoding="utf-8"))
    response = requests.get(
        f"{BASE}/{method}", timeout=60,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()
    raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


def contains_class(value: str, label: str) -> bool:
    return label in {x.strip() for x in str(value).split(";") if x.strip()}


def rank_desc(series: pd.Series) -> pd.Series:
    return series.rank(method="min", ascending=False).astype("Int64")


def top_counts(values: list[str], limit: int = 8) -> str:
    counts = Counter(v for v in values if v)
    return ";".join(f"{k}={v}" for k, v in counts.most_common(limit))


def component_metrics(graph: nx.Graph) -> tuple[pd.DataFrame, dict[str, int], list[set[str]]]:
    components = sorted(nx.connected_components(graph), key=lambda x: (-len(x), sorted(x)[0]))
    component_id = {gene: i + 1 for i, genes in enumerate(components) for gene in genes}
    degree = dict(graph.degree())
    strength = dict(graph.degree(weight="score"))
    betweenness = nx.betweenness_centrality(graph, weight="distance", normalized=True)
    closeness: dict[str, float] = {}
    eigenvector: dict[str, float] = {}
    for genes in components:
        sub = graph.subgraph(genes).copy()
        if len(genes) == 1:
            gene = next(iter(genes))
            closeness[gene] = math.nan
            eigenvector[gene] = math.nan
            continue
        closeness.update(nx.closeness_centrality(sub, distance="distance"))
        try:
            eigenvector.update(nx.eigenvector_centrality(sub, weight="score", max_iter=5000, tol=1e-10))
        except nx.PowerIterationFailedConvergence:
            for gene in genes:
                eigenvector[gene] = math.nan
    metrics = pd.DataFrame({
        "gene": list(graph.nodes()),
        "component_id": [component_id[x] for x in graph.nodes()],
        "component_size": [len(components[component_id[x] - 1]) for x in graph.nodes()],
        "degree": [degree[x] for x in graph.nodes()],
        "weighted_degree_strength": [strength[x] for x in graph.nodes()],
        "betweenness": [betweenness[x] for x in graph.nodes()],
        "closeness_within_component": [closeness[x] for x in graph.nodes()],
        "eigenvector_within_component": [eigenvector[x] for x in graph.nodes()],
    })
    return metrics, component_id, components


def make_graph(genes: list[str], edges: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    for row in edges.itertuples(index=False):
        graph.add_edge(row.gene_A, row.gene_B, score=float(row.combined_score), distance=1.0 / float(row.combined_score))
    return graph


def network_source(nodes: pd.DataFrame, edges: pd.DataFrame, positions: dict[str, np.ndarray]) -> pd.DataFrame:
    node_rows = nodes.copy()
    node_rows["record_type"] = "node"
    node_rows["x"] = node_rows["gene"].map(lambda x: positions[x][0])
    node_rows["y"] = node_rows["gene"].map(lambda x: positions[x][1])
    node_rows["xend"] = np.nan
    node_rows["yend"] = np.nan
    node_rows["gene_A"] = ""
    node_rows["gene_B"] = ""
    node_rows["combined_score"] = np.nan
    edge_rows = pd.DataFrame({
        "record_type": "edge",
        "gene": "",
        "x": edges["gene_A"].map(lambda x: positions[x][0]),
        "y": edges["gene_A"].map(lambda x: positions[x][1]),
        "xend": edges["gene_B"].map(lambda x: positions[x][0]),
        "yend": edges["gene_B"].map(lambda x: positions[x][1]),
        "gene_A": edges["gene_A"],
        "gene_B": edges["gene_B"],
        "combined_score": edges["combined_score"],
    })
    for column in node_rows.columns:
        if column not in edge_rows.columns:
            edge_rows[column] = np.nan
    return pd.concat([edge_rows[node_rows.columns], node_rows], ignore_index=True)


def packed_component_layout(graph: nx.Graph) -> dict[str, np.ndarray]:
    """Keep the largest component legible and place smaller components outside it."""
    components = sorted(nx.connected_components(graph), key=lambda x: (-len(x), sorted(x)[0]))
    positions: dict[str, np.ndarray] = {}
    largest = graph.subgraph(components[0]).copy()
    positions.update(nx.spring_layout(largest, seed=SEED, weight="score", iterations=700, scale=1.25))
    peripheral = components[1:]
    radius = 2.45
    for index, geneset in enumerate(peripheral):
        angle = 2 * math.pi * index / max(1, len(peripheral))
        center = np.array([radius * math.cos(angle), radius * math.sin(angle)])
        sub = graph.subgraph(geneset).copy()
        if len(geneset) == 1:
            local = {next(iter(geneset)): np.array([0.0, 0.0])}
        else:
            local = nx.spring_layout(sub, seed=SEED + index + 1, weight="score", iterations=300, scale=0.18)
        for gene, coordinate in local.items():
            positions[gene] = center + coordinate
    return positions


def main() -> None:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd().resolve()
    conv_path = root / "results/07_candidate_convergence/CANDIDATE_CONVERGENCE_TABLE.tsv"
    pde_genome_path = root / "results/05_deseq2/PDE2A_GENOMEWIDE_POSITION.tsv"
    pathway_map_path = root / "results/07_candidate_convergence/PATHWAY_TO_FUNCTIONAL_MODULE.tsv"
    for path in (conv_path, pde_genome_path, pathway_map_path):
        require_file(path)

    out = root / "results/08_network_prioritization"
    raw = out / "raw_STRING_responses"
    fig_root = root / "results/07_figures/manuscript_candidates/Network_prioritization"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    fig_root.mkdir(parents=True, exist_ok=True)

    conv = pd.read_csv(conv_path, sep="\t")
    genome_pde = pd.read_csv(pde_genome_path, sep="\t")
    if len(conv) != 84 or conv["gene_name"].duplicated().any():
        raise ValueError("Frozen Checkpoint 7 candidate table must contain 84 unique genes")
    if len(genome_pde) != 1 or genome_pde.iloc[0]["gene_name"] != "PDE2A":
        raise ValueError("Frozen genome-wide PDE2A rank row is invalid")
    conv["signaling_regulatory_module_member"] = conv["functional_modules"].fillna("").map(
        lambda x: "Signaling / regulatory" in {y.strip() for y in x.split(";")}
    )
    conv["signaling_regulatory_functional_class"] = conv["functional_class"].fillna("").map(
        lambda x: contains_class(x, "signaling/regulatory")
    )
    module_count = int(conv["signaling_regulatory_module_member"].sum())
    functional_edge_down_count = int(((conv["direction"] == "Edge-down") & conv["signaling_regulatory_functional_class"]).sum())
    if module_count != 17 or functional_edge_down_count != 35:
        raise ValueError(f"Checkpoint 7 definition counts changed: module={module_count}, class_edge_down={functional_edge_down_count}")

    pde84 = conv.loc[conv["gene_name"] == "PDE2A"].iloc[0]
    pg = genome_pde.iloc[0]
    reconciliation = pd.DataFrame([
        ["signaling_regulatory_module_member", "DEG belongs to >=1 significant GSEA pathway mapped to Signaling / regulatory", module_count, "84 FDR-significant DEGs", 84, "DEFINITION_CONFIRMED"],
        ["signaling_regulatory_functional_class", "Edge-down DEG has signaling/regulatory in frozen multi-label functional_class", functional_edge_down_count, "66 Edge-down FDR-significant DEGs", 66, "DEFINITION_CONFIRMED"],
        ["PDE2A_padj_rank_genomewide", "padj ascending", int(pg["rank_padj_ascending"]), "20,413 tested genes", 20413, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_abs_Wald_rank_genomewide", "absolute Wald statistic descending", int(pg["rank_abs_Wald_stat_descending"]), "20,413 tested genes", 20413, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_abs_shrunken_LFC_rank_genomewide", "absolute shrunken LFC descending", int(pg["rank_abs_shrunken_LFC_descending"]), "20,413 tested genes", 20413, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_padj_rank_FDR84", "padj ascending", int(pde84["padj_rank"]), "84 FDR-significant DEGs", 84, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_abs_Wald_rank_FDR84", "absolute Wald statistic descending", int(pde84["absolute_Wald_rank"]), "84 FDR-significant DEGs", 84, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_abs_shrunken_LFC_rank_FDR84", "absolute shrunken LFC descending", int(pde84["absolute_shrunken_LFC_rank"]), "84 FDR-significant DEGs", 84, "DENOMINATOR_CONFIRMED"],
        ["PDE2A_leading_edge_count_rank_FDR84", "significant leading-edge count descending", int(pde84["leading_edge_count_rank"]), "84 FDR-significant DEGs", 84, "DENOMINATOR_CONFIRMED"],
    ], columns=["field_name", "definition", "value", "rank_or_count_universe", "denominator", "audit_status"])
    write_tsv(reconciliation, out / "CHECKPOINT7_DEFINITION_RECONCILIATION.tsv")

    version_payload = api_get("version", raw / "STRING_version.json")
    string_version = str(version_payload[0].get("string_version", "unknown")) if version_payload else "unknown"
    genes = conv["gene_name"].tolist()
    identifiers = "\r".join(genes)
    mapping_payload = api_post("get_string_ids", {"identifiers": identifiers, "species": SPECIES, "limit": 5}, raw / "STRING_mapping.json")
    by_query: dict[int, list[dict]] = defaultdict(list)
    for record in mapping_payload:
        by_query[int(record["queryIndex"])].append(record)
    string_id_inputs: dict[str, list[str]] = defaultdict(list)
    for idx, records in by_query.items():
        for record in records:
            string_id_inputs[record["stringId"]].append(genes[idx])
    mapping_rows = []
    input_to_string: dict[str, str] = {}
    string_to_input: dict[str, str] = {}
    for idx, gene in enumerate(genes):
        records = by_query.get(idx, [])
        if not records:
            mapping_rows.append([gene, "", "", "UNMAPPED", "NOT_APPLICABLE", 0, "No STRING mapping returned"])
            continue
        record = records[0]
        sid = record["stringId"]
        preferred = record.get("preferredName", "")
        ambiguity = "UNAMBIGUOUS"
        if len(records) > 1:
            ambiguity = "MULTIPLE_RECORDS_RETURNED"
        if len(set(string_id_inputs[sid])) > 1:
            ambiguity = "DUPLICATE_STRING_ID_FOR_MULTIPLE_INPUTS"
        status = "MAPPED_EXACT" if preferred.upper() == gene.upper() else "MAPPED_ALIAS"
        input_to_string[gene] = sid
        string_to_input[sid] = gene
        mapping_rows.append([gene, sid, preferred, status, ambiguity, len(records), record.get("annotation", "")])
    mapping = pd.DataFrame(mapping_rows, columns=[
        "input_gene", "STRING_mapped_ID", "STRING_preferred_name", "mapping_status",
        "mapping_ambiguity", "mapping_record_count", "STRING_annotation",
    ])
    write_tsv(mapping, out / "PPI_STRING_MAPPING_AUDIT.tsv")
    pde_map = mapping.loc[mapping["input_gene"] == "PDE2A"].iloc[0]
    if not str(pde_map["mapping_status"]).startswith("MAPPED") or pde_map["mapping_ambiguity"] != "UNAMBIGUOUS":
        raise RuntimeError("PDE2A did not map unambiguously to STRING; stopping before network analysis")

    mapped_genes = mapping.loc[mapping["mapping_status"].str.startswith("MAPPED"), "input_gene"].tolist()
    mapped_identifiers = "\r".join(mapped_genes)
    network_payloads = {}
    for threshold in (400, 700):
        network_payloads[threshold] = api_post(
            "network",
            {"identifiers": mapped_identifiers, "species": SPECIES, "required_score": threshold,
             "network_type": "functional", "add_nodes": 0},
            raw / f"STRING_network_score{threshold}.json",
        )

    def clean_edges(payload: list[dict], threshold: int) -> pd.DataFrame:
        rows = []
        for record in payload:
            sid_a, sid_b = record["stringId_A"], record["stringId_B"]
            if sid_a not in string_to_input or sid_b not in string_to_input:
                raise RuntimeError("STRING returned an external node despite add_nodes=0")
            a, b = string_to_input[sid_a], string_to_input[sid_b]
            if a == b:
                continue
            a, b = sorted((a, b))
            score = float(record["score"])
            if score * 1000 + 1e-9 < threshold:
                continue
            rows.append({
                "gene_A": a, "gene_B": b, "STRING_ID_A": input_to_string[a], "STRING_ID_B": input_to_string[b],
                "combined_score": score, "nscore": record.get("nscore", np.nan),
                "fscore": record.get("fscore", np.nan), "pscore": record.get("pscore", np.nan),
                "ascore": record.get("ascore", np.nan), "escore": record.get("escore", np.nan),
                "dscore": record.get("dscore", np.nan), "tscore": record.get("tscore", np.nan),
            })
        if not rows:
            return pd.DataFrame(columns=["gene_A", "gene_B", "STRING_ID_A", "STRING_ID_B", "combined_score",
                                         "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"])
        frame = pd.DataFrame(rows).sort_values("combined_score", ascending=False)
        return frame.drop_duplicates(["gene_A", "gene_B"], keep="first").sort_values(["gene_A", "gene_B"]).reset_index(drop=True)

    edges400 = clean_edges(network_payloads[400], 400)
    edges700 = clean_edges(network_payloads[700], 700)
    write_tsv(edges400, out / "PPI_84DEG_STRING_edges_score400.tsv")
    write_tsv(edges700, out / "PPI_84DEG_STRING_edges_score700.tsv")
    graph400 = make_graph(genes, edges400)
    graph700 = make_graph(genes, edges700)
    metrics, component_id, components = component_metrics(graph400)

    node_columns = ["gene_name", "direction", "padj", "shrunken_log2FC", "Wald_stat", "baseMean",
                    "patient_direction_consistency", "total_significant_leading_edge_count",
                    "functional_module_count", "functional_modules", "functional_class",
                    "signaling_regulatory_module_member", "signaling_regulatory_functional_class"]
    nodes = conv[node_columns].rename(columns={
        "gene_name": "gene", "total_significant_leading_edge_count": "leading_edge_count",
        "patient_direction_consistency": "patient_consistency",
    }).merge(mapping[["input_gene", "STRING_mapped_ID", "mapping_status", "mapping_ambiguity"]],
             left_on="gene", right_on="input_gene", how="left").drop(columns="input_gene")
    nodes = nodes.merge(metrics, on="gene", how="left")
    write_tsv(nodes, out / "PPI_84DEG_STRING_nodes.tsv")
    write_tsv(metrics, out / "PPI_84DEG_NETWORK_METRICS.tsv")

    # Fixed weighted Louvain on the largest component only, if sufficiently large.
    community_id: dict[str, int | None] = {g: None for g in genes}
    communities: list[set[str]] = []
    largest = components[0]
    if len(largest) >= 10:
        lcc = graph400.subgraph(largest).copy()
        communities = list(nx.community.louvain_communities(lcc, weight="score", seed=SEED))
        communities = sorted(communities, key=lambda x: (-len(x), sorted(x)[0]))
        for idx, geneset in enumerate(communities, start=1):
            for gene in geneset:
                community_id[gene] = idx
    community_table = pd.DataFrame({
        "gene": genes,
        "community_id": [community_id[g] for g in genes],
        "community_size": [len(communities[community_id[g] - 1]) if community_id[g] is not None else np.nan for g in genes],
        "community_algorithm": ["weighted_Louvain_on_largest_component" if community_id[g] is not None else "not_in_largest_component" for g in genes],
    })
    write_tsv(community_table, out / "PPI_NETWORK_COMMUNITIES.tsv")
    nodes = nodes.merge(community_table[["gene", "community_id", "community_size"]], on="gene", how="left")

    context_rows = []
    for cid, geneset in enumerate(communities, start=1):
        sub = nodes[nodes["gene"].isin(geneset)]
        modules = [m.strip() for value in sub["functional_modules"].fillna("") for m in value.split(";") if m.strip()]
        classes = [m.strip() for value in sub["functional_class"].fillna("") for m in value.split(";") if m.strip()]
        context_rows.append({
            "community_id": cid, "community_size": len(geneset), "member_genes": ";".join(sorted(geneset)),
            "major_functional_modules": top_counts(modules), "major_functional_classes": top_counts(classes),
            "edge_down_n": int((sub["direction"] == "Edge-down").sum()),
            "edge_up_n": int((sub["direction"] == "Edge-up").sum()),
        })
    community_context = pd.DataFrame(context_rows)
    write_tsv(community_context, out / "PPI_COMMUNITY_FUNCTIONAL_CONTEXT.tsv")

    # PDE2A context is generated after all global topology calculations.
    pde_node = nodes.loc[nodes["gene"] == "PDE2A"].iloc[0]
    neighbors = sorted(graph400.neighbors("PDE2A"))
    context_records = [{
        "record_type": "summary", "gene": "PDE2A", "STRING_score": np.nan,
        "direction": pde_node["direction"], "shrunken_LFC": pde_node["shrunken_log2FC"],
        "padj": pde_node["padj"], "leading_edge_count": pde_node["leading_edge_count"],
        "functional_modules": pde_node["functional_modules"], "component_id": pde_node["component_id"],
        "community_id": pde_node["community_id"], "degree": pde_node["degree"],
        "weighted_degree_strength": pde_node["weighted_degree_strength"],
        "betweenness": pde_node["betweenness"], "closeness": pde_node["closeness_within_component"],
        "eigenvector_centrality": pde_node["eigenvector_within_component"],
        "retained_with_neighbor_at_score700": graph700.degree("PDE2A") > 0,
    }]
    for neighbor in neighbors:
        n = nodes.loc[nodes["gene"] == neighbor].iloc[0]
        context_records.append({
            "record_type": "direct_neighbor", "gene": neighbor,
            "STRING_score": graph400["PDE2A"][neighbor]["score"], "direction": n["direction"],
            "shrunken_LFC": n["shrunken_log2FC"], "padj": n["padj"],
            "leading_edge_count": n["leading_edge_count"], "functional_modules": n["functional_modules"],
            "component_id": n["component_id"], "community_id": n["community_id"],
            "degree": n["degree"], "weighted_degree_strength": n["weighted_degree_strength"],
            "betweenness": n["betweenness"], "closeness": n["closeness_within_component"],
            "eigenvector_centrality": n["eigenvector_within_component"],
            "retained_with_neighbor_at_score700": graph700.has_edge("PDE2A", neighbor),
        })
    pde_context = pd.DataFrame(context_records)
    write_tsv(pde_context, out / "PDE2A_PPI_CONTEXT.tsv")

    # Frozen functional-class selection only; no extra gates.
    regulatory = nodes[(nodes["direction"] == "Edge-down") & nodes["signaling_regulatory_functional_class"]].copy()
    regulatory = regulatory.merge(community_table[["gene", "community_id"]], on="gene", how="left", suffixes=("", "_community"))
    regulatory["enzyme_yes_no"] = regulatory["functional_class"].map(lambda x: contains_class(x, "enzyme"))
    regulatory["receptor_yes_no"] = regulatory["functional_class"].map(lambda x: contains_class(x, "receptor"))
    regulatory["structural_marker_yes_no"] = regulatory["functional_class"].map(lambda x: contains_class(x, "structural/marker"))
    reg_out = regulatory[[
        "gene", "padj", "shrunken_log2FC", "Wald_stat", "patient_consistency", "leading_edge_count",
        "functional_module_count", "functional_modules", "functional_class", "enzyme_yes_no",
        "receptor_yes_no", "structural_marker_yes_no", "degree", "betweenness", "community_id",
        "signaling_regulatory_module_member", "signaling_regulatory_functional_class",
    ]].rename(columns={"degree": "PPI_degree", "betweenness": "PPI_betweenness", "community_id": "PPI_community"})
    write_tsv(reg_out, out / "EDGE_DOWN_REGULATORY_CANDIDATES.tsv")

    evidence = reg_out.copy()
    n_reg = len(evidence)
    evidence["FDR_significant"] = True
    evidence["enzyme"] = evidence["enzyme_yes_no"]
    evidence["signaling_regulatory"] = evidence["signaling_regulatory_functional_class"]
    evidence["PPI_component"] = evidence["gene"].map(component_id)
    evidence["padj_rank"] = evidence["padj"].rank(method="min", ascending=True).astype("Int64")
    evidence["absolute_shrunken_LFC_rank"] = rank_desc(evidence["shrunken_log2FC"].abs())
    evidence["leading_edge_count_rank"] = rank_desc(evidence["leading_edge_count"])
    evidence["PPI_degree_rank"] = rank_desc(evidence["PPI_degree"])
    evidence["rank_universe"] = f"Edge-down signaling/regulatory functional-class candidates (n={n_reg})"
    evidence["rank_denominator"] = n_reg
    write_tsv(evidence, out / "CANDIDATE_EVIDENCE_MATRIX.tsv")

    pde_le = int(pde84["total_significant_leading_edge_count"])
    required_candidates = {"ADCY4", "KDR", "CAMK2B", "GNG11", "PDE1A", "PDE2A", "NR4A1"}
    auto_candidates = set(conv.loc[(conv["direction"] == "Edge-down") & (conv["total_significant_leading_edge_count"] >= pde_le), "gene_name"])
    compare_genes = sorted((required_candidates | auto_candidates) & set(genes))
    high_conv = nodes[nodes["gene"].isin(compare_genes)][[
        "gene", "padj", "shrunken_log2FC", "patient_consistency", "leading_edge_count",
        "functional_module_count", "functional_class", "degree", "community_id",
    ]].rename(columns={"degree": "PPI_degree", "community_id": "PPI_community"})
    high_conv["in_prespecified_observed_set"] = high_conv["gene"].isin(required_candidates)
    high_conv["meets_automatic_rule"] = high_conv["gene"].isin(auto_candidates)
    high_conv["automatic_rule"] = f"Edge-down FDR DEG with leading_edge_count >= PDE2A ({pde_le})"
    high_conv = high_conv.sort_values(["leading_edge_count", "padj", "gene"], ascending=[False, True, True])
    write_tsv(high_conv, out / "HIGH_CONVERGENCE_EDGE_DOWN_COMPARISON.tsv")

    # Fixed layouts and figure source data.
    positions = packed_component_layout(graph400)
    global_source = network_source(nodes, edges400, positions)
    layout_path = out / "PPI_84DEG_FIXED_LAYOUT.tsv"
    write_tsv(pd.DataFrame({"gene": genes, "x": [positions[g][0] for g in genes], "y": [positions[g][1] for g in genes]}), layout_path)

    def write_figure_source(name: str, data: pd.DataFrame) -> None:
        directory = fig_root / name
        directory.mkdir(parents=True, exist_ok=True)
        write_tsv(data, directory / "source_data.tsv")
        plot_script = root / "scripts/24_plot_checkpoint8_network.R"
        require_file(plot_script)
        (directory / "plotting_script.R").write_text(plot_script.read_text(encoding="utf-8"), encoding="utf-8")

    write_figure_source("Figure_PPI_01_84DEG_Global", global_source)
    write_figure_source("Figure_PPI_02_PDE2A_Context", global_source)
    if neighbors:
        local_genes = ["PDE2A"] + neighbors
        local_graph = graph400.subgraph(local_genes).copy()
        local_positions = nx.spring_layout(local_graph, seed=SEED, weight="score", iterations=500)
        local_nodes = nodes[nodes["gene"].isin(local_genes)].copy()
        local_edges = edges400[(edges400["gene_A"].isin(local_genes)) & (edges400["gene_B"].isin(local_genes))].copy()
        write_figure_source("Figure_PPI_03_PDE2A_LocalNeighborhood", network_source(local_nodes, local_edges, local_positions))

    matrix_fields = ["shrunken_log2FC", "patient_consistency", "leading_edge_count", "functional_module_count", "PPI_degree", "enzyme_yes_no"]
    ordered = evidence.sort_values(["leading_edge_count", "padj", "gene"], ascending=[False, True, True]).reset_index(drop=True)
    matrix_rows = []
    for order_idx, row in ordered.iterrows():
        for field in matrix_fields:
            raw_value = float(row[field]) if field != "enzyme_yes_no" else float(bool(row[field]))
            values = ordered[field].astype(float) if field != "enzyme_yes_no" else ordered[field].astype(int).astype(float)
            spread = values.max() - values.min()
            display = 0.5 if spread == 0 else (raw_value - values.min()) / spread
            matrix_rows.append({
                "gene": row["gene"], "evidence_field": field, "raw_value": raw_value,
                "within_field_display_scale": display, "display_order": order_idx + 1,
                "sorting_rule": "leading_edge_count descending; padj ascending; gene ascending",
                "note": "Within-column display scaling only; no values are summed or weighted",
            })
    write_figure_source("Figure_Prioritization_01_EvidenceMatrix", pd.DataFrame(matrix_rows))

    # Run copied R scripts after all source data are frozen.
    for directory in sorted(fig_root.iterdir()):
        if directory.is_dir() and (directory / "plotting_script.R").is_file():
            import subprocess
            subprocess.run(["Rscript", str(directory / "plotting_script.R"), str(directory)], check=True)

    components700 = sorted(nx.connected_components(graph700), key=lambda x: -len(x))
    mapped_n = int(mapping["mapping_status"].str.startswith("MAPPED").sum())
    unmapped_n = int((mapping["mapping_status"] == "UNMAPPED").sum())
    ambiguous_n = int((mapping["mapping_ambiguity"] != "UNAMBIGUOUS").sum() - unmapped_n)
    isolates400 = sorted(nx.isolates(graph400))
    isolates700 = sorted(nx.isolates(graph700))
    top_hubs = metrics.sort_values(["degree", "weighted_degree_strength", "gene"], ascending=[False, False, True]).head(10)
    pde_degree_rank = int(metrics["degree"].rank(method="min", ascending=False)[metrics["gene"] == "PDE2A"].iloc[0])
    pde_is_hub = bool(pde_degree_rank <= 10)
    pde_community_context = community_context.loc[community_context["community_id"] == pde_node["community_id"]]
    pde_comm_text = "not_in_largest_component" if pde_community_context.empty else pde_community_context.iloc[0]["major_functional_modules"]
    audit = [
        "# Checkpoint 8 PPI and transparent prioritization audit", "", f"Date: {pd.Timestamp.now().date()}", "",
        "## Scope and source", "",
        f"- STRING version: {string_version}; database date: not returned by the STRING version endpoint.",
        "- STRING taxonomy: Homo sapiens (9606); network_type=functional; add_nodes=0.",
        "- Primary network: combined score >=0.400. Sensitivity network: combined score >=0.700.",
        "- Primary input was all 84 frozen FDR-significant DEGs; no external node was added.", "",
        "## Mapping and topology", "",
        f"A. STRING mapping: mapped={mapped_n}/84; unmapped={unmapped_n}; ambiguous={ambiguous_n}. PDE2A mapping={pde_map['STRING_mapped_ID']} ({pde_map['mapping_status']}, {pde_map['mapping_ambiguity']}).",
        f"B. Score400: nodes=84; edges={graph400.number_of_edges()}; components={len(components)}; largest component={len(components[0])}; isolated nodes={len(isolates400)} ({';'.join(isolates400) if isolates400 else 'none'}).",
        f"C. Score700 sensitivity: nodes=84; edges={graph700.number_of_edges()}; components={len(components700)}; largest component={len(components700[0])}; isolated nodes={len(isolates700)}.",
        "D. Top degree nodes (topology only, not candidate definition): " + "; ".join(f"{r.gene}(degree={r.degree})" for r in top_hubs.itertuples()) + ".",
        f"E. Weighted Louvain was run once on the score400 largest component using seed={SEED}; communities={len(communities)}.",
        f"F. PDE2A: component={pde_node['component_id']}; community={pde_node['community_id']}; degree={pde_node['degree']}; strength={pde_node['weighted_degree_strength']:.3f}; degree rank={pde_degree_rank}/84; direct neighbors={';'.join(neighbors) if neighbors else 'none'}; retained with >=1 score700 edge={graph700.degree('PDE2A') > 0}; top-10-degree hub={pde_is_hub}.",
        f"G. PDE2A community functional context: {pde_comm_text}.",
        f"H. Edge-down signaling/regulatory functional-class candidates: {n_reg}. This differs from signaling/regulatory pathway-module membership ({module_count}/84).",
        "I. High-convergence comparison was descriptive and contains the prespecified naturally observed genes plus every Edge-down DEG with leading-edge count >= PDE2A; no PDE-family extension was performed.",
        "J. PDE2A network context was evaluated only after global network construction. It retains gene-level significance, moderate/large negative effect, 5/5 consistency, 11 direction-consistent leading edges, and enzymatic/regulatory annotation; network topology is reported independently.",
        "K. No machine learning, no composite score, no DEG threshold modification, no public-dataset intersection was performed.",
        "L. Centrality metrics describe topology only. Isolates retain degree/strength/betweenness=0 and component-wise closeness/eigenvector=NA.", "",
        "## Development figures", "",
        "- PNG 300 dpi plus source_data.tsv and plotting_script.R only; no PDF/TIFF.",
        "- Global and PDE2A-context networks use identical nodes, edges, positions, colors, and size rules; the context version adds only PDE2A label/outline.",
    ]
    (out / "QC_CHECKPOINT8_PPI_PRIORITIZATION_AUDIT.md").write_text("\n".join(audit) + "\n", encoding="utf-8")
    methods = {
        "python_version": sys.version.split()[0], "networkx_version": nx.__version__,
        "pandas_version": pd.__version__, "requests_version": requests.__version__,
        "STRING_version": string_version, "STRING_taxonomy": SPECIES,
        "STRING_network_type": "functional", "STRING_add_nodes": 0,
        "primary_required_score": 400, "sensitivity_required_score": 700,
        "community_algorithm": "networkx weighted Louvain on largest component",
        "random_seed": SEED,
    }
    (out / "CHECKPOINT8_METHODS_AND_VERSIONS.json").write_text(json.dumps(methods, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
