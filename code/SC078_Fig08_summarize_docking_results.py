from __future__ import annotations

import csv
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "PDE2A_docking_candidate_clinical_tables_20260807"
OUT.mkdir(parents=True, exist_ok=True)

PURPLE = "7C3AED"
PURPLE_LIGHT = "EDE9FE"
GRAY = "F3F4F6"
INK = "111827"
MUTED = "4B5563"
RED = "991B1B"


candidate_rows = [
    {
        "候选分子": "CGP-60474",
        "项目内筛选来源": "L1000CDS2逆转签名：strict 2个条件、wide 4个条件；retain_tier1",
        "已知主要机制": "CDK研究抑制剂（CDK1/2为主），不是已确证PDE2A抑制剂",
        "与PDE2A的直接证据": "not_tested：项目内无PDE2A酶学IC50/Ki、无靶标占有或竞争结合",
        "当前对接状态": "已对接1Z1L；但原候选审计表标记docking_eligible=No",
        "建议用途": "仅作表型逆转/探索性结合候选；必须与直接PDE2A抑制剂做实验对照",
    },
    {
        "候选分子": "Alvocidib (flavopiridol)",
        "项目内筛选来源": "L1000CDS2逆转签名：strict 1个条件、wide 4个条件；retain_tier1",
        "已知主要机制": "多靶点CDK抑制剂；项目定义为CDK translational lead",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested；本项目审计表标记docking_eligible=No",
        "建议用途": "CGP-60474的同类机制对照，检验逆转效应是否来自CDK轴",
    },
    {
        "候选分子": "BMS-536924",
        "项目内筛选来源": "L1000CDS2逆转签名：strict 1个条件、wide 1个条件；retain_tier2",
        "已知主要机制": "IGF-1R/胰岛素受体轴研究抑制剂",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested；本项目审计表标记docking_eligible=No",
        "建议用途": "探索性非CDK对照；需先补HNSCC敏感性和靶点依据",
    },
    {
        "候选分子": "Etoposide",
        "项目内筛选来源": "GDSC1+CTRP2表达-药敏投射；PDE2A-loss选择性Tier 1",
        "已知主要机制": "TOP2抑制/ DNA损伤；FDA批准肿瘤药",
        "与PDE2A的直接证据": "not_applicable：项目结论明确不支持直接PDE2A结合",
        "当前对接状态": "not_tested；不应为证明PDE2A直接结合而对接",
        "建议用途": "临床可及性比较物；检验PDE2A状态是否预测药敏，而非直接靶向PDE2A",
    },
    {
        "候选分子": "Bosutinib",
        "项目内筛选来源": "FDA再利用筛选Tier 2；跨数据库KD敏感，但存在暴露不足",
        "已知主要机制": "ABL/SRC/TEC激酶抑制剂；FDA批准CML药物",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested",
        "建议用途": "次级临床比较物；因预测浓度高于标示总Cmax，不宜主推",
    },
    {
        "候选分子": "Bleomycin",
        "项目内筛选来源": "FDA再利用筛选Tier 2；仅部分计算支持",
        "已知主要机制": "DNA断裂诱导；FDA批准抗肿瘤药",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_recommended：高分子量、复杂糖肽，不是本轮常规小分子对接优先对象",
        "建议用途": "药敏验证备选，不进入PDE2A直接结合主线",
    },
    {
        "候选分子": "Midostaurin",
        "项目内筛选来源": "FDA再利用筛选Tier 2；缺少独立通过的CTRP模型",
        "已知主要机制": "多激酶抑制剂；FDA批准AML/系统性肥大细胞增多症药物",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested",
        "建议用途": "低优先级药敏比较物，不作为PDE2A抑制剂",
    },
    {
        "候选分子": "Pyrimethamine",
        "项目内筛选来源": "FDA再利用筛选Tier 2；仅部分计算支持",
        "已知主要机制": "二氢叶酸还原酶抑制剂；FDA批准抗寄生虫药",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested",
        "建议用途": "低优先级再利用候选；需先做HN6剂量反应",
    },
    {
        "候选分子": "Pemetrexed",
        "项目内筛选来源": "FDA再利用筛选Tier 2；仅部分计算支持",
        "已知主要机制": "抗叶酸/胸苷酸合成酶轴；FDA批准肿瘤药",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested",
        "建议用途": "低优先级临床比较物，不进入直接结合主线",
    },
    {
        "候选分子": "Lenalidomide",
        "项目内筛选来源": "FDA再利用筛选Tier 2；仅部分计算支持",
        "已知主要机制": "CRBN依赖蛋白降解/免疫调节；FDA批准药",
        "与PDE2A的直接证据": "not_tested",
        "当前对接状态": "not_tested",
        "建议用途": "低优先级再利用比较物，不作为PDE2A抑制剂",
    },
]


docking_rows = [
    {
        "配体": "CGP-60474",
        "受体/方法": "PDE2A 1Z1L；AutoDock Vina 1.2.7；grid center=24.6/25.2/13.1 Å；exhaustiveness=9；20 modes",
        "对接分数": "best -8.166 kcal/mol；range -8.166至-7.179",
        "初始预测位点（重原子≤4.0 Å）": "TYR655, ASP808, GLN859, HIS656, ILE866, LEU770, HIS660, GLN812, ILE822, TYR827, PHE862, LEU809, ILE826",
        "100 ns MD持续接触（>50%帧）": "ASP808 97.9%; PHE830 92.3%; LEU770 89.7%; THR768 79.4%; ILE826 77.9%; TYR655 77.2%; ASN704 69.5%; HIS656 59.6%; PHE862 54.8%; ASP769 54.3%",
        "验证状态与结论": "未做共晶配体redocking；受体删除Zn/Mg/PO4；配体RMSD(20-100 ns)=8.44±0.45 Å。只能写predicted pose，不能写稳定原始结合构象或真实结合自由能。",
    },
    {
        "配体": "BAY 60-7550",
        "受体/方法": "项目内not_tested；文献共晶PDB 4HTX可用于阳性redocking",
        "对接分数": "not_tested",
        "初始预测位点（重原子≤4.0 Å）": "not_tested；文献结构锚点包括GLN859 glutamine switch与诱导型疏水口袋",
        "100 ns MD持续接触（>50%帧）": "not_tested",
        "验证状态与结论": "应作为标准化重对接的首选阳性对照；不能把文献IC50换算成Vina分数。",
    },
    {
        "配体": "PF-05180999",
        "受体/方法": "项目内not_tested",
        "对接分数": "not_tested",
        "初始预测位点（重原子≤4.0 Å）": "not_tested",
        "100 ns MD持续接触（>50%帧）": "not_tested",
        "验证状态与结论": "可作为不同化学骨架的高效PDE2A参考抑制剂；需在同一受体、同一网格、同一参数下重算后才可横向比较。",
    },
    {
        "配体": "TAK-915",
        "受体/方法": "项目内not_tested；有PDE2A共晶PDB 5VP0",
        "对接分数": "not_tested",
        "初始预测位点（重原子≤4.0 Å）": "not_tested",
        "100 ns MD持续接触（>50%帧）": "not_tested",
        "验证状态与结论": "可作为临床期、结构已知的第二阳性对照；优先用于共晶姿势恢复验证。",
    },
    {
        "配体": "Lu AF64280",
        "受体/方法": "项目内not_tested",
        "对接分数": "not_tested",
        "初始预测位点（重原子≤4.0 Å）": "not_tested",
        "100 ns MD持续接触（>50%帧）": "not_tested",
        "验证状态与结论": "选择性PDE2A工具化合物；目前仅可做药理效力参照，不能补写结合能。",
    },
    {
        "配体": "EHNA",
        "受体/方法": "项目内not_tested",
        "对接分数": "not_tested",
        "初始预测位点（重原子≤4.0 Å）": "not_tested",
        "100 ns MD持续接触（>50%帧）": "not_tested",
        "验证状态与结论": "历史PDE2抑制剂但强抑制ADA；只适合非选择性/历史对照，不宜作为首要阳性对照。",
    },
]


comparator_rows = [
    {
        "化合物": "CGP-60474",
        "直接PDE2A效力": "not_tested",
        "已知药理/选择性": "CDK1/cyclin B IC50 17 nM；CDK2/cyclin A 80 nM；CDK2/cyclin E 50 nM",
        "商品化状态": "是，科研级试剂（Cayman 29631）；非人用药",
        "临床阶段": "无已确认临床开发；不是FDA批准药",
        "本项目中的意义": "逆转签名/CDK探针。若不能在酶学和PDE2A依赖细胞实验中复现直接抑制，不得称PDE2A抑制剂。",
    },
    {
        "化合物": "BAY 60-7550",
        "直接PDE2A效力": "人重组PDE2 IC50 4.7±1.0 nM；牛心PDE2 2.0±0.7 nM",
        "已知药理/选择性": "约50倍优于PDE1；对多数其他PDE >100倍；共晶PDB 4HTX",
        "商品化状态": "是，科研级试剂（例如MCE HY-14992）",
        "临床阶段": "临床前工具化合物；未获批",
        "本项目中的意义": "首选阳性药理对照：验证PDE2A抑制是否能复制CGP-60474的cAMP/cGMP和细胞表型。",
    },
    {
        "化合物": "PF-05180999 (PF-999)",
        "直接PDE2A效力": "IC50 1.6 nM",
        "已知药理/选择性": ">2000倍选择性（相对其余PDE面板）；可升高脑cGMP",
        "商品化状态": "是，科研级试剂（Cayman 30374/MCE HY-111371）",
        "临床阶段": "Phase 1：NCT01429740、NCT01530529完成；NCT01981486在入组前撤回",
        "本项目中的意义": "第二化学骨架阳性对照；提供人体暴露历史，但不能外推到HNSCC疗效。",
    },
    {
        "化合物": "TAK-915",
        "直接PDE2A效力": "IC50 0.61 nM",
        "已知药理/选择性": ">4100倍相对PDE1A、>20000倍相对其他PDE；共晶PDB 5VP0",
        "商品化状态": "是，科研级试剂（MCE HY-103493）",
        "临床阶段": "Phase 1完成：NCT02461160；PET靶标占有研究NCT02584569完成",
        "本项目中的意义": "高选择性、临床期参照；适合确认PDE2A特异性和做结构重对接。",
    },
    {
        "化合物": "Lu AF64280",
        "直接PDE2A效力": "hPDE2A Ki 20 nM",
        "已知药理/选择性": "至少50倍优于PDE9A/PDE10A；对其他全长人PDE >250倍",
        "商品化状态": "not_verified（本轮未确认稳定供应商）",
        "临床阶段": "临床前工具化合物",
        "本项目中的意义": "文献药理参照；采购和批次可得性确认后再纳入实验。",
    },
    {
        "化合物": "EHNA",
        "直接PDE2A效力": "按体系约0.8-4 μM",
        "已知药理/选择性": "同时极强抑制腺苷脱氨酶（ADA，约2 nM），解释易混杂",
        "商品化状态": "是，科研级试剂（MCE HY-103160B/盐酸盐HY-103160A）",
        "临床阶段": "临床前/历史工具化合物；未获批",
        "本项目中的意义": "仅作历史或非选择性对照，不能替代BAY 60-7550/TAK-915。",
    },
    {
        "化合物": "Etoposide",
        "直接PDE2A效力": "not_applicable（不是PDE2A抑制剂）",
        "已知药理/选择性": "TOP2抑制；项目中为PDE2A-loss选择性计算候选",
        "商品化状态": "是；FDA批准肿瘤药",
        "临床阶段": "已上市，但不是PDE2A靶向治疗，也无本项目HNSCC疗效证据",
        "本项目中的意义": "临床可及性比较物：用于PDE2A状态-药敏关系，不用于证明直接结合。",
    },
]


references = [
    "Duan Q, et al. L1000CDS2: LINCS L1000 characteristic direction signatures search engine. npj Syst Biol Appl. 2016;2:16015. https://doi.org/10.1038/npjsba.2016.15",
    "Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell. 2017;171:1437-1452.e17. https://doi.org/10.1016/j.cell.2017.10.049",
    "Iorio F, et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell. 2016;166:740-754. https://doi.org/10.1016/j.cell.2016.06.017",
    "Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. Brief Bioinform. 2021;22:bbab260. https://doi.org/10.1093/bib/bbab260",
    "Eberhardt J, et al. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. J Chem Inf Model. 2021;61:3891-3898. https://doi.org/10.1021/acs.jcim.1c00203",
    "Boess FG, et al. Inhibition of phosphodiesterase 2 increases neuronal cGMP, synaptic plasticity and memory performance. Neuropharmacology. 2004;47:1081-1092. https://doi.org/10.1016/j.neuropharm.2004.07.040",
    "Zhu J, Yang Q, Dai D, Huang Q. X-ray crystal structure of phosphodiesterase 2 in complex with a highly selective, nanomolar inhibitor reveals a binding-induced pocket important for selectivity. J Am Chem Soc. 2013;135:11708-11711. https://doi.org/10.1021/ja404449g",
    "Helal CJ, et al. Identification of a Potent, Highly Selective, and Brain Penetrant Phosphodiesterase 2A Inhibitor Clinical Candidate. J Med Chem. 2018;61:1001-1018. https://doi.org/10.1021/acs.jmedchem.7b01466",
    "Mikami S, et al. Discovery of Clinical Candidate TAK-915: A Highly Potent, Selective, and Brain-Penetrating PDE2A Inhibitor. J Med Chem. 2017;60:7677-7702. https://doi.org/10.1021/acs.jmedchem.7b00807",
    "Redrobe JP, et al. In vitro and in vivo characterisation of Lu AF64280. Psychopharmacology. 2014;231:3151-3167. https://doi.org/10.1007/s00213-014-3492-7",
    "ClinicalTrials.gov: PF-05180999 NCT01429740, NCT01530529, NCT01981486; TAK-915 NCT02461160, NCT02584569. Accessed 2026-08-07.",
    "Commercial availability checked 2026-08-07 using Cayman Chemical and MedChemExpress product pages. 'Commercial' means research-use reagent only, not regulatory approval.",
]


def write_tsv(name: str, rows: list[dict[str, str]]) -> None:
    path = OUT / name
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def md_table(rows: list[dict[str, str]]) -> str:
    headers = list(rows[0])
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        vals = [str(row[h]).replace("|", "\\|").replace("\n", "<br>") for h in headers]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_markdown() -> None:
    text = f"""# PDE2A分子对接候选、结合位点与临床转化证据表

检索与核对日期：2026-08-07。本文档为项目现有计算结果与定向文献核对，不是系统综述。

## 核心结论

1. CGP-60474是PDE2A-KD/OE转录逆转筛选得到的CDK研究候选，不是已验证的PDE2A抑制剂。项目内候选审计表甚至将其标记为`docking_eligible=No`，因此对接只能称为探索性假设生成。
2. 当前唯一真实完成的对接是CGP-60474-PDE2A 1Z1L，Vina最佳分数为-8.166 kcal/mol。该分数不是实验结合自由能；受体删除了催化Zn/Mg/PO4且无共晶redocking，不能据此声称直接抑制。
3. 与商品化抑制剂的正确比较应分两层：BAY 60-7550/TAK-915/PF-05180999用于PDE2A特异性和靶标占有验证；etoposide用于临床可及的PDE2A状态-药敏比较。两类比较不能混写。
4. 本地Drugs@FDA快照（更新至2026-07-31）未检出BAY 60-7550、PF-05180999、TAK-915、Lu AF64280、EHNA或CGP-60474。现有选择性PDE2A抑制剂仍应写成科研级或临床期化合物，不能写成获批PDE2A药物，更不能写成HNSCC临床治疗。

## Table 1. 项目候选小分子及筛选依据

{md_table(candidate_rows)}

表注：L1000CDS2依据PDE2A扰动表达签名寻找方向相反的小分子扰动；这种方法支持“表型逆转候选”，不支持“直接PDE2A结合”。GDSC/CTRP表达-药敏模型支持“PDE2A状态相关药敏候选”，同样不支持直接靶向。[1-4]

## Table 2. 对接分数、预测位点与MD接触证据

{md_table(docking_rows)}

表注：初始位点由当前top-1 PDBQT构象中配体-受体重原子距离≤4.0 Å重新计算；MD接触为全部3334帧中重原子距离阈值定义的占有率。Vina分数不能解释为实验ΔG，且不同受体/参数下的分数不得横向排名。[5]

## Table 3. 商品化/临床期PDE2A抑制剂与临床比较物

{md_table(comparator_rows)}

表注：IC50和Ki受酶构建体、底物和检测体系影响，只有同一实验体系的数值才是严格头对头比较。“商品化”仅表示科研级试剂可购买，不表示获批、可用于患者或具有HNSCC临床证据。[6-11]

## 分子筛选方法的规范写法

PDE2A敲低和过表达转录组先经扰动方向质控，再从可逆差异基因构建strict和wide两套查询签名。两套签名分别提交L1000CDS2，优先保留在多个细胞/条件中重复出现且在两套签名中方向一致的逆转化合物。CGP-60474分别在strict和wide查询中出现2和4个LINCS条件，因此作为CDK机制探针保留；alvocidib和BMS-536924作为备选。该流程依据Connectivity Map/L1000的“疾病或扰动签名被药物反向调节”思想，提供表型相关性而非靶标结合证据。[1,2]

独立的临床再利用层使用GDSC/CTRP细胞系表达-药敏资源构建质量门控的岭回归投射，并要求PDE2A敲低敏感、过表达方向相反及跨数据库一致。该层筛得etoposide为Tier 1，但其含义是“PDE2A-loss选择性脆弱性候选”，不是“PDE2A抑制剂”。[3,4]

## 与商品化抑制剂比较的最小实验闭环

1. **直接酶学**：在同一批人重组PDE2A和同一底物体系测CGP-60474、BAY 60-7550、TAK-915或PF-05180999的完整剂量反应，报告IC50、95% CI、Hill slope和独立重复数。
2. **细胞靶标效应**：在HN6 NC、PDE2A-KD、Vector和PDE2A-OE中，于短时程测cAMP/cGMP；BAY 60-7550和TAK-915/PF-05180999应作为两个不同骨架的阳性对照。
3. **表型效应**：72 h细胞活力/克隆形成比较CGP-60474、直接PDE2A抑制剂和etoposide，报告IC50/Emax及`drug × PDE2A-state`交互；不能只比单一浓度。
4. **机制拆分**：同步测p-RB、细胞周期和CDK轴。若CGP-60474抑制增殖但不抑制PDE2A酶活、不能被PDE2A过表达救援，结论应回到“CDK介导的表型逆转”，删除直接PDE2A结合叙述。
5. **进入文章的停止规则**：只有在共晶配体redocking成功、保留/合理参数化金属、参考抑制剂同协议对接、至少两条独立MD重复和酶学/细胞靶标实验均支持时，才升级为“candidate direct PDE2A inhibitor”。否则维持“exploratory docking candidate”。

## 可直接用于正文的结论

**英文版：** Transcriptomic reversal screening nominated CGP-60474 as a phenotype-oriented CDK probe rather than an established PDE2A inhibitor. Exploratory docking predicted a top AutoDock Vina score of -8.166 kcal/mol in the metal-depleted 1Z1L receptor model, while the ligand deviated substantially from the initial pose during the 100-ns simulation. These findings support hypothesis generation but do not establish direct PDE2A inhibition. Biochemical and cellular benchmarking against commercially available selective PDE2A inhibitors, preferably BAY 60-7550 plus TAK-915 or PF-05180999, is therefore required. Etoposide should be interpreted separately as an FDA-approved, PDE2A-loss-associated drug-sensitivity candidate rather than a PDE2A-binding comparator.

**中文版：** 转录组逆转筛选将CGP-60474定位为表型相关的CDK研究探针，而非已确证的PDE2A抑制剂。探索性对接在去除金属离子的1Z1L受体模型中得到-8.166 kcal/mol的最佳Vina分数，但配体在100 ns模拟中明显偏离初始姿势。因此，现有结果仅支持提出直接结合假设，尚不能证明PDE2A抑制。后续必须以商品化选择性PDE2A抑制剂BAY 60-7550，并加入TAK-915或PF-05180999作为第二化学骨架对照，完成酶学和细胞靶标效应验证。Etoposide应单独解释为FDA批准的PDE2A-loss相关药敏候选，而不是PDE2A结合对照。

## 参考文献

"""
    text += "\n".join(f"{i}. {ref}" for i, ref in enumerate(references, 1)) + "\n"
    (OUT / "PDE2A_docking_selection_and_clinical_interpretation.md").write_text(text, encoding="utf-8")


def write_xlsx() -> None:
    wb = Workbook()
    wb.remove(wb.active)
    datasets = [
        ("Table1_candidates", candidate_rows),
        ("Table2_docking_sites", docking_rows),
        ("Table3_comparators", comparator_rows),
    ]
    thin = Side(style="thin", color="D1D5DB")
    for title, rows in datasets:
        ws = wb.create_sheet(title)
        headers = list(rows[0])
        ws.append(headers)
        for row in rows:
            ws.append([row[h] for h in headers])
        for cell in ws[1]:
            cell.font = Font(name="Microsoft YaHei", bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor=PURPLE)
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for row in ws.iter_rows(min_row=2):
            for cell in row:
                cell.font = Font(name="Microsoft YaHei", size=10, color=INK)
                cell.alignment = Alignment(vertical="top", wrap_text=True)
                cell.border = Border(left=thin, right=thin, top=thin, bottom=thin)
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        widths = [20, 34, 30, 34, 34, 42]
        for idx, width in enumerate(widths[: len(headers)], 1):
            ws.column_dimensions[get_column_letter(idx)].width = width
        ws.row_dimensions[1].height = 36
        for idx in range(2, ws.max_row + 1):
            ws.row_dimensions[idx].height = 72
    ref_ws = wb.create_sheet("References")
    ref_ws.append(["No.", "Reference"])
    for i, ref in enumerate(references, 1):
        ref_ws.append([i, ref])
    ref_ws.column_dimensions["A"].width = 8
    ref_ws.column_dimensions["B"].width = 120
    ref_ws.freeze_panes = "A2"
    for cell in ref_ws[1]:
        cell.font = Font(name="Arial", bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor=PURPLE)
    for row in ref_ws.iter_rows():
        for cell in row:
            cell.alignment = Alignment(vertical="top", wrap_text=True)
    wb.save(OUT / "PDE2A_docking_candidate_clinical_tables.xlsx")


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=90, start=110, bottom=90, end=110) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for m, v in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{m}"))
        if node is None:
            node = OxmlElement(f"w:{m}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(v))
        node.set(qn("w:type"), "dxa")


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def set_table_widths(table, widths: list[float]) -> None:
    table.autofit = False
    tbl_pr = table._tbl.tblPr
    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    total_dxa = int(sum(widths) * 1440)
    tbl_w.set(qn("w:w"), str(total_dxa))
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), "110")
    tbl_ind.set(qn("w:type"), "dxa")
    grid = table._tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(int(width * 1440)))
        grid.append(col)
    for row in table.rows:
        for cell, width in zip(row.cells, widths):
            cell.width = Inches(width)
            tc_w = cell._tc.get_or_add_tcPr().find(qn("w:tcW"))
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                cell._tc.get_or_add_tcPr().append(tc_w)
            tc_w.set(qn("w:w"), str(int(width * 1440)))
            tc_w.set(qn("w:type"), "dxa")


def set_run_font(run, size=9, bold=False, color=INK) -> None:
    run.font.name = "Arial"
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), "Arial")
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), "Arial")
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    run.font.size = Pt(size)
    run.bold = bold
    run.font.color.rgb = RGBColor.from_string(color)


def add_table(doc: Document, rows: list[dict[str, str]], widths: list[float]) -> None:
    headers = list(rows[0])
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    set_table_widths(table, widths)
    header = table.rows[0]
    set_repeat_table_header(header)
    for cell, label in zip(header.cells, headers):
        set_cell_shading(cell, PURPLE)
        set_cell_margins(cell)
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(0)
        set_run_font(p.add_run(label), size=8.5, bold=True, color="FFFFFF")
    for ridx, row in enumerate(rows, 1):
        cells = table.add_row().cells
        for cell, header_name in zip(cells, headers):
            set_cell_margins(cell)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            if ridx % 2 == 0:
                set_cell_shading(cell, "FAFAFA")
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.08
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            set_run_font(p.add_run(row[header_name]), size=8.2, color=INK)


def add_heading(doc: Document, text: str, level: int = 1) -> None:
    p = doc.add_paragraph(style=f"Heading {level}")
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.left_indent = Inches(0)
    p.paragraph_format.right_indent = Inches(0)
    p.paragraph_format.first_line_indent = Inches(0)
    p.paragraph_format.space_before = Pt(12 if level == 1 else 8)
    p.paragraph_format.space_after = Pt(5)
    run = p.add_run(text)
    set_run_font(run, size=15 if level == 1 else 12, bold=True, color=PURPLE)


def add_body(doc: Document, text: str, bold_prefix: str | None = None) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(5)
    p.paragraph_format.line_spacing = 1.15
    if bold_prefix and text.startswith(bold_prefix):
        set_run_font(p.add_run(bold_prefix), size=10.5, bold=True)
        text = text[len(bold_prefix) :]
    set_run_font(p.add_run(text), size=10.5)


def add_numbered(doc: Document, items: list[str]) -> None:
    for idx, item in enumerate(items, 1):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.25)
        p.paragraph_format.first_line_indent = Inches(-0.25)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing = 1.12
        set_run_font(p.add_run(f"{idx}. "), size=10, bold=True, color=PURPLE)
        set_run_font(p.add_run(item), size=10)


def write_docx() -> None:
    doc = Document()
    section = doc.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = Inches(11), Inches(8.5)
    section.top_margin = Inches(0.62)
    section.bottom_margin = Inches(0.62)
    section.left_margin = Inches(0.55)
    section.right_margin = Inches(0.55)
    section.header_distance = Inches(0.3)
    section.footer_distance = Inches(0.3)

    normal = doc.styles["Normal"]
    normal.font.name = "Arial"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor.from_string(INK)

    title = doc.add_paragraph()
    title.paragraph_format.space_after = Pt(4)
    run = title.add_run("PDE2A分子对接候选、结合位点与临床转化证据表")
    set_run_font(run, size=21, bold=True, color=PURPLE)
    subtitle = doc.add_paragraph()
    subtitle.paragraph_format.space_after = Pt(12)
    set_run_font(subtitle.add_run("项目证据审计与文章可用表格 | 核对日期：2026-08-07"), size=10.5, color=MUTED)

    add_body(
        doc,
        "结论：CGP-60474是L1000逆转签名得到的CDK研究候选，不是已验证的PDE2A抑制剂。现有-8.166 kcal/mol Vina分数来自去除Zn/Mg/PO4且未做redocking验证的1Z1L模型，只能支持探索性预测。商品化PDE2A抑制剂应作为实验阳性对照，而不是用其文献IC50与Vina分数直接比较。",
        "结论：",
    )

    add_heading(doc, "Table 1. 项目候选小分子及筛选依据")
    add_table(doc, candidate_rows, [1.00, 1.65, 1.45, 1.60, 1.60, 2.50])
    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    note.paragraph_format.space_after = Pt(5)
    set_run_font(note.add_run("表注："), size=8.5, bold=True, color=PURPLE)
    set_run_font(note.add_run("L1000逆转与GDSC/CTRP药敏投射分别支持表型逆转和状态相关药敏，不等于直接靶向PDE2A。[1-4]"), size=8.5, color=MUTED)

    doc.add_page_break()
    add_heading(doc, "Table 2. 对接分数、预测位点与MD接触证据")
    add_table(doc, docking_rows, [0.95, 1.40, 1.05, 1.90, 2.15, 2.35])
    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    set_run_font(note.add_run("表注："), size=8.5, bold=True, color=PURPLE)
    set_run_font(note.add_run("初始位点按top-1构象重原子距离≤4.0 Å计算。Vina score不是实验结合自由能；未在同一协议计算的分数不得横向排名。[5]"), size=8.5, color=MUTED)

    doc.add_page_break()
    add_heading(doc, "Table 3. 商品化/临床期PDE2A抑制剂与临床比较物")
    add_table(doc, comparator_rows, [1.00, 1.30, 1.55, 1.30, 1.65, 3.00])
    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(4)
    set_run_font(note.add_run("表注："), size=8.5, bold=True, color=PURPLE)
    set_run_font(note.add_run("IC50/Ki受实验体系影响；商品化仅指科研级试剂可购买，不表示获批或具有HNSCC临床证据。[6-11]"), size=8.5, color=MUTED)

    doc.add_page_break()
    add_heading(doc, "筛选方法与临床解释")
    add_body(doc, "PDE2A敲低和过表达转录组经方向质控后构建strict/wide可逆签名，并用L1000CDS2寻找方向相反且跨条件复现的化合物。CGP-60474在两套查询中分别出现2和4个条件，适合作为CDK机制探针，但该算法不证明直接PDE2A结合。[1,2]")
    add_body(doc, "独立的FDA再利用层依据GDSC/CTRP表达-药敏数据与质量门控模型识别PDE2A状态相关药敏。Etoposide为Tier 1计算候选，其含义是PDE2A-loss选择性脆弱性，而不是PDE2A抑制剂。[3,4]")
    add_heading(doc, "最小实验闭环", level=2)
    add_numbered(
        doc,
        [
            "同一人重组PDE2A体系测CGP-60474、BAY 60-7550及TAK-915或PF-05180999的完整剂量反应，报告IC50、95% CI、Hill slope和独立重复数。",
            "在HN6 NC、PDE2A-KD、Vector和PDE2A-OE中短时测cAMP/cGMP，确认药物是否产生PDE2A依赖靶标效应。",
            "72 h活力/克隆形成比较CGP-60474、直接PDE2A抑制剂和etoposide，报告IC50/Emax及drug × PDE2A-state交互。",
            "同步检测p-RB和细胞周期。若CGP-60474仅表现CDK效应而无PDE2A酶学/救援证据，应删除直接PDE2A抑制叙述。",
            "重做结构计算时使用4HTX/5VP0共晶配体redocking、保留并合理参数化金属、同协议对接全部参考物，并增加独立MD重复。",
        ],
    )
    add_heading(doc, "可直接用于正文的结论", level=2)
    add_body(doc, "Transcriptomic reversal screening nominated CGP-60474 as a phenotype-oriented CDK probe rather than an established PDE2A inhibitor. Exploratory docking predicted a top AutoDock Vina score of -8.166 kcal/mol in the metal-depleted 1Z1L receptor model, while the ligand deviated substantially from the initial pose during the 100-ns simulation. These findings support hypothesis generation but do not establish direct PDE2A inhibition. Biochemical and cellular benchmarking against commercially available selective PDE2A inhibitors, preferably BAY 60-7550 plus TAK-915 or PF-05180999, is therefore required. Etoposide should be interpreted separately as an FDA-approved, PDE2A-loss-associated drug-sensitivity candidate rather than a PDE2A-binding comparator.")

    doc.add_page_break()
    add_heading(doc, "参考文献")
    for idx, ref in enumerate(references, 1):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.25)
        p.paragraph_format.first_line_indent = Inches(-0.25)
        p.paragraph_format.space_after = Pt(3)
        set_run_font(p.add_run(f"{idx}. {ref}"), size=8.8)

    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    set_run_font(footer.add_run("PDE2A docking evidence audit | 2026-08-07"), size=8, color=MUTED)
    doc.save(OUT / "PDE2A_docking_candidate_clinical_tables.docx")


def main() -> None:
    write_tsv("Table1_candidate_small_molecules.tsv", candidate_rows)
    write_tsv("Table2_docking_energy_binding_sites.tsv", docking_rows)
    write_tsv("Table3_reference_inhibitor_clinical_comparison.tsv", comparator_rows)
    write_markdown()
    write_xlsx()
    write_docx()


if __name__ == "__main__":
    main()
