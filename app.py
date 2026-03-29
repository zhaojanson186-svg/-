import streamlit as st
import re
import pandas as pd
import io
import requests

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="抗体序列核心分析中台", page_icon="🧬", layout="wide")

st.title("🧬 高通量抗体序列处理与成药性分析中台 (纯净版)")
st.info("💡 专注核心生信管线：极速进行序列提纯、CDR/FR 区域智能解析、高危 PTM (翻译后修饰) 扫描、以及 CMC 理化性质聚类。")

# ==========================================
# 2. 深度 CMC 评估与序列分析引擎
# ==========================================
def calculate_pi(seq):
    seq = seq.upper()
    pKa_acidic = {'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1}
    pKa_basic = {'K': 10.8, 'R': 12.5, 'H': 6.5}
    def net_charge(pH):
        charge = 1.0 / (1.0 + 10**(pH - 8.0)) - 1.0 / (1.0 + 10**(3.1 - pH))
        for aa, pka in pKa_acidic.items(): charge -= seq.count(aa) / (1.0 + 10**(pka - pH))
        for aa, pka in pKa_basic.items(): charge += seq.count(aa) / (1.0 + 10**(pH - pka))
        return charge
    pH_min, pH_max = 0.0, 14.0
    for _ in range(100):
        pH_mid = (pH_min + pH_max) / 2
        if net_charge(pH_mid) > 0: pH_min = pH_mid
        else: pH_max = pH_mid
    return round(pH_mid, 2)

def calculate_gravy(seq):
    if not seq or seq == "未识别": return 0.0
    kd_scale = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    return round(sum(kd_scale.get(aa, 0.0) for aa in seq.upper()) / len(seq), 3)

def detect_unpaired_cysteine(seq):
    cys_positions = [i+1 for i, aa in enumerate(seq.upper()) if aa == 'C']
    count = len(cys_positions)
    if count % 2 != 0: return f"🚨 高危: 奇数({count})个 Cys @{cys_positions}"
    elif count > 2 and count % 2 == 0: return f"⚠️ 警告: 额外配对({count})个 Cys @{cys_positions}"
    return "✅ 正常 (2个 Cys)"

def guess_germline(seq):
    seq = seq.upper()[:50]
    if re.search(r'GGGSVQ', seq) or re.search(r'W[FY]RQAPGKERE', seq): return "Camelid VHH (纳米抗体)"
    if re.search(r'SGGGLVQ', seq): return "Human IGHV3 (高人源化)"
    if re.search(r'SGAEVKKPG', seq): return "Human IGHV1/5 (高人源化)"
    if re.search(r'SGSELKKPG', seq): return "Human IGHV7 (高人源化)"
    if re.search(r'SGPGLVKPSG', seq): return "Human IGHV4 (高人源化)"
    if re.search(r'SGPEVKKPG', seq): return "Human IGHV2 (高人源化)"
    if re.search(r'QQSG[AP]E[LV]V', seq) or re.search(r'QQSDA', seq) or re.search(r'GSLKLS', seq): return "Murine IGHV (鼠源)"
    if re.search(r'VQL[VQE]QSG', seq) or re.search(r'VQL[LVE]ESG', seq): return "IGHV (亚族未定)"
    if re.search(r'SP[SS][SF]LSASVG', seq): return "Human IGKV1/3 (高人源化)"
    if re.search(r'SPLSLPVTPG', seq): return "Human IGKV2 (高人源化)"
    if re.search(r'SP[DS]SLA[VS]SLG', seq): return "Human IGKV4 (高人源化)"
    if re.search(r'SP[SA]YLAASP', seq) or re.search(r'FMSTSVG', seq): return "Murine IGKV (鼠源)"
    if re.search(r'[DE][IV][VQAM][ML]TQS', seq): return "IGKV (亚族未定)"
    if re.search(r'QPPS[AS]SG', seq) or re.search(r'Q[PP]SVS[VAS]P', seq): return "Human IGLV (Lambda)"
    if re.search(r'LTQP', seq): return "IGLV (Lambda未定)"
    if seq.startswith('Q') or seq.startswith('E'): return "疑似重链 (高度变异)"
    if seq.startswith('D') or seq.startswith('A'): return "疑似轻链 (高度变异)"
    return "未知架构"

def get_region_finder(seq, cdrs, domain_type):
    if "Fc" in domain_type: return lambda i: "Fc区"
    c1, c2, c3 = cdrs.get("CDR1",""), cdrs.get("CDR2",""), cdrs.get("CDR3","")
    idx1 = seq.find(c1) if c1 != "未识别" else -1
    idx2 = seq.find(c2, max(0, idx1)) if c2 != "未识别" else -1
    idx3 = seq.find(c3, max(0, idx2)) if c3 != "未识别" else -1
    
    def region_of(i):
        if idx1 != -1 and i < idx1: return "FR1"
        if idx1 != -1 and i < idx1 + len(c1): return "CDR1"
        if idx1 != -1 and idx2 != -1 and i < idx2: return "FR2"
        if idx2 != -1 and i < idx2 + len(c2): return "CDR2"
        if idx2 != -1 and idx3 != -1 and i < idx3: return "FR3"
        if idx3 != -1 and i < idx3 + len(c3): return "CDR3"
        if idx3 != -1 and i >= idx3 + len(c3): return "FR4"
        return "可变区"
    return region_of

def detect_ptms_detailed(seq, cdrs, domain_type):
    region_finder = get_region_finder(seq, cdrs, domain_type)
    ptm_rules = {"N-糖基化": r"N[^P][ST]", "脱氨基": r"N[GSN]", "异构化": r"D[GS]", "酸断裂": r"DP", "氧化": r"M"}
    found_ptms = []
    for ptm_name, pattern in ptm_rules.items():
        for match in re.finditer(pattern, seq):
            found_ptms.append(f"[{region_finder(match.start())}] {ptm_name}({match.group()}) @{match.start()+1}")
    found_ptms.sort(key=lambda x: int(re.search(r'@(\d+)', x).group(1)) if re.search(r'@(\d+)', x) else 0)
    return " | ".join(found_ptms) if found_ptms else "✅ 无常见高危 PTM"

def extract_cdrs_via_api(seq):
    api_url = "https://api.antibody-informatics.org/v1/anarci/annotate"
    payload = {"sequence": seq, "scheme": "imgt"}
    try:
        response = requests.post(api_url, json=payload, headers={"Content-Type": "application/json"}, timeout=3)
        if response.status_code == 200:
            data = response.json()
            return { "CDR1": data.get("CDR1", "未识别"), "CDR2": data.get("CDR2", "未识别"), "CDR3": data.get("CDR3", "未识别") }
        else: raise Exception("API 异常")
    except Exception:
        if "M" in seq or "L" in seq:
             return extract_vh_cdrs_regex(seq) if seq.startswith("E") or seq.startswith("Q") else extract_vl_cdrs_regex(seq)
        return extract_vh_cdrs_regex(seq)

def extract_vh_cdrs_regex(vh_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    cdr3_match = re.search(r"Y[YFCA]C[A-Z]{1,3}(.*?)W[GS][A-Z][GTSVI]", vh_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    cdr1_match = re.search(r"C[A-Z]{2,6}(.{5,16})W[A-Z][RQK]", vh_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    cdr2_match = re.search(r"(?:EW[IVMASTL][A-Z]|KWM[A-Z]|REG[VLIA][A-Z]|RWV[A-Z])(.{8,30}?)[RKQ][VFSTIAM][TVIAMFSC][A-Z]?", vh_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    return cdrs

def extract_vl_cdrs_regex(vl_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    cdr3_match = re.search(r"Y[YFCA]C(.*?)(?:F[GSA][A-Z][GTV]|FGC)", vl_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    cdr1_match = re.search(r"C(.{8,18})W[YFL]", vl_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq)
    if not cdr2_match: cdr2_match = re.search(r"W[YFL].{10,22}?([A-Z]{7})G[A-Z]{1,2}[RFS]", vl_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    return cdrs

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名序列_1"] = re.sub(r'\s+', '', text).upper()
        return sequences
    for part in text.split(">"):
        if not part.strip(): continue
        lines = part.strip().split("\n")
        name, seq = lines[0].strip(), "".join(lines[1:]).replace(" ", "").upper()
        if name and seq: sequences[name] = seq
    return sequences

# ==========================================
# 3. 核心功能：数据预处理 (Excel 转 FASTA)
# ==========================================
st.markdown("---")
st.markdown("### 🗂️ 模块一：序列预处理清洗器")
st.caption("如果你有从本地测序公司或高通量筛选结果拿到的 Excel 表格，可以在此一键提取并标准化为 FASTA 格式。")

uploaded_excel = st.file_uploader("📤 请上传包含重/轻链序列的 Excel 或 CSV 文件", type=['xlsx', 'xls', 'csv'], key="excel_uploader")
if uploaded_excel is not None:
    if st.button("🔄 一键清洗并转换为 FASTA 格式", type="primary"):
        try:
            if uploaded_excel.name.endswith('.csv'): df_raw = pd.read_csv(uploaded_excel, header=None)
            else: df_raw = pd.read_excel(uploaded_excel, header=None)
            
            fasta_lines = []
            invalid_names = ['nan', 'None', '', 'Protein name', 'Humanized variants', 'Heavy chain', 'Light chain', 'Sequence']
            for index, row in df_raw.iterrows():
                name = str(row.iloc[0]).strip() if len(row) > 0 else ""
                hc = str(row.iloc[1]).strip() if len(row) > 1 else ""
                lc = str(row.iloc[2]).strip() if len(row) > 2 else ""
                
                if name in invalid_names or "sequence" in name.lower() or name.startswith('Unnamed:'): continue
                if hc and hc not in ['nan', 'None'] and len(hc) > 10: fasta_lines.append(f">{name}_Heavy_Chain\n{hc.upper()}")
                if lc and lc not in ['nan', 'None'] and len(lc) > 10: fasta_lines.append(f">{name}_Light_Chain\n{lc.upper()}")
            
            fasta_content = '\n'.join(fasta_lines)
            if fasta_content:
                st.success(f"✅ 清洗完成！成功提取出 {fasta_content.count('>')} 条链序列。你可以直接全选下方文本复制到【模块二】中进行深度解析。")
                st.text_area("📋 提取出的标准化 FASTA 序列 (支持全选复制):", value=fasta_content, height=200)
            else:
                st.warning("⚠️ 未能从文件中提取出有效序列，请确保第1列是名称，第2/3列是氨基酸序列。")
        except Exception as e:
            st.error(f"❌ 读取文件错误: {e}")

# ==========================================
# 4. 核心功能：高通量解析与成药性打分
# ==========================================
st.markdown("---")
st.markdown("### 🔬 模块二：序列多维解析与 CMC 排雷引擎")

col1, col2 = st.columns([3, 1])
with col1:
    raw_input = st.text_area("📥 请在此粘贴需要评估的候选序列 (支持混合文本或 FASTA 格式):", height=250, key="main_input")
with col2:
    st.markdown("##### ⚙️ 引擎设置")
    engine_choice = st.radio("选择底层提取引擎：", 
                             ("🚀 本地正则引擎 (Regex)\n极速/零延迟/断网可用", 
                              "☁️ 外部云端对齐 (API)\n更高精度/依赖网络"))
    
    analyze_btn = st.button("🚀 启动深度解析", type="primary", use_container_width=True)

if analyze_btn:
    use_api = "API" in engine_choice
    
    if raw_input:
        seq_dict = parse_fasta(raw_input)
        st.success(f"✅ 成功读取 {len(seq_dict)} 条输入序列，正在调用底层引擎扫库...")
        
        all_results = []
        vh_pattern = re.compile(r"([EQ].{100,135}VTVSS)")
        vl_pattern = re.compile(r"([DEQA].{95,125}(?:VEIK|LEIK|TVLG|VTVL|FGC))")
        fc_pattern = re.compile(r"(CPPCP.*?LSPGK)")
        
        progress_bar = st.progress(0)
        total_seqs = len(seq_dict)
        
        for idx, (seq_name, clean_seq) in enumerate(seq_dict.items()):
            for i, vh in enumerate(vh_pattern.findall(clean_seq)):
                cdrs = extract_cdrs_via_api(vh) if use_api else extract_vh_cdrs_regex(vh)
                comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
                all_results.append({
                    "序列名称 (ID)": seq_name, "区域": f"VH_{i+1}", "类型": "重链/纳米抗体",
                    "同源 Germline": guess_germline(vh), "孤立Cys 雷达": detect_unpaired_cysteine(vh),
                    "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vh),
                    "PTM 风险预警": detect_ptms_detailed(vh, cdrs, "VH"),
                    "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vh
                })
            for i, vl in enumerate(vl_pattern.findall(clean_seq)):
                cdrs = extract_cdrs_via_api(vl) if use_api else extract_vl_cdrs_regex(vl)
                comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
                all_results.append({
                    "序列名称 (ID)": seq_name, "区域": f"VL_{i+1}", "类型": "轻链",
                    "同源 Germline": guess_germline(vl), "孤立Cys 雷达": detect_unpaired_cysteine(vl),
                    "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vl),
                    "PTM 风险预警": detect_ptms_detailed(vl, cdrs, "VL"),
                    "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vl
                })
            for i, fc in enumerate(fc_pattern.findall(clean_seq)):
                all_results.append({
                    "序列名称 (ID)": seq_name, "区域": f"Fc_{i+1}", "类型": "Fc区",
                    "同源 Germline": "IgG Fc", "孤立Cys 雷达": "-", "CDR_GRAVY": "-",
                    "pI (等电点)": calculate_pi(fc), "PTM 风险预警": detect_ptms_detailed(fc, {}, "Fc"),
                    "CDR1": "-", "CDR2": "-", "CDR3": "-", "完整序列": fc
                })
            progress_bar.progress((idx + 1) / total_seqs)
        
        if all_results:
            df = pd.DataFrame(all_results)
            
            # --- 数据高亮与美化 ---
            def highlight_alerts(val):
                val_str = str(val)
                if '🚨' in val_str or '高危' in val_str or '脱氨基' in val_str or '异构化' in val_str:
                    return 'background-color: #ffcccc; color: #990000; font-weight: bold'
                elif '⚠️' in val_str or '氧化' in val_str:
                    return 'background-color: #fff2cc; color: #b26b00'
                return ''
            
            st.markdown("#### 📊 CMC 序列解析总表")
            st.dataframe(df.style.map(highlight_alerts, subset=['孤立Cys 雷达', 'PTM 风险预警']), use_container_width=True)
            
            # --- 序列聚类统计 ---
            df_v = df[df['类型'].isin(['重链/纳米抗体', '轻链'])]
            cluster_v = df_v.groupby('完整序列').agg(
                链类型=('类型', 'first'), 相同序列数=('序列名称 (ID)', 'count'),
                来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x.unique())),
                PTM风险=('PTM 风险预警', 'first'), CDR3=('CDR3', 'first')
            ).reset_index().sort_values(by=['链类型', '相同序列数'], ascending=[True, False])
            
            df_valid_cdr3 = df_v[~df_v['CDR3'].str.contains('未识别|失败', na=False)]
            cluster_cdr3 = df_valid_cdr3.groupby('CDR3').agg(
                链类型=('类型', 'first'), 共享该CDR3数量=('序列名称 (ID)', 'count'),
                来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x.unique())), 代表完整序列=('完整序列', 'first')
            ).reset_index().sort_values(by=['链类型', '共享该CDR3数量'], ascending=[True, False])
            
            # --- 导出模块 ---
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.to_excel(writer, index=False, sheet_name='解析总表_Total_Report')
                if not cluster_v.empty: cluster_v.to_excel(writer, index=False, sheet_name='全序列聚类去重')
                if not cluster_cdr3.empty: cluster_cdr3.to_excel(writer, index=False, sheet_name='基于CDR3同源聚类')
            
            st.markdown("#### 📥 报告导出")
            st.download_button("一键下载综合分析与聚类报告 (.xlsx)", data=buffer.getvalue(), file_name="Antibody_Sequence_Analysis_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
        else:
            st.warning("⚠️ 未能在输入的文本中识别出标准抗体片段，请检查序列是否完整。")
    else:
        st.error("请先输入序列文本！")
