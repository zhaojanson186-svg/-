import streamlit as st
import re
import pandas as pd
import io
import requests
import concurrent.futures
import time

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="抗体核心分析中台 V16", page_icon="🧬", layout="wide")

st.title("🧬 高通量抗体序列处理与 CMC 质控中台 (V16 最终版)")
st.info("💡 更新日志：PTM 引擎已聚焦工业级致命缺陷；重构了同名序列防覆盖解析系统；新增 Fv 组合唯一性 (Unique Clone) 聚类排重功能。")

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
    seq = seq.upper()[:60]
    if re.search(r'GGGSVQ', seq) or re.search(r'W[FY]RQAPGKERE', seq): return "Camelid VHH (纳米抗体)"
    if re.search(r'SGGGLVQ', seq): return "Human IGHV3 (高人源化)"
    if re.search(r'SGAEVKKPG', seq): return "Human IGHV1/5 (高人源化)"
    if re.search(r'SGSELKKPG', seq): return "Human IGHV7 (高人源化)"
    if re.search(r'SGPGLVKPSG', seq): return "Human IGHV4 (高人源化)"
    if re.search(r'SGPEVKKPG', seq): return "Human IGHV2 (高人源化)"
    if re.search(r'QQSG[AP]E[LV]V', seq) or re.search(r'QQSDA', seq) or re.search(r'GSLKLS', seq): return "Murine IGHV (鼠源)"
    if re.search(r'VK[IL]SC', seq) or re.search(r'V[KR]LSC', seq) or re.search(r'VTMSCK', seq): return "Murine IGHV (鼠源)"
    if re.search(r'VQL[VQE]QSG', seq) or re.search(r'VQL[LVE]ESG', seq): return "IGHV (亚族未定)"
    if re.search(r'SP[SS][SF]LSASVG', seq): return "Human IGKV1/3 (高人源化)"
    if re.search(r'SPLSLPVTPG', seq): return "Human IGKV2 (高人源化)"
    if re.search(r'SP[DS]SLA[VS]SLG', seq): return "Human IGKV4 (高人源化)"
    if re.search(r'SP[SA]YLAASP', seq) or re.search(r'FMSTSVG', seq): return "Murine IGKV (鼠源)"
    if re.search(r'VTITCRAS', seq) or re.search(r'VSISCKAS', seq) or re.search(r'VTMTCSAS', seq): return "Murine IGKV (鼠源)"
    if re.search(r'[DE][IV][VQAM][ML]TQS', seq): return "IGKV (亚族未定)"
    if re.search(r'QPPS[AS]SG', seq) or re.search(r'Q[PP]SVS[VAS]P', seq): return "Human IGLV (Lambda)"
    if re.search(r'LTQP', seq): return "IGLV (Lambda未定)"
    if seq.startswith('Q') or seq.startswith('E') or seq.startswith('G') or seq.startswith('D'): return "疑似重链 (高度变异)"
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
    """【PTM 核心瘦身】：只提示发生在 CDR 上的，且具有致命成药性缺陷的位点"""
    if "Fc" in domain_type: return "✅ 无 CDR PTM (Fc区)"
    
    region_finder = get_region_finder(seq, cdrs, domain_type)
    
    # 极简致命规则：去除了一般的 M 氧化和慢速脱氨/异构化，只抓大雷
    ptm_rules = {
        "N-糖基化": r"N[^P][ST]", 
        "极速脱氨基(NG)": r"NG", 
        "极速异构化(DG)": r"DG", 
        "酸断裂点(DP)": r"DP"
    }
    
    found_ptms = []
    for ptm_name, pattern in ptm_rules.items():
        for match in re.finditer(pattern, seq):
            region_name = region_finder(match.start())
            # 严格过滤：只收录位于 CDR 区域内的修饰
            if region_name.startswith("CDR"):
                found_ptms.append(f"[{region_name}] {ptm_name}({match.group()}) @{match.start()+1}")
                
    found_ptms.sort(key=lambda x: int(re.search(r'@(\d+)', x).group(1)) if re.search(r'@(\d+)', x) else 0)
    return " | ".join(found_ptms) if found_ptms else "✅ 无高危 PTM"

# ==========================================
# 3. 提取引擎 (Cys 锚定精准切片)
# ==========================================
def extract_cdrs_via_api(seq, chain_type="VH"):
    api_url = "https://api.antibody-informatics.org/v1/anarci/annotate"
    payload = {"sequence": seq, "scheme": "imgt"}
    try:
        time.sleep(0.1)
        response = requests.post(api_url, json=payload, headers={"Content-Type": "application/json"}, timeout=5)
        if response.status_code == 200:
            data = response.json()
            api_germline = f"{data.get('Species', '')} {data.get('V_gene', '')}".strip()
            return { 
                "CDR1": data.get("CDR1", "未识别"), 
                "CDR2": data.get("CDR2", "未识别"), 
                "CDR3": data.get("CDR3", "未识别"),
                "API_Germline": api_germline if api_germline else None
            }
        else: raise Exception("API 异常")
    except Exception:
        if chain_type == "VH": return extract_vh_cdrs_regex(seq)
        else: return extract_vl_cdrs_regex(seq)

def extract_vh_cdrs_regex(vh_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别", "API_Germline": None}
    # Cys 锚定防吞噬，限定 4-30 氨基酸长度
    cdr3_match = re.search(r"C([A-Z]{4,30})(?=W[GS][A-Z][GTSVI]|WG[QA]G|W[GS]G)", vh_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    
    cdr1_match = re.search(r"C[A-Z]{2,6}(.{5,16}?)W[VILFMA][A-Z]", vh_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    
    cdr2_match = re.search(r"(?:[ELKDR][A-Z]{0,1}W[IVLMST][A-Z]{1,2}|REG[VLIA][A-Z]|RWV[A-Z])(.{8,30}?)[RKQ][VFSILAM][TVILAMFSC][A-Z]?", vh_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    return cdrs

def extract_vl_cdrs_regex(vl_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别", "API_Germline": None}
    # Cys 锚定防吞噬，限定 4-25 氨基酸长度
    cdr3_match = re.search(r"C([A-Z]{4,25})(?=F[GSA][A-Z]G|F[GSA][A-Z][GTV]|FGC)", vl_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    
    cdr1_match = re.search(r"C(.{8,18}?)W[YFL]", vl_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    
    cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq)
    if not cdr2_match: cdr2_match = re.search(r"W[YFL].{10,22}?([A-Z]{7})G[A-Z]{1,2}[RFS]", vl_seq)
    if cdr2_match: cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq).group(1) if re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq) else cdr2_match.group(1)
    return cdrs

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名待测序列_1"] = re.sub(r'\s+', '', text).upper()
        return sequences
        
    name_counter = {}
    for part in text.split(">"):
        if not part.strip(): continue
        lines = part.strip().split("\n")
        raw_name = lines[0].strip()
        seq = "".join(lines[1:]).replace(" ", "").upper()
        if raw_name and seq:
            if raw_name in name_counter:
                name_counter[raw_name] += 1
                # 如果名称带有 _VH/_VL 等后缀，确保在重命名时后缀依然在最后，方便后续配对
                match = re.search(r'([-_](VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light))$', raw_name, flags=re.IGNORECASE)
                if match:
                    suffix = match.group(1)
                    base = raw_name[:-len(suffix)]
                    name = f"{base}_Dup{name_counter[raw_name]}{suffix}"
                else:
                    name = f"{raw_name}_Dup{name_counter[raw_name]}"
            else:
                name_counter[raw_name] = 1
                name = raw_name
                
            sequences[name] = seq
    return sequences

def process_single_seq(seq_name, clean_seq, use_api):
    results = []
    upper_name = seq_name.upper()
    is_vh = any(k in upper_name for k in ['VH', 'HC', 'HEAVY', '重链'])
    is_vl = any(k in upper_name for k in ['VL', 'LC', 'LIGHT', '轻链'])
    is_fc = 'FC' in upper_name
    
    if is_vh:
        cdrs = extract_cdrs_via_api(clean_seq, "VH") if use_api else extract_vh_cdrs_regex(clean_seq)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        germline_val = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(clean_seq)
        results.append({
            "序列名称 (ID)": seq_name, "区域": "VH_1", "类型": "重链/纳米抗体",
            "同源 Germline": germline_val, "孤立Cys 雷达": detect_unpaired_cysteine(clean_seq),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(clean_seq),
            "PTM 风险预警": detect_ptms_detailed(clean_seq, cdrs, "VH"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": clean_seq
        })
    elif is_vl:
        cdrs = extract_cdrs_via_api(clean_seq, "VL") if use_api else extract_vl_cdrs_regex(clean_seq)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        germline_val = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(clean_seq)
        results.append({
            "序列名称 (ID)": seq_name, "区域": "VL_1", "类型": "轻链",
            "同源 Germline": germline_val, "孤立Cys 雷达": detect_unpaired_cysteine(clean_seq),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(clean_seq),
            "PTM 风险预警": detect_ptms_detailed(clean_seq, cdrs, "VL"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": clean_seq
        })
    elif is_fc:
        results.append({
            "序列名称 (ID)": seq_name, "区域": "Fc_1", "类型": "Fc区",
            "同源 Germline": "IgG Fc", "孤立Cys 雷达": "-", "CDR_GRAVY": "-",
            "pI (等电点)": calculate_pi(clean_seq), "PTM 风险预警": detect_ptms_detailed(clean_seq, {}, "Fc"),
            "CDR1": "-", "CDR2": "-", "CDR3": "-", "完整序列": clean_seq
        })
    else:
        cdrs = extract_cdrs_via_api(clean_seq, "VH") if use_api else extract_vh_cdrs_regex(clean_seq)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        germline_val = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(clean_seq)
        results.append({
            "序列名称 (ID)": seq_name, "区域": "VH_1", "类型": "重链/纳米抗体",
            "同源 Germline": germline_val, "孤立Cys 雷达": detect_unpaired_cysteine(clean_seq),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(clean_seq),
            "PTM 风险预警": detect_ptms_detailed(clean_seq, cdrs, "VH"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": clean_seq
        })
    return results

# ==========================================
# 4. 模块一：序列预处理清洗器
# ==========================================
st.markdown("---")
st.markdown("### 🗂️ 模块一：序列预处理清洗器")
uploaded_excel = st.file_uploader("📤 请上传包含重/轻链序列的 Excel 或 CSV 文件", type=['xlsx', 'xls', 'csv'], key="excel_uploader")
if uploaded_excel is not None:
    if st.button("🔄 一键清洗并转换为 FASTA 格式", type="primary"):
        try:
            if uploaded_excel.name.endswith('.csv'): df_raw = pd.read_csv(uploaded_excel, header=None)
            else: df_raw = pd.read_excel(uploaded_excel, header=None)
            fasta_lines = []
            invalid_names = ['nan', 'None', '', 'Protein name', 'Sequence']
            for index, row in df_raw.iterrows():
                name = str(row.iloc[0]).strip() if len(row) > 0 else ""
                hc = str(row.iloc[1]).strip() if len(row) > 1 else ""
                lc = str(row.iloc[2]).strip() if len(row) > 2 else ""
                if name in invalid_names or "sequence" in name.lower(): continue
                if hc and hc not in ['nan', 'None'] and len(hc) > 10: fasta_lines.append(f">{name}_VH\n{hc.upper()}")
                if lc and lc not in ['nan', 'None'] and len(lc) > 10: fasta_lines.append(f">{name}_VL\n{lc.upper()}")
            fasta_content = '\n'.join(fasta_lines)
            if fasta_content:
                st.success(f"✅ 成功清洗出 {fasta_content.count('>')} 条链序列。")
                st.text_area("📋 提取出的标准化 FASTA 序列:", value=fasta_content, height=200)
            else:
                st.warning("⚠️ 未能提取出有效序列。")
        except Exception as e:
            st.error(f"❌ 读取错误: {e}")

# ==========================================
# 5. 模块二：高通量解析与成药性打分
# ==========================================
st.markdown("---")
st.markdown("### 🔬 模块二：多线程序列解析与 CMC 排雷引擎")

col1, col2 = st.columns([3, 1])
with col1:
    raw_input = st.text_area("📥 请在此粘贴需要评估的候选序列 (支持 FASTA 格式):", height=250, key="main_input")
with col2:
    st.markdown("##### ⚙️ 引擎设置")
    engine_choice = st.radio("选择底层提取引擎：", ("🚀 本地正则引擎 (Regex)\n极速/推荐", "☁️ 外部云端对齐 (API)\n高精度"))
    analyze_btn = st.button("🚀 启动并发解析", type="primary", use_container_width=True)

if analyze_btn:
    use_api = "API" in engine_choice
    if raw_input:
        seq_dict = parse_fasta(raw_input)
        total_seqs = len(seq_dict)
        st.success(f"✅ 成功读取 {total_seqs} 条输入序列，正在启动并发扫库...")
        
        all_results = []
        progress_bar = st.progress(0)
        workers = 3 if use_api else 8
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_seq = {executor.submit(process_single_seq, name, seq, use_api): name for name, seq in seq_dict.items()}
            completed = 0
            for future in concurrent.futures.as_completed(future_to_seq):
                try:
                    res = future.result()
                    all_results.extend(res)
                except Exception as e:
                    st.error(f"序列解析失败: {e}")
                completed += 1
                progress_bar.progress(completed / total_seqs)
        
        if all_results:
            df = pd.DataFrame(all_results)
            
            def highlight_alerts(val):
                val_str = str(val)
                # 针对极其致命缺陷高红警告
                if '🚨' in val_str or '糖基化' in val_str or 'NG' in val_str or 'DG' in val_str or 'DP' in val_str or '冗余' in val_str:
                    return 'background-color: #ffcccc; color: #990000; font-weight: bold'
                # 针对较弱的缺陷给予黄底警告
                elif '⚠️' in val_str:
                    return 'background-color: #fff2cc; color: #b26b00'
                return ''
            
            st.markdown("#### 📊 CMC 序列解析总表")
            st.dataframe(df.style.map(highlight_alerts, subset=['孤立Cys 雷达', 'PTM 风险预警']), use_container_width=True)
            
            # 使用标准的 flags 参数，避免正则表达式模式错误
            def extract_base_name(name):
                return re.sub(r'[-_](VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light)$', '', name, flags=re.IGNORECASE).strip()
            
            df['归属分子名'] = df['序列名称 (ID)'].apply(extract_base_name)
            
            paired_data = []
            for name, group in df.groupby('归属分子名'):
                vh_rows = group[group['类型'].str.contains('重链')]
                vl_rows = group[group['类型'].str.contains('轻链')]
                
                if not vh_rows.empty and not vl_rows.empty:
                    vh_seq = vh_rows.iloc[0]['完整序列']
                    vl_seq = vl_rows.iloc[0]['完整序列']
                    vh_pi = vh_rows.iloc[0]['pI (等电点)']
                    vl_pi = vl_rows.iloc[0]['pI (等电点)']
                    delta_pi = abs(vh_pi - vl_pi)
                    
                    warning_flag = "⚠️ 关注: ΔpI > 2.0" if delta_pi > 2.0 else "✅ 正常 (电荷分布对称)"
                    
                    vh_ptm = vh_rows.iloc[0]['PTM 风险预警']
                    vl_ptm = vl_rows.iloc[0]['PTM 风险预警']
                    combined_ptm = []
                    if vh_ptm and "✅" not in vh_ptm:
                        combined_ptm.append(f"VH: {vh_ptm}")
                    if vl_ptm and "✅" not in vl_ptm:
                        combined_ptm.append(f"VL: {vl_ptm}")
                    ptm_summary = " | ".join(combined_ptm) if combined_ptm else "✅ Fv 无高危 PTM"

                    paired_data.append({
                        "核心分子名 (ID)": name, 
                        "VH_Seq": vh_seq, "VL_Seq": vl_seq,
                        "重链_pI": vh_pi, "轻链_pI": vl_pi,
                        "ΔpI": round(delta_pi, 2), "Fv质控状态": warning_flag,
                        "PTM风险汇总": ptm_summary
                    })
            
            df_paired = pd.DataFrame(paired_data)
            if not df_paired.empty:
                # 核心升级：生成 VH + VL 的双链联合指纹，进行全局去重与 Unique 判定
                df_paired['Fv_Fingerprint'] = df_paired['VH_Seq'] + "||" + df_paired['VL_Seq']
                
                fv_cluster = df_paired.groupby('Fv_Fingerprint').agg(
                    包含相同配对数=('核心分子名 (ID)', 'count'),
                    代表分子名=('核心分子名 (ID)', 'first'),
                    合并来源分子名=('核心分子名 (ID)', lambda x: ', '.join(x.unique())),
                    重链_pI=('重链_pI', 'first'),
                    轻链_pI=('轻链_pI', 'first'),
                    ΔpI=('ΔpI', 'first'),
                    Fv质控状态=('Fv质控状态', 'first'),
                    PTM风险汇总=('PTM风险汇总', 'first')
                ).reset_index()
                
                # 基于频次打标签
                fv_cluster['唯一性 (Unique)'] = fv_cluster['包含相同配对数'].apply(
                    lambda x: "✅ 唯一 (Unique)" if x == 1 else f"🚨 冗余克隆 (Dup x{x})"
                )
                
                # 重新整理列名供前端显示，先排序再切片，避免 KeyError
                df_paired_final = fv_cluster.sort_values(
                    by=['包含相同配对数', '代表分子名'], ascending=[False, True]
                )[[
                    '唯一性 (Unique)', '包含相同配对数', '代表分子名', '合并来源分子名', 
                    '重链_pI', '轻链_pI', 'ΔpI', 'Fv质控状态', 'PTM风险汇总'
                ]]
                
                # 更新到 df_paired 供后续导出下载使用
                df_paired = df_paired_final
                
                st.markdown("#### ⚖️ 分子水平 Fv 质控与配对唯一性 (Unique Fv)")
                st.dataframe(df_paired_final.style.map(highlight_alerts, subset=['Fv质控状态', 'PTM风险汇总', '唯一性 (Unique)']), use_container_width=True)
            
            df_v = df[df['类型'].isin(['重链/纳米抗体', '轻链'])]
            df_valid_cdr3 = df_v[~df_v['CDR3'].str.contains('未识别|失败', na=False)].copy()
            cluster_cdr3 = pd.DataFrame()
            if not df_valid_cdr3.empty:
                df_valid_cdr3['CDR3长度 (AA)'] = df_valid_cdr3['CDR3'].apply(len)
                cluster_cdr3 = df_valid_cdr3.groupby(['类型', 'CDR3']).agg(
                    出现频次=('序列名称 (ID)', 'count'), CDR3长度=('CDR3长度 (AA)', 'first'),
                    同源Germline=('同源 Germline', 'first'), 来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x.unique()))
                ).reset_index().sort_values(by=['类型', '出现频次'], ascending=[True, False])
                
                st.markdown("#### 🧬 CDR3 核心指纹与多样性统计")
                st.dataframe(cluster_cdr3, use_container_width=True, hide_index=True)

            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.drop(columns=['归属分子名']).to_excel(writer, index=False, sheet_name='单链_解析总表')
                if not df_paired.empty: df_paired.to_excel(writer, index=False, sheet_name='配对_Fv双链汇总')
                if not cluster_cdr3.empty: cluster_cdr3.to_excel(writer, index=False, sheet_name='CDR3_多样性统计')
            
            st.markdown("#### 📥 报告导出")
            st.download_button("一键下载综合分析与聚类报告 (.xlsx)", data=buffer.getvalue(), file_name="Antibody_Sequence_Analysis_Report_V16_Final.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
        else:
            st.warning("⚠️ 未能识别出有效片段。")
    else:
        st.error("请先输入序列文本！")
