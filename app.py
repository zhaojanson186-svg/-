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

st.title("🧬 高通量抗体序列处理与 CMC 质控中台 (V16 进阶版)")
st.info("💡 核心升级：新增 CDR3 核心指纹多样性分析与长度聚类引擎。完美兼容复杂鼠源抗体（Murine）及罕见框架突变的高鲁棒性解析。")

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
    """增强版 Germline 推断（支持小鼠、羊驼及多种人源化变种突变）"""
    seq = seq.upper()[:60] # 拉长扫描窗口
    # 纳米抗体
    if re.search(r'GGGSVQ', seq) or re.search(r'W[FY]RQAPGKERE', seq): return "Camelid VHH (纳米抗体)"
    
    # 鼠源特异性强探针 (Murine)
    if re.search(r'QQSG[AP]E[LV]V', seq) or re.search(r'QQSDA', seq) or re.search(r'GSLKLS', seq): return "Murine IGHV (鼠源)"
    if re.search(r'VK[IL]SC', seq) or re.search(r'V[KR]LSC', seq) or re.search(r'VTMSCK', seq): return "Murine IGHV (鼠源)"
    if re.search(r'GGLVKPGGSL', seq) or re.search(r'GPELVRPGAS', seq) or re.search(r'GPGILQPSQT', seq): return "Murine IGHV (鼠源)"
    
    if re.search(r'SP[SA]YLAASP', seq) or re.search(r'FMSTSVG', seq) or re.search(r'FMSTTIG', seq): return "Murine IGKV (鼠源)"
    if re.search(r'VTITCRAS', seq) or re.search(r'VSISCKAS', seq) or re.search(r'VTMTCSAS', seq) or re.search(r'VTMTC', seq): return "Murine IGKV (鼠源)"
    if re.search(r'QIVLTQSP', seq) or re.search(r'DIVMTQSQ', seq) or re.search(r'DIVMTQSP', seq): return "Murine IGKV (鼠源)"
    if re.search(r'QIVLSQSP', seq) or re.search(r'DIQMTQ', seq) or re.search(r'DTVLTQSP', seq): return "Murine IGKV (鼠源)"
    if re.search(r'DILMTQSP', seq) or re.search(r'DIVLTQSP', seq) or re.search(r'DIVITQSP', seq): return "Murine IGKV (鼠源)"

    # 高人源化探针 (Humanized)
    if re.search(r'SGGGLVQ', seq): return "Human IGHV3 (高人源化)"
    if re.search(r'SGAEVKKPG', seq): return "Human IGHV1/5 (高人源化)"
    if re.search(r'SGSELKKPG', seq): return "Human IGHV7 (高人源化)"
    if re.search(r'SGPGLVKPSG', seq): return "Human IGHV4 (高人源化)"
    if re.search(r'SGPEVKKPG', seq): return "Human IGHV2 (高人源化)"
    
    if re.search(r'SP[SS][SF]LSASVG', seq): return "Human IGKV1/3 (高人源化)"
    if re.search(r'SPLSLPVTPG', seq): return "Human IGKV2 (高人源化)"
    if re.search(r'SP[DS]SLA[VS]SLG', seq): return "Human IGKV4 (高人源化)"
    if re.search(r'QPPS[AS]SG', seq) or re.search(r'Q[PP]SVS[VAS]P', seq): return "Human IGLV (Lambda)"
    
    # 兜底匹配
    if re.search(r'VQL[VQE]QSG', seq) or re.search(r'VQL[LVE]ESG', seq): return "IGHV (亚族未定)"
    if re.search(r'[DE][IV][VQAM][ML]TQS', seq): return "IGKV (亚族未定)"
    if re.search(r'LTQP', seq): return "IGLV (Lambda未定)"
    
    if seq.startswith('Q') or seq.startswith('E') or seq.startswith('G') or seq.startswith('D') or seq.startswith('V'): return "疑似重链 (高变异)"
    if seq.startswith('D') or seq.startswith('A') or seq.startswith('Q'): return "疑似轻链 (高变异)"
    
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

# ==========================================
# 3. 提取引擎 (修复降级 Bug + 高级正则增强 + IMGT 对齐)
# ==========================================
def extract_cdrs_via_api(seq, chain_type="VH"):
    """【核心修复】：传入 chain_type，确保 API 失败时精准降级，绝不混用正则"""
    api_url = "https://api.antibody-informatics.org/v1/anarci/annotate"
    payload = {"sequence": seq, "scheme": "imgt"}
    try:
        # 添加微小延迟，防止并发过高被服务器掐断
        time.sleep(0.1)
        response = requests.post(api_url, json=payload, headers={"Content-Type": "application/json"}, timeout=5)
        if response.status_code == 200:
            data = response.json()
            # 提取底层引擎基于 IMGT 数据库比对出的种属和基因型
            species = str(data.get("species", "")).capitalize()
            v_gene = str(data.get("v_gene", ""))
            api_germline = f"{species} {v_gene}".strip() if (species or v_gene) else None
            
            return { 
                "CDR1": data.get("CDR1", "未识别"), 
                "CDR2": data.get("CDR2", "未识别"), 
                "CDR3": data.get("CDR3", "未识别"),
                "API_Germline": api_germline # 返回 IMGT 数据库识别结果
            }
        else: raise Exception("API 异常")
    except Exception:
        # 降级时不再猜首字母，直接根据传入的真实类型匹配
        if chain_type == "VH":
            return extract_vh_cdrs_regex(seq)
        else:
            return extract_vl_cdrs_regex(seq)

def extract_vh_cdrs_regex(vh_seq):
    """优化后的重链 CDR 正则提取引擎"""
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    
    # CDR3 提取 (兼容更多结尾变体)
    cdr3_match = re.search(r"Y[YFCAV][CST][A-Z]{1,3}(.*?)W[GS][A-Z][GTSVI]", vh_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    
    # CDR1 提取 (拓宽下游锚点，兼容 WVIQ 等罕见序列)
    cdr1_match = re.search(r"C[A-Z]{2,6}(.{5,16}?)W[VILFMA][A-Z]", vh_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    
    # CDR2 提取 (兼容 LDWIA 等异常突变的 FR2 结尾)
    cdr2_match = re.search(r"(?:[ELKDR][A-Z]{0,1}W[IVLMST][A-Z]{1,2}|REG[VLIA][A-Z]|RWV[A-Z])(.{8,30}?)[RKQ][VFSILAM][TVILAMFSC][A-Z]?", vh_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    return cdrs

def extract_vl_cdrs_regex(vl_seq):
    """优化后的轻链 CDR 正则提取引擎"""
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    
    # CDR3 提取 (兼容突变的 YVC 锚点和 FGGGT 等非标准 FR4 结尾)
    cdr3_match = re.search(r"Y[YFCAVS][CST](.*?)(?:F[GSA][A-Z]G|F[GSA][A-Z][GTV]|FGC)", vl_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    
    # CDR1 提取
    cdr1_match = re.search(r"C(.{8,18}?)W[YFL]", vl_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    
    # CDR2 提取
    cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq)
    if not cdr2_match: cdr2_match = re.search(r"W[YFL].{10,22}?([A-Z]{7})G[A-Z]{1,2}[RFS]", vl_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    return cdrs

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名待测序列_1"] = re.sub(r'\s+', '', text).upper()
        return sequences
    for part in text.split(">"):
        if not part.strip(): continue
        lines = part.strip().split("\n")
        name, seq = lines[0].strip(), "".join(lines[1:]).replace(" ", "").upper()
        if name and seq: sequences[name] = seq
    return sequences

# 增强型长序列捕获：兼容小鼠片段中特有的 G/D 头部起始，以及 VTVSA / LELK 等非经典收尾
vh_pattern = re.compile(r"([EQGVD].{90,140}(?:VTVSS|VTVSA))")
vl_pattern = re.compile(r"([DEQA].{90,130}(?:VEIK|LEIK|LELK|REIK|TVLG|VTVL|FGC))")
fc_pattern = re.compile(r"(CPPCP.*?LSPGK)")

def process_single_seq(seq_name, clean_seq, use_api):
    results = []
    
    for i, vh in enumerate(vh_pattern.findall(clean_seq)):
        cdrs = extract_cdrs_via_api(vh, "VH") if use_api else extract_vh_cdrs_regex(vh)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        final_germline = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(vh)
        
        results.append({
            "序列名称 (ID)": seq_name, "区域": f"VH_{i+1}", "类型": "重链/纳米抗体",
            "同源 Germline": final_germline, "孤立Cys 雷达": detect_unpaired_cysteine(vh),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vh),
            "PTM 风险预警": detect_ptms_detailed(vh, cdrs, "VH"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vh
        })
    for i, vl in enumerate(vl_pattern.findall(clean_seq)):
        cdrs = extract_cdrs_via_api(vl, "VL") if use_api else extract_vl_cdrs_regex(vl)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        final_germline = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(vl)
        
        results.append({
            "序列名称 (ID)": seq_name, "区域": f"VL_{i+1}", "类型": "轻链",
            "同源 Germline": final_germline, "孤立Cys 雷达": detect_unpaired_cysteine(vl),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vl),
            "PTM 风险预警": detect_ptms_detailed(vl, cdrs, "VL"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vl
        })
    for i, fc in enumerate(fc_pattern.findall(clean_seq)):
        results.append({
            "序列名称 (ID)": seq_name, "区域": f"Fc_{i+1}", "类型": "Fc区",
            "同源 Germline": "IgG Fc", "孤立Cys 雷达": "-", "CDR_GRAVY": "-",
            "pI (等电点)": calculate_pi(fc), "PTM 风险预警": detect_ptms_detailed(fc, {}, "Fc"),
            "CDR1": "-", "CDR2": "-", "CDR3": "-", "完整序列": fc
        })
    return results

# ==========================================
# 4. 模块一：序列预处理清洗器
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
# 5. 模块二：高通量解析与成药性打分
# ==========================================
st.markdown("---")
st.markdown("### 🔬 模块二：多线程序列解析与 CMC 排雷引擎")

col1, col2 = st.columns([3, 1])
with col1:
    raw_input = st.text_area("📥 请在此粘贴需要评估的候选序列 (最好使用带 > 的 FASTA 格式，以便识别配对关系):", height=250, key="main_input")
with col2:
    st.markdown("##### ⚙️ 引擎设置")
    engine_choice = st.radio("选择底层提取引擎：", 
                             ("🚀 本地正则引擎 (Regex)\n极速/零延迟/断网可用", 
                              "☁️ 外部云端对齐 (API)\n融合 IMGT 智能判断"))
    
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
            
            # --- 智能名称清洗：自动剥离重轻链后缀，提取核心分子ID进行配对 ---
            def extract_base_name(name):
                # 正则匹配末尾的 _VH, -VL, _Heavy_chain 等标识并剔除
                return re.sub(r'[-_](?i)(VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light)$', '', name).strip()
            df['归属分子名'] = df['序列名称 (ID)'].apply(extract_base_name)
            
            # --- 数据高亮与美化 ---
            def highlight_alerts(val):
                val_str = str(val)
                if '🚨' in val_str or '高危' in val_str or '脱氨基' in val_str or '异构化' in val_str or '重复' in val_str:
                    return 'background-color: #ffcccc; color: #990000; font-weight: bold'
                elif '⚠️' in val_str or '氧化' in val_str or '冗余' in val_str:
                    return 'background-color: #fff2cc; color: #b26b00'
                return ''
            
            st.markdown("#### 📊 1. CMC 序列解析总表")
            st.dataframe(df.drop(columns=['归属分子名']).style.map(highlight_alerts, subset=['孤立Cys 雷达', 'PTM 风险预警']), use_container_width=True)
            
            # ==========================================
            # 新增核心功能：CDR3 核心指纹与多样性聚类
            # ==========================================
            st.markdown("#### 🧬 2. CDR3 核心指纹与多样性聚类 (Unique Analysis)")
            st.caption("抗原识别的核心在于 CDR3。本模块已对高置信度 CDR3 进行提纯去重，分析其出现频次与长度特征，帮助锁定核心克隆簇 (Clonal Families)。")
            
            df_v = df[df['类型'].isin(['重链/纳米抗体', '轻链'])]
            # 过滤掉未识别的无效数据，并创建副本以避免 SettingWithCopyWarning
            df_valid_cdr3 = df_v[~df_v['CDR3'].str.contains('未识别|失败', na=False)].copy()
            
            if not df_valid_cdr3.empty:
                # 计算每个提取出的 CDR3 长度
                df_valid_cdr3['CDR3长度 (AA)'] = df_valid_cdr3['CDR3'].apply(len)
                
                # 基于 CDR3 序列进行高阶聚类分析
                cluster_cdr3 = df_valid_cdr3.groupby(['类型', 'CDR3']).agg(
                    出现频次=('序列名称 (ID)', 'count'),
                    CDR3长度=('CDR3长度 (AA)', 'first'),
                    同源Germline=('同源 Germline', 'first'),
                    来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x.unique())),
                    代表完整序列=('完整序列', 'first')
                ).reset_index().sort_values(by=['类型', '出现频次'], ascending=[True, False])
                
                # 在前端将重链和轻链分开展示，UI更清晰
                c1, c2 = st.columns(2)
                with c1:
                    st.markdown("##### 🔵 重链 CDR3 多样性")
                    vh_cdr3 = cluster_cdr3[cluster_cdr3['类型'].str.contains('重')].drop(columns=['类型'])
                    st.dataframe(vh_cdr3, use_container_width=True, hide_index=True)
                with c2:
                    st.markdown("##### 🟢 轻链 CDR3 多样性")
                    vl_cdr3 = cluster_cdr3[cluster_cdr3['类型'].str.contains('轻')].drop(columns=['类型'])
                    st.dataframe(vl_cdr3, use_container_width=True, hide_index=True)
            else:
                st.info("尚未在当前数据中提取到有效的 CDR3 指纹。")

            # --- 全序列级别聚类 ---
            cluster_v = df_v.groupby('完整序列').agg(
                链类型=('类型', 'first'), 相同序列数=('序列名称 (ID)', 'count'),
                来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x.unique())),
                PTM风险=('PTM 风险预警', 'first'), CDR3=('CDR3', 'first')
            ).reset_index().sort_values(by=['链类型', '相同序列数'], ascending=[True, False])

            # --- 核心升级：完整抗体分子 (Fv) 智能组装与唯一性 (Unique) 判断 ---
            paired_data = []
            grouped = df.groupby('归属分子名')
            for name, group in grouped:
                vh_rows = group[group['类型'].str.contains('重链')]
                vl_rows = group[group['类型'].str.contains('轻链')]
                
                # 当【归属同一个分子名】下同时存在 VH 和 VL 时，进行双链组装
                if not vh_rows.empty and not vl_rows.empty:
                    vh_pi = vh_rows.iloc[0]['pI (等电点)']
                    vl_pi = vl_rows.iloc[0]['pI (等电点)']
                    vh_seq = vh_rows.iloc[0]['完整序列']
                    vl_seq = vl_rows.iloc[0]['完整序列']
                    vh_cdr3 = vh_rows.iloc[0]['CDR3']
                    vl_cdr3 = vl_rows.iloc[0]['CDR3']
                    
                    delta_pi = abs(vh_pi - vl_pi)
                    warning_flag = "🚨 高危: ΔpI > 2.0 (易发生配对错乱)" if delta_pi > 2.0 else "✅ 正常 (电荷分布对称)"
                    
                    paired_data.append({
                        "核心分子名 (ID)": name,
                        "重链 (VH) pI": vh_pi,
                        "轻链 (VL) pI": vl_pi,
                        "ΔpI (电荷不对称性)": round(delta_pi, 2),
                        "Fv 链间质控状态": warning_flag,
                        "VH_完整序列": vh_seq,
                        "VL_完整序列": vl_seq,
                        "VH_CDR3": vh_cdr3,
                        "VL_CDR3": vl_cdr3
                    })
            
            df_paired = pd.DataFrame(paired_data)
            
            if not df_paired.empty:
                # 拼接 VH 和 VL 序列生成该双链分子的“绝对序列指纹”
                df_paired['Fv_指纹'] = df_paired['VH_完整序列'] + "||" + df_paired['VL_完整序列']
                
                # 依据分子指纹进行聚类计算，找出 Unique Fv 分子
                cluster_fv = df_paired.groupby('Fv_指纹').agg(
                    出现频次=('核心分子名 (ID)', 'count'),
                    是否唯一=('核心分子名 (ID)', lambda x: "✅ 唯一分子 (Unique)" if len(x) == 1 else "⚠️ 冗余克隆 (Duplicate)"),
                    同源分子重叠名单=('核心分子名 (ID)', lambda x: ', '.join(x.unique())),
                    重链_CDR3=('VH_CDR3', 'first'),
                    轻链_CDR3=('VL_CDR3', 'first'),
                    ΔpI_电荷不对称性=('ΔpI (电荷不对称性)', 'first'),
                    成药性质控=('Fv 链间质控状态', 'first')
                ).reset_index().drop(columns=['Fv_指纹']).sort_values(by=['出现频次', '是否唯一'], ascending=[False, True])
                
                st.markdown("#### 🧩 3. 完整抗体双链 (Fv) 智能组装与唯一性聚类")
                st.caption("系统已自动剔除前缀/后缀，将同属一个分子的 VH 和 VL 序列进行组装。下表直接为您提取了 **Unique 完整克隆库** 及其对应的双链配对成药性风险。")
                st.dataframe(cluster_fv.style.map(highlight_alerts, subset=['成药性质控', '是否唯一']), use_container_width=True, hide_index=True)

            else:
                st.info("💡 完整分子组装雷达待命：本次解析未检测到配对的重/轻链。")

            # --- 导出模块 ---
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.drop(columns=['归属分子名']).to_excel(writer, index=False, sheet_name='1_单链_成药性解析总表')
                if not df_valid_cdr3.empty: cluster_cdr3.to_excel(writer, index=False, sheet_name='2_CDR3_多样性与长度聚类')
                if not df_paired.empty: cluster_fv.to_excel(writer, index=False, sheet_name='3_完整分子(Fv)_唯一性分析')
                if not cluster_v.empty: cluster_v.to_excel(writer, index=False, sheet_name='4_单链_全序列同源去重')
            
            st.markdown("#### 📥 4. 报告导出")
            st.download_button("一键下载综合分析与聚类报告 (.xlsx)", data=buffer.getvalue(), file_name="Antibody_Sequence_Analysis_Report_V16.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
        else:
            st.warning("⚠️ 未能在输入的文本中识别出标准抗体片段，请检查序列是否完整。")
    else:
        st.error("请先输入序列文本！")
