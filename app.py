# ... existing code ...
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
# ... existing code ...
def detect_ptms_detailed(seq, cdrs, domain_type):
# ... existing code ...
# ==========================================
# 3. 提取引擎 (修复降级 Bug + 高级正则增强)
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
            # 【核心升级】：提取底层引擎基于 IMGT 数据库比对出的种属和基因型
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
# ... existing code ...
def extract_vl_cdrs_regex(vl_seq):
# ... existing code ...
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
        # 传入 "VH" 标识
        cdrs = extract_cdrs_via_api(vh, "VH") if use_api else extract_vh_cdrs_regex(vh)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        # 【新增逻辑】：如果 API 返回了 IMGT 数据库识别的种属，则优先使用，否则降级使用我们手写的超强正则
        final_germline = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(vh)
        
        results.append({
            "序列名称 (ID)": seq_name, "区域": f"VH_{i+1}", "类型": "重链/纳米抗体",
            "同源 Germline": final_germline, "孤立Cys 雷达": detect_unpaired_cysteine(vh),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vh),
            "PTM 风险预警": detect_ptms_detailed(vh, cdrs, "VH"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vh
        })
    for i, vl in enumerate(vl_pattern.findall(clean_seq)):
        # 传入 "VL" 标识
        cdrs = extract_cdrs_via_api(vl, "VL") if use_api else extract_vl_cdrs_regex(vl)
        comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
        # 【新增逻辑】：轻链同样获取 IMGT 基因型
        final_germline = cdrs.get("API_Germline") if cdrs.get("API_Germline") else guess_germline(vl)
        
        results.append({
            "序列名称 (ID)": seq_name, "区域": f"VL_{i+1}", "类型": "轻链",
            "同源 Germline": final_germline, "孤立Cys 雷达": detect_unpaired_cysteine(vl),
            "CDR_GRAVY": calculate_gravy(comb_cdr), "pI (等电点)": calculate_pi(vl),
            "PTM 风险预警": detect_ptms_detailed(vl, cdrs, "VL"),
            "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vl
        })
    for i, fc in enumerate(fc_pattern.findall(clean_seq)):
# ... existing code ...
