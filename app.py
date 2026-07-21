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
            # 【新增】：提取底层引擎基于 IMGT 数据库比对出的种属和基因型
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
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别", "API_Germline": None}
# ... existing code ...
    cdr2_match = re.search(r"(?:[ELKDR][A-Z]{0,1}W[IVLMST][A-Z]{1,2}|REG[VLIA][A-Z]|RWV[A-Z])(.{8,30}?)[RKQ][VFSILAM][TVILAMFSC][A-Z]?", vh_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    return cdrs

def extract_vl_cdrs_regex(vl_seq):
    """优化后的轻链 CDR 正则提取引擎"""
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别", "API_Germline": None}
# ... existing code ...
    cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq)
    if not cdr2_match: cdr2_match = re.search(r"W[YFL].{10,22}?([A-Z]{7})G[A-Z]{1,2}[RFS]", vl_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    return cdrs

def parse_fasta(text):
# ... existing code ...
vl_pattern = re.compile(r"([DEQA].{95,125}(?:VEIK|LEIK|TVLG|VTVL|FGC))")
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
        # 【新增逻辑】
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
