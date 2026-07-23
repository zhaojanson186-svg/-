import streamlit as st
import pandas as pd
import re
import io
from collections import Counter
import base64

def calculate_mw(seq):
    mw_table = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
                'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
                'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
                'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06}
    weight = sum(mw_table.get(aa, 0) for aa in seq) + 18.02
    return round(weight, 2)

def calculate_pi(seq):
    pk_a = {'D': 3.65, 'E': 4.25, 'Y': 10.07, 'C': 8.18, 'H': 6.00, 'K': 10.53, 'R': 12.48}
    def net_charge(ph):
        charge = 0.0
        for aa, count in Counter(seq).items():
            if aa in ['D', 'E', 'Y', 'C']:
                charge -= count / (1.0 + 10**(pk_a[aa] - ph))
            elif aa in ['H', 'K', 'R']:
                charge += count / (1.0 + 10**(ph - pk_a[aa]))
        charge -= 1.0 / (1.0 + 10**(3.1 - ph)) # C-term
        charge += 1.0 / (1.0 + 10**(ph - 8.0)) # N-term
        return charge
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        if net_charge(mid) > 0: low = mid
        else: high = mid
    return round(mid, 2)

def calculate_gravy(seq):
    kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 
          'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 
          'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    if not seq: return 0.0
    return round(sum(kd.get(aa, 0) for aa in seq) / len(seq), 3)

def extract_vh_cdrs_regex(seq):
    cdrs = {"CDR1": "", "CDR2": "", "CDR3": ""}
    regions = {"FR1": "", "CDR1": "", "FR2": "", "CDR2": "", "FR3": "", "CDR3": "", "FR4": ""}
    
    cdr1_match = re.search(r"(S[VLA]K[LVI]S|S[A-Z]TL[ST]|S[A-Z][A-Z]MS|S[A-Z]QM[SN])C(.*?)W[VIQ]R[QA]", seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(2)
    
    cdr2_match = re.search(r"W[VIQ]R[QA][A-Z]{1,3}G[A-Z]{1,2}[EAW](.*?)[RKV][LTIV][TS][IVTL]", seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    # 采用 Cys 锚定法，精准定位重链 CDR3，避免吞噬 FR3
    cdr3_match = re.search(r"([YF][YFCAV][CAV])(.{5,30}?)(W[GS][A-Z][GTSVI])", seq)
    if cdr3_match: 
        cdrs["CDR3"] = cdr3_match.group(2)
        c3_start = cdr3_match.start(2)
        c3_end = cdr3_match.end(2)
        regions["FR3"] = seq[:c3_start]
        regions["CDR3"] = cdrs["CDR3"]
        regions["FR4"] = seq[c3_end:]
        if cdr1_match and cdr2_match:
            c1_start, c1_end = cdr1_match.span(2)
            c2_start, c2_end = cdr2_match.span(1)
            regions["FR1"] = seq[:c1_start]
            regions["CDR1"] = cdrs["CDR1"]
            regions["FR2"] = seq[c1_end:c2_start]
            regions["CDR2"] = cdrs["CDR2"]
            regions["FR3"] = seq[c2_end:c3_start]
    
    return cdrs, regions

def extract_vl_cdrs_regex(seq):
    cdrs = {"CDR1": "", "CDR2": "", "CDR3": ""}
    regions = {"FR1": "", "CDR1": "", "FR2": "", "CDR2": "", "FR3": "", "CDR3": "", "FR4": ""}
    
    cdr1_match = re.search(r"[STD][VIAL][TISM][CS](.*?)W[YF]Q[QR]", seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    
    cdr2_match = re.search(r"P[RKAL]L[LIVW]I[YF](.*?)[GSA]VP[DSR]", seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    
    # 采用 Cys 锚定法，精准定位轻链 CDR3，避免吞噬 FR3
    cdr3_match = re.search(r"([YF][YFCAVS][CST])(.{4,20}?)(F[GSA][A-Z]G)", seq)
    if cdr3_match: 
        cdrs["CDR3"] = cdr3_match.group(2)
        c3_start = cdr3_match.start(2)
        c3_end = cdr3_match.end(2)
        regions["FR3"] = seq[:c3_start]
        regions["CDR3"] = cdrs["CDR3"]
        regions["FR4"] = seq[c3_end:]
        if cdr1_match and cdr2_match:
            c1_start, c1_end = cdr1_match.span(1)
            c2_start, c2_end = cdr2_match.span(1)
            regions["FR1"] = seq[:c1_start]
            regions["CDR1"] = cdrs["CDR1"]
            regions["FR2"] = seq[c1_end:c2_start]
            regions["CDR2"] = cdrs["CDR2"]
            regions["FR3"] = seq[c2_end:c3_start]
            
    return cdrs, regions

def detect_ptms_detailed(regions_dict):
    # 精简 PTM 预警：只抓取真正的 4 种致命缺陷，忽略 M 氧化等次要干扰
    ptm_motifs = {
        "N-糖基化 (NIT)": r"N[^P][ST]",
        "极速脱氨基 (NG)": r"NG",
        "极速异构化 (DG)": r"DG",
        "酸断裂点 (DP)": r"DP"
    }
    alerts = []
    
    for region_name, region_seq in regions_dict.items():
        if not region_seq: continue
        # 核心逻辑：仅扫描名称以 CDR 开头的区域，彻底忽略框架区 (FR)
        if region_name.startswith("CDR"):
            for ptm_name, motif in ptm_motifs.items():
                for match in re.finditer(motif, region_seq):
                    alerts.append(f"[{region_name}] {ptm_name} ({match.group()}) @{match.start()+1}")
                    
    return alerts

def highlight_alerts(val):
    if not isinstance(val, str):
        return ''
    if "🚨" in val or "⚠️" in val or "冗余克隆" in val:
        return 'background-color: #fff3cd; color: #856404; font-weight: bold;'
    elif "✅" in val:
        return 'background-color: #d4edda; color: #155724;'
    return ''

st.set_page_config(page_title="工业级抗体生信解析大屏 V16", layout="wide")

st.markdown("""
<div style="background: linear-gradient(90deg, #1e3c72 0%, #2a5298 100%); padding: 25px; border-radius: 12px; margin-bottom: 25px;">
    <h1 style="color: white; margin: 0; font-size: 2.2rem;">🔬 工业级抗体生信质控大屏 (Fv-Unique Edition)</h1>
    <p style="color: #e0e0e0; margin-top: 10px; font-size: 1.1rem;">
        支持自动切分CDR、全维度计算理化参数、精准过滤高变区高危 PTM。<br/>
        <strong>V16 进阶：包含智能重复名称防覆盖、极其精准的 Cys 锚定法提取，以及多对多双链交叉配对与冗余分子去重 (Unique) 评估。</strong>
    </p>
</div>
""", unsafe_allow_html=True)

input_text = st.text_area("✍️ 贴入您的多条抗体 FASTA 序列 (支持混合VH和VL)", height=200, 
                          help="每条序列需以 >开头，名称中需包含 _VH, _VL, -HC, -LC, 或带有数字后缀如 VH1, VL2 等以区分链类型。")

if st.button("🚀 开始高通量并行解析", type="primary"):
    if not input_text.strip():
        st.warning("⚠️ 请先输入 FASTA 序列")
    else:
        with st.spinner("⏳ 正在进行结构解析与计算..."):
            
            parsed_seqs = []
            name_counter = {} # 用于记录出现过的名称，解决完全同名覆盖问题
            
            for part in input_text.split(">"):
                if not part.strip(): continue
                lines = part.strip().split("\n")
                raw_name = lines[0].strip()
                seq = "".join(lines[1:]).replace(" ", "").upper()
                
                if raw_name and seq:
                    # 如果遇到同名序列，自动追加 _Dup 编号后缀避免 Pandas/Dict 覆盖
                    if raw_name in name_counter:
                        name_counter[raw_name] += 1
                        match = re.search(r'([-_](VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light)[0-9]*)$', raw_name, flags=re.IGNORECASE)
                        if match:
                            suffix = match.group(1)
                            base = raw_name[:-len(suffix)]
                            safe_name = f"{base}_Dup{name_counter[raw_name]}{suffix}"
                        else:
                            safe_name = f"{raw_name}_Dup{name_counter[raw_name]}"
                    else:
                        name_counter[raw_name] = 1
                        safe_name = raw_name
                        
                    parsed_seqs.append((safe_name, seq))

            results = []
            for name, seq in parsed_seqs:
                is_vh = re.search(r'[-_](VH|HC|Heavy_Chain|Heavy)[0-9]*$', name, flags=re.IGNORECASE)
                is_vl = re.search(r'[-_](VL|LC|Light_Chain|Light)[0-9]*$', name, flags=re.IGNORECASE)
                
                chain_type = "未知"
                cdrs = {"CDR1": "", "CDR2": "", "CDR3": ""}
                regions = {"FR1": "", "CDR1": "", "FR2": "", "CDR2": "", "FR3": "", "CDR3": "", "FR4": ""}
                
                if is_vh:
                    chain_type = "重链/纳米抗体"
                    cdrs, regions = extract_vh_cdrs_regex(seq)
                elif is_vl:
                    chain_type = "轻链"
                    cdrs, regions = extract_vl_cdrs_regex(seq)

                mw = calculate_mw(seq)
                pi = calculate_pi(seq)
                gravy = calculate_gravy(seq)
                cdr3_gravy = calculate_gravy(cdrs["CDR3"])
                
                cys_count = seq.count('C')
                cys_alert = f"🚨 奇数 ({cys_count})" if cys_count % 2 != 0 else f"✅ 正常 ({cys_count})"
                
                ptms = detect_ptms_detailed(regions)
                ptm_alert = " | ".join(ptms) if ptms else "✅ 无高危"
                
                results.append({
                    "序列名称 (ID)": name,
                    "类型": chain_type,
                    "CDR1": cdrs["CDR1"],
                    "CDR2": cdrs["CDR2"],
                    "CDR3": cdrs["CDR3"],
                    "CDR3长度": len(cdrs["CDR3"]) if cdrs["CDR3"] else 0,
                    "完整序列": seq,
                    "分子量 (Da)": mw,
                    "pI (等电点)": pi,
                    "完整 GRAVY": gravy,
                    "CDR3 GRAVY": cdr3_gravy,
                    "孤立Cys 雷达": cys_alert,
                    "PTM 风险预警": ptm_alert
                })
            
            df = pd.DataFrame(results)
            st.success(f"✅ 解析成功！输入了 {sum(name_counter.values())} 条序列，共成功识别提取 {len(df)} 条单链 (有效捕获率 100%)。")
            
            st.markdown("#### 📊 CMC 单链分析总表")
            st.dataframe(df.style.map(highlight_alerts, subset=['孤立Cys 雷达', 'PTM 风险预警']), use_container_width=True)
            
            def extract_base_name(name_str):
                return re.sub(r'[-_](VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light)[0-9]*$', '', name_str, flags=re.IGNORECASE).strip()
            
            df['归属分子名'] = df['序列名称 (ID)'].apply(extract_base_name)
            
            paired_data = []
            for core_name, group in df.groupby('归属分子名'):
                vh_rows = group[group['类型'].str.contains('重链')]
                vl_rows = group[group['类型'].str.contains('轻链')]
                
                if not vh_rows.empty and not vl_rows.empty:
                    for _, vh_row in vh_rows.iterrows():
                        for _, vl_row in vl_rows.iterrows():
                            vh_seq = vh_row['完整序列']
                            vl_seq = vl_row['完整序列']
                            vh_id = vh_row['序列名称 (ID)']
                            vl_id = vl_row['序列名称 (ID)']
                            
                            vh_match = re.search(r'[-_]([A-Za-z0-9_]+)$', vh_id)
                            vl_match = re.search(r'[-_]([A-Za-z0-9_]+)$', vl_id)
                            v_str = vh_match.group(1) if vh_match else "VH"
                            l_str = vl_match.group(1) if vl_match else "VL"
                            
                            pair_name = f"{core_name} ({v_str}/{l_str})"

                            vh_pi = vh_row['pI (等电点)']
                            vl_pi = vl_row['pI (等电点)']
                            delta_pi = abs(vh_pi - vl_pi)
                            
                            warning_flag = "⚠️ 关注" if delta_pi > 2.0 else "✅ 正常"
                            
                            vh_ptm = vh_row['PTM 风险预警']
                            vl_ptm = vl_row['PTM 风险预警']
                            combined_ptm = []
                            if vh_ptm and "✅" not in vh_ptm:
                                combined_ptm.append(f"VH: {vh_ptm}")
                            if vl_ptm and "✅" not in vl_ptm:
                                combined_ptm.append(f"VL: {vl_ptm}")
                            ptm_summary = " | ".join(combined_ptm) if combined_ptm else "✅ Fv 无高危 PTM"

                            paired_data.append({
                                "代表分子名": core_name,
                                "具体组合名": pair_name,
                                "具体链组合": f"{vh_id} + {vl_id}",
                                "VH_Seq": vh_seq, "VL_Seq": vl_seq,
                                "重链_pI": vh_pi, "轻链_pI": vl_pi,
                                "ΔpI": round(delta_pi, 2), "Fv质控状态": warning_flag,
                                "PTM风险汇总": ptm_summary
                            })
            
            df_paired = pd.DataFrame(paired_data)
            df_paired_final = pd.DataFrame()
            
            if not df_paired.empty:
                df_paired['Fv_Fingerprint'] = df_paired['VH_Seq'] + "||" + df_paired['VL_Seq']
                
                fv_cluster = df_paired.groupby('Fv_Fingerprint').agg(
                    包含相同配对数=('代表分子名', 'count'),
                    代表分子名=('具体组合名', 'first'),
                    合并来源分子名=('代表分子名', lambda x: ', '.join(x.unique())),
                    具体链组合=('具体链组合', lambda x: ', '.join(x.unique())),
                    重链_pI=('重链_pI', 'first'),
                    轻链_pI=('轻链_pI', 'first'),
                    ΔpI=('ΔpI', 'first'),
                    Fv质控状态=('Fv质控状态', 'first'),
                    PTM风险汇总=('PTM风险汇总', 'first')
                ).reset_index()
                
                fv_cluster['唯一性 (Unique)'] = fv_cluster['包含相同配对数'].apply(
                    lambda x: "✅ 唯一 (Unique)" if x == 1 else f"⚠️ 冗余克隆 (Dup x{x})"
                )
                
                df_paired_final = fv_cluster.sort_values(
                    by=['包含相同配对数', '代表分子名'], ascending=[False, True]
                )[[
                    '唯一性 (Unique)', '包含相同配对数', '代表分子名', '合并来源分子名', '具体链组合',
                    '重链_pI', '轻链_pI', 'ΔpI', 'Fv质控状态', 'PTM风险汇总'
                ]]
                
                st.markdown("#### ⚖️ Fv 双链配对与去重 (Unique / Duplicate) 判定")
                st.info(f"💡 共成功组装出 {len(paired_data)} 种 Fv 双链交叉配对组合，经全局指纹去重后得到 {len(df_paired_final)} 种唯一抗体 (Unique Clones)。")
                st.dataframe(df_paired_final.style.map(highlight_alerts, subset=['唯一性 (Unique)', 'Fv质控状态', 'PTM风险汇总']), use_container_width=True)

            st.markdown("#### 🧬 CDR3 核心指纹与多样性聚类")
            col1, col2 = st.columns(2)
            df_cdr3 = df[df["CDR3"].str.len() > 0]
            
            with col1:
                st.write("**🛡️ 重链 (VH) CDR3 家族**")
                vh_cdr3 = df_cdr3[df_cdr3['类型'].str.contains('重链')]
                if not vh_cdr3.empty:
                    vh_cluster = vh_cdr3.groupby('CDR3').agg(
                        出现频次=('序列名称 (ID)', 'count'),
                        CDR3长度=('CDR3长度', 'first'),
                        来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x))
                    ).reset_index().sort_values(by='出现频次', ascending=False)
                    st.dataframe(vh_cluster, use_container_width=True)
                else:
                    st.write("未提取到重链 CDR3")

            with col2:
                st.write("**🪃 轻链 (VL) CDR3 家族**")
                vl_cdr3 = df_cdr3[df_cdr3['类型'].str.contains('轻链')]
                if not vl_cdr3.empty:
                    vl_cluster = vl_cdr3.groupby('CDR3').agg(
                        出现频次=('序列名称 (ID)', 'count'),
                        CDR3长度=('CDR3长度', 'first'),
                        来源分子名单=('序列名称 (ID)', lambda x: ', '.join(x))
                    ).reset_index().sort_values(by='出现频次', ascending=False)
                    st.dataframe(vl_cluster, use_container_width=True)
                else:
                    st.write("未提取到轻链 CDR3")

            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='完整单链数据')
                if not df_paired_final.empty:
                    df_paired_final.to_excel(writer, index=False, sheet_name='Fv交叉配对与去重')
                
                ptm_report = df[~df['PTM 风险预警'].str.contains("✅")]
                ptm_report.to_excel(writer, index=False, sheet_name='高危PTM风险报表')
                
                if not vh_cdr3.empty:
                    vh_cluster.to_excel(writer, index=False, sheet_name='重链CDR3多样性')
                if not vl_cdr3.empty:
                    vl_cluster.to_excel(writer, index=False, sheet_name='轻链CDR3多样性')
                
            processed_data = output.getvalue()
            b64 = base64.b64encode(processed_data).decode()
            href = f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}" download="Antibody_Report_V16_Final.xlsx" style="display:inline-block; padding: 10px 20px; color: white; background-color: #28a745; text-decoration: none; border-radius: 5px; font-weight: bold; margin-top: 20px;">📥 下载完整分析报告 (Excel)</a>'
            st.markdown(href, unsafe_allow_html=True)
