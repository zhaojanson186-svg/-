import streamlit as st
import pandas as pd
import re
import io
from collections import defaultdict

PI_DICT = {'D': 2.77, 'E': 3.22, 'C': 8.18, 'Y': 10.46, 'H': 6.00, 'K': 10.53, 'R': 12.48}
GRAVY_DICT = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def calculate_pi(sequence):
    if not sequence: return None
    charge = 0.0
    for aa in sequence:
        if aa in PI_DICT:
            if aa in ['K', 'R', 'H']: charge += 1
            elif aa in ['D', 'E', 'C', 'Y']: charge -= 1
    return round(7.0 + (charge * 0.1), 2) 

def calculate_gravy(sequence):
    if not sequence: return None
    total = sum(GRAVY_DICT.get(aa, 0) for aa in sequence)
    return round(total / len(sequence), 2)

def extract_base_name(name):
    # 标准 flags=re.IGNORECASE，且 \d*$ 兼容 VH1, VL2 等带数字的后缀
    base = re.sub(r'[-_](VH|VL|HC|LC|Heavy_Chain|Light_Chain|Heavy|Light)\d*$', '', name, flags=re.IGNORECASE).strip()
    return base

def detect_ptms_detailed(sequence):
    if not sequence: return "未知"
    alerts = []
    
    # 极简高危 PTM 规则：只抓最致命的缺陷
    ptm_motifs = {
        'N-糖基化 (NXT/NXS)': r'N[^P][ST]',
        '极速脱氨基 (NG)': r'NG',
        '极速异构化 (DG)': r'DG',
        '酸断裂点 (DP)': r'DP'
    }
    
    regions = {
        'FR1': sequence[:25], 'CDR1': sequence[25:35], 'FR2': sequence[35:50],
        'CDR2': sequence[50:65], 'FR3': sequence[65:95], 'CDR3': sequence[95:115],
        'FR4': sequence[115:]
    }
    
    # 只对 CDR 区域的高危 PTM 报警
    for region_name, region_seq in regions.items():
        if region_name.startswith("CDR"):  
            for ptm_name, pattern in ptm_motifs.items():
                for match in re.finditer(pattern, region_seq):
                    global_pos = len(''.join(list(regions.values())[:list(regions.keys()).index(region_name)])) + match.start() + 1
                    alerts.append(f"[{region_name}] {ptm_name} @{global_pos}")
                    
    return " | ".join(alerts) if alerts else "无高危 PTM"

def parse_fasta(fasta_text):
    sequences = []
    current_id = ""
    current_seq = ""
    id_counts = defaultdict(int)
    
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if current_id and current_seq:
                id_counts[current_id] += 1
                final_id = current_id if id_counts[current_id] == 1 else f"{current_id}_Dup{id_counts[current_id]}"
                sequences.append({"id": final_id, "seq": current_seq})
            current_id = line[1:].strip()
            current_seq = ""
        else:
            current_seq += line.upper()
            
    if current_id and current_seq:
        id_counts[current_id] += 1
        final_id = current_id if id_counts[current_id] == 1 else f"{current_id}_Dup{id_counts[current_id]}"
        sequences.append({"id": final_id, "seq": current_seq})
        
    parsed_data = []
    
    for item in sequences:
        seq_id = item["id"]
        seq = item["seq"]
        
        # 兼容带数字后缀的轻重链判定 (VH1, VL2)
        chain_type = "未知"
        if re.search(r'[-_](VH|HC|Heavy)\d*$', seq_id, flags=re.IGNORECASE):
            chain_type = "重链 (VH)"
        elif re.search(r'[-_](VL|LC|Light)\d*$', seq_id, flags=re.IGNORECASE):
            chain_type = "轻链 (VL)"
            
        cdr3_seq = "解析失败"
        cdr3_len = 0
        match = re.search(r'C([A-Z]{3,25}?)[WF]G.G', seq)
        if match:
            cdr3_seq = match.group(1)
            cdr3_len = len(cdr3_seq)
            
        parsed_data.append({
            '序列名称 (ID)': seq_id,
            '链类型': chain_type,
            '归属分子名': extract_base_name(seq_id),
            '全长序列': seq,
            '理论等电点 (pI)': calculate_pi(seq),
            'CDR3序列 (预估)': cdr3_seq,
            'CDR3长度': cdr3_len,
            '全长疏水指数 (GRAVY)': calculate_gravy(seq),
            '高危 PTM 风险预警': detect_ptms_detailed(seq)
        })
        
    return pd.DataFrame(parsed_data)

st.set_page_config(page_title="工业级抗体生信大屏_V16", layout="wide")
st.title("🧬 工业级抗体序列质控与 Fv 配对中台 (V16 最终版)")

with st.sidebar:
    st.header("📥 数据输入区")
    fasta_input = st.text_area("请粘贴 FASTA 序列 (支持数百条序列混合):", height=300)
    process_btn = st.button("🚀 开始极速分析", type="primary")

if process_btn and fasta_input:
    with st.spinner("正在进行智能去重、序列切片与 PTM 雷达扫描..."):
        df_all = parse_fasta(fasta_input)
        
        st.subheader(f"📊 单链全景质控总表 (共识别 {len(df_all)} 条序列)")
        st.dataframe(df_all, use_container_width=True)
        
        st.subheader("🔗 Fv 双链组装与聚类排重评估")
        
        paired_data = []
        for base_name, group in df_all.groupby('归属分子名'):
            vhs = group[group['链类型'] == '重链 (VH)']
            vls = group[group['链类型'] == '轻链 (VL)']
            
            if not vhs.empty and not vls.empty:
                for _, vh in vhs.iterrows():
                    for _, vl in vls.iterrows():
                        vh_id = vh['序列名称 (ID)']
                        vl_id = vl['序列名称 (ID)']
                        
                        # 兼容多链命名 (例如：展示为 F0630-24C7 (VH1/VL2))
                        combo_name = base_name if (len(vhs)==1 and len(vls)==1) else f"{base_name} ({vh_id}/{vl_id})"
                        
                        fv_fingerprint = vh['全长序列'] + "||" + vl['全长序列']
                        
                        vh_pi = vh['理论等电点 (pI)']
                        vl_pi = vl['理论等电点 (pI)']
                        delta_pi = round(abs(vh_pi - vl_pi), 2) if vh_pi and vl_pi else None
                        
                        vh_ptm = vh['高危 PTM 风险预警']
                        vl_ptm = vl['高危 PTM 风险预警']
                        ptm_summary = ""
                        if vh_ptm != "无高危 PTM": ptm_summary += f"VH: {vh_ptm} | "
                        if vl_ptm != "无高危 PTM": ptm_summary += f"VL: {vl_ptm}"
                        ptm_summary = ptm_summary.strip(" | ") if ptm_summary else "☑️ Fv 无高危 PTM"
                        
                        paired_data.append({
                            '组合分子名': combo_name,
                            '具体链组合': f"重链: {vh_id} | 轻链: {vl_id}",
                            '指纹 (Fingerprint)': fv_fingerprint,
                            '重链_pI': vh_pi,
                            '轻链_pI': vl_pi,
                            'ΔpI': delta_pi,
                            'PTM风险汇总': ptm_summary
                        })
                        
        if paired_data:
            df_paired_raw = pd.DataFrame(paired_data)
            cluster_data = []
            
            for fp, group in df_paired_raw.groupby('指纹 (Fingerprint)'):
                count = len(group)
                rep_name = group.iloc[0]['组合分子名']
                chain_details = group.iloc[0]['具体链组合']
                merged_names = ", ".join(group['组合分子名'].tolist())
                unique_flag = "☑️ 唯一 (Unique)" if count == 1 else f"⚠️ 冗余 (Dup x{count})"
                
                delta_pi = group.iloc[0]['ΔpI']
                qc_status = "✅ 正常"
                # ΔpI 过大降级为关注提示，不报红
                if delta_pi and delta_pi > 2.0:
                    qc_status = "⚠️ 关注: ΔpI过大"
                    
                cluster_data.append({
                    '唯一性 (Unique)': unique_flag,
                    '包含相同配对数': count,
                    '代表分子名': rep_name,
                    '具体链组合': chain_details,
                    '合并来源分子名': merged_names,
                    '重链_pI': group.iloc[0]['重链_pI'],
                    '轻链_pI': group.iloc[0]['轻链_pI'],
                    'ΔpI': delta_pi,
                    'Fv质控状态': qc_status,
                    'PTM风险汇总': group.iloc[0]['PTM风险汇总']
                })
                
            # 先排序，再进行列切片，彻底避免 KeyError
            df_paired_final = pd.DataFrame(cluster_data).sort_values(
                by=['包含相同配对数', '代表分子名'], ascending=[False, True]
            )[[
                '唯一性 (Unique)', '包含相同配对数', '代表分子名', '具体链组合', '合并来源分子名', 
                '重链_pI', '轻链_pI', 'ΔpI', 'Fv质控状态', 'PTM风险汇总'
            ]]
            
            def highlight_fv(row):
                colors = [''] * len(row)
                if '⚠️ 冗余' in str(row['唯一性 (Unique)']):
                    colors[row.index.get_loc('唯一性 (Unique)')] = 'background-color: #ffebb5; color: black;'
                if '⚠️ 关注' in str(row['Fv质控状态']):
                    colors[row.index.get_loc('Fv质控状态')] = 'background-color: #e2e3e5; color: #383d41;' # 降级为灰色信息
                if 'VH' in str(row['PTM风险汇总']) or 'VL' in str(row['PTM风险汇总']):
                    colors[row.index.get_loc('PTM风险汇总')] = 'background-color: #f8d7da; color: #721c24;'
                return colors

            st.dataframe(df_paired_final.style.apply(highlight_fv, axis=1), use_container_width=True)
            
            if len(df_all) > 0:
                output = io.BytesIO()
                # 引擎替换为系统标配的 openpyxl，修复导出报错
                with pd.ExcelWriter(output, engine='openpyxl') as writer:
                    df_all.to_excel(writer, index=False, sheet_name='完整单链数据')
                    df_paired_final.to_excel(writer, index=False, sheet_name='Fv组装与排重')
                
                st.download_button(
                    label="💾 下载完整生信分析报告 (Excel)",
                    data=output.getvalue(),
                    file_name="工业级抗体大屏分析报告_V16版.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
        else:
            st.info("尚未识别到可以成功成对的 VH/VL 序列。")
