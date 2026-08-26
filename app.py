import streamlit as st
import pandas as pd
import re
import io
import difflib

# ==================== 新增 1：连接 Google Sheets 的入库引擎 ====================
import gspread
from datetime import datetime

@st.cache_resource(ttl=3600)
def init_gspread():
    """纯净版初始化连接 (用于写入)"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None, "未在后台 Secrets 检测到配置"
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client, "OK"
    except Exception as e:
        return None, str(e)

def save_batch_to_database(target, notes, full_df):
    """批量将全量无冗余抗体记录追加到云端表格"""
    client, error_msg = init_gspread()
    if not client:
        return False, f"🤖 机器人启动失败，原因: {error_msg}"
    
    try:
        sheet = client.open("Antibody_Database").sheet1
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M")
        
        # 准备批量写入的数据结构
        rows_to_insert = []
        for _, row in full_df.iterrows():
            clone_id = row['代表分子名']
            vh_seq = row['_VH_Seq']
            vl_seq = row['_VL_Seq']
            gravy = "N/A" # app.py 暂未计算整体 GRAVY，先用 N/A 占位
            hic_score = "N/A" # 同上
            qc_status = row['Fv质控状态']
            
            # 组装数据，顺序必须严格对应表头：
            # [克隆ID, 靶点, 入库时间, 重链序列, 轻链序列, GRAVY, Max_HIC_Score, 质控状态, 备注]
            new_row = [
                clone_id, 
                target, 
                current_time, 
                vh_seq, 
                vl_seq, 
                gravy, 
                hic_score, 
                qc_status, 
                notes
            ]
            rows_to_insert.append(new_row)
        
        if rows_to_insert:
            # ====== 修复覆盖问题的核心防线 ======
            existing_rows = sheet.get_all_values()
            next_empty_row = len(existing_rows) + 1
            
            # 使用 update 进行绝对坐标的精准打击
            sheet.update(f"A{next_empty_row}", rows_to_insert)
            # ====================================
            return True, f"成功入库 {len(rows_to_insert)} 条序列！"
        else:
            return False, "没有可以入库的数据。"
    
    except gspread.exceptions.SpreadsheetNotFound:
        return False, "机器人找不到名叫 Antibody_Database 的表格，请检查名字和共享权限。"
    except Exception as e:
        return False, str(e)
# ==============================================================================

# 设置全局页面配置
st.set_page_config(page_title="工业级抗体生信大屏 V16", page_icon="🧬", layout="wide")

# 初始化 Session State
if 'analysis_started' not in st.session_state:
    st.session_state.analysis_started = False

# 采用 ExPASy 主流 pKa 参考值
PI_DICT = {
    'D': -3.90, 'E': -4.07, 'C': -8.18, 'Y': -10.46,
    'H': 6.04, 'K': 10.54, 'R': 12.48,
    'N_term': 8.0, 'C_term': -3.1
}

def calculate_pi(sequence):
    """基于 Henderson-Hasselbalch 方程和二分法迭代计算精确的 pI 值"""
    if not sequence or not isinstance(sequence, str):
        return None
    
    sequence = sequence.upper()
    counts = {aa: sequence.count(aa) for aa in PI_DICT.keys() if len(aa) == 1}
    
    def net_charge(pH):
        charge = 0.0
        charge += counts.get('H', 0) / (1.0 + 10**(pH - PI_DICT['H']))
        charge += counts.get('K', 0) / (1.0 + 10**(pH - PI_DICT['K']))
        charge += counts.get('R', 0) / (1.0 + 10**(pH - PI_DICT['R']))
        charge += 1.0 / (1.0 + 10**(pH - PI_DICT['N_term']))
        
        charge -= counts.get('D', 0) / (1.0 + 10**(PI_DICT['D'] * -1 - pH))
        charge -= counts.get('E', 0) / (1.0 + 10**(PI_DICT['E'] * -1 - pH))
        charge -= counts.get('C', 0) / (1.0 + 10**(PI_DICT['C'] * -1 - pH))
        charge -= counts.get('Y', 0) / (1.0 + 10**(PI_DICT['Y'] * -1 - pH))
        charge -= 1.0 / (1.0 + 10**(PI_DICT['C_term'] * -1 - pH))
        return charge
        
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2.0
        if net_charge(mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)

def scan_ptm(seq):
    ptms = []
    for m in re.finditer(r'N[^P][ST]', seq): ptms.append(f"N-糖基化 (NIT) @{m.start()+1}")
    for m in re.finditer(r'NG', seq): ptms.append(f"脱酰胺 (NG) @{m.start()+1}")
    for m in re.finditer(r'DG', seq): ptms.append(f"异构化 (DG) @{m.start()+1}")
    for m in re.finditer(r'DP', seq): ptms.append(f"酸断裂 (DP) @{m.start()+1}")
        
    c_count = seq.count('C')
    if c_count > 2:
        c_positions = [i+1 for i, res in enumerate(seq) if res == 'C']
        ptms.append(f"🚨游离Cys风险: 发现 {c_count} 个 @ {c_positions}")
    elif c_count == 1:
        c_positions = [i+1 for i, res in enumerate(seq) if res == 'C']
        ptms.append(f"🚨结构缺陷: 仅 1 个半胱氨酸 @ {c_positions}")
        
    return " | ".join(ptms) if ptms else "无高危 PTM"

@st.cache_data
def parse_fasta(fasta_text):
    sequences = []
    current_id = ""
    current_seq = []
    for line in fasta_text.strip().split('\n'):
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if current_id:
                sequences.append({'ID': current_id, 'Sequence': "".join(current_seq)})
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        sequences.append({'ID': current_id, 'Sequence': "".join(current_seq)})
    return pd.DataFrame(sequences)

st.title("🔬 工业级抗体序列生信大屏中台 (V16)")
st.markdown("集成了 **pI 计算**、**CDR3 富集聚类**、**双链组合**、**发散性挑选** 以及 **全量云端知识库直连**。")

fasta_input = st.text_area("请粘贴全量 FASTA 格式数据 (支持批量混合粘贴重/轻链):", height=200,
                          help="请以 >分子名_链类型 (如 >F0630-4D1-VH, >F0630-6B4-VL2) 的格式命名")

col1, col2 = st.columns([1, 5])
with col1:
    if st.button("🚀 开始极速分析", use_container_width=True):
        st.session_state.analysis_started = True

if st.session_state.analysis_started and fasta_input:
    with st.spinner('🚀 正在进行工业级序列解析与组装...'):
        df = parse_fasta(fasta_input)
        
        if df.empty:
            st.error("未检测到有效数据，请检查格式。")
        else:
            processed_data = []
            for _, row in df.iterrows():
                seq_id = row['ID']
                seq = row['Sequence']
                
                # 👇 升级了这里的正则表达式，增加 ([_ \-]?\d*) 允许 VH 后面带有横杠或下划线的数字后缀 (如 VH-1, VL_2)
                match = re.search(r"^(.*)[_ \-](VH|VL|VHH|LC|HC)([_ \-]?\d*)$", seq_id, re.IGNORECASE)
                if match:
                    base_name = match.group(1).strip()
                    chain_type = match.group(2).upper()
                    suffix = match.group(3).strip()
                    if chain_type in ['VH', 'VHH', 'HC']: c_type = '重链 (VH)'
                    elif chain_type in ['VL', 'LC']: c_type = '轻链 (VL)'
                    else: c_type = '未知'
                    chain_label = f"{chain_type}{suffix}"
                else:
                    base_name = seq_id
                    c_type = '未知'
                    chain_label = '未知'
                
                cdr3_seq = "N/A"
                cys_matches = [m.start() for m in re.finditer(r'C', seq)]
                trp_matches = [m.start() for m in re.finditer(r'W', seq)]
                
                if len(cys_matches) >= 2 and trp_matches:
                    cys_2 = cys_matches[-1] if cys_matches[-1] < 110 else (cys_matches[-2] if len(cys_matches)>1 else -1)
                    valid_trps = [w for w in trp_matches if w > cys_2]
                    if cys_2 != -1 and valid_trps:
                        trp = valid_trps[0]
                        if 3 < (trp - cys_2) < 30:
                            cdr3_seq = seq[cys_2 + 1 : trp]
                
                processed_data.append({
                    '序列ID': seq_id, '分子基名': base_name, '链标签': chain_label,
                    '链类型': c_type, '预测CDR3': cdr3_seq, '长度': len(seq),
                    '理论pI': calculate_pi(seq), '高危 PTM 风险预警': scan_ptm(seq),
                    '全长序列': seq
                })
                
            df_all = pd.DataFrame(processed_data)
            
            st.subheader("🧬 单链 CDR3 克隆聚类与丰度分析")
            cdr3_cluster = df_all[df_all['预测CDR3'] != "N/A"].groupby(['链类型', '预测CDR3']).agg(
                克隆数=('序列ID', 'count'), 包含序列ID=('序列ID', lambda x: ", ".join(x)),
                分子基名集合=('分子基名', lambda x: ", ".join(set(x)))
            ).reset_index()
            
            cdr3_cluster = cdr3_cluster.sort_values(by=['链类型', '克隆数'], ascending=[False, False])
            cdr3_cluster.insert(0, '丰度状态', cdr3_cluster['克隆数'].apply(lambda x: '🔥 优势富集 (x{})'.format(x) if x > 1 else '唯一'))
            df_cdr3_final = cdr3_cluster
            
            def highlight_cdr3(row):
                if '🔥' in str(row['丰度状态']): return ['background-color: #ffe8cc; color: #d97706; font-weight: bold'] * len(row)
                return [''] * len(row)
            st.dataframe(df_cdr3_final.style.apply(highlight_cdr3, axis=1), use_container_width=True)
            
            df_vh = df_all[df_all['链类型'] == '重链 (VH)']
            df_vl = df_all[df_all['链类型'] == '轻链 (VL)']
            
            paired_data = []
            if not df_vh.empty and not df_vl.empty:
                common_bases = set(df_vh['分子基名']).intersection(set(df_vl['分子基名']))
                for base in common_bases:
                    vh_candidates = df_vh[df_vh['分子基名'] == base].to_dict('records')
                    vl_candidates = df_vl[df_vl['分子基名'] == base].to_dict('records')
                    for vh in vh_candidates:
                        for vl in vl_candidates:
                            vh_id = vh['序列ID']
                            vl_id = vl['序列ID']
                            combo_name = f"{base} ({vh['链标签']}/{vl['链标签']})" if len(vh_candidates)>1 or len(vl_candidates)>1 else base
                            vh_pi = vh['理论pI']
                            vl_pi = vl['理论pI']
                            delta_pi = round(abs(vh_pi - vl_pi), 2) if vh_pi and vl_pi else None
                            fv_fingerprint = f"{vh['全长序列']}||{vl['全长序列']}"
                            
                            vh_ptm, vl_ptm = vh['高危 PTM 风险预警'], vl['高危 PTM 风险预警']
                            ptm_summary = ""
                            if vh_ptm != "无高危 PTM": ptm_summary += f"VH: {vh_ptm} | "
                            if vl_ptm != "无高危 PTM": ptm_summary += f"VL: {vl_ptm}"
                            ptm_summary = ptm_summary.strip(" | ") if ptm_summary else "☑️ Fv 无高危 PTM"
                            
                            paired_data.append({
                                '组合分子名': combo_name, '具体链组合': f"重链: {vh_id} | 轻链: {vl_id}",
                                '指纹 (Fingerprint)': fv_fingerprint, '重链_pI': vh_pi, '轻链_pI': vl_pi,
                                'ΔpI': delta_pi, 'PTM风险汇总': ptm_summary,
                                '_VH_ID': vh_id, '_VL_ID': vl_id, '_VH_Seq': vh['全长序列'], '_VL_Seq': vl['全长序列']
                            })
                            
            df_paired_full = None 
            df_paired_final = None
            df_diverse_final = None
            
            if paired_data:
                st.subheader("🔗 Fv 完整双链配对与冗余聚类结果")
                df_paired_raw = pd.DataFrame(paired_data)
                fv_cluster = df_paired_raw.groupby('指纹 (Fingerprint)').agg(
                    count=('组合分子名', 'count'), rep_name=('组合分子名', 'first'),
                    merged_names=('组合分子名', lambda x: ", ".join(x)), chain_details=('具体链组合', lambda x: " /// ".join(set(x)))
                ).reset_index()
                
                cluster_data = []
                for _, row in fv_cluster.iterrows():
                    fp, count, rep_name, merged_names, chain_details = row['指纹 (Fingerprint)'], row['count'], row['rep_name'], row['merged_names'], row['chain_details']
                    group = df_paired_raw[df_paired_raw['指纹 (Fingerprint)'] == fp]
                    unique_flag = "✅ 唯一 (Unique)" if count == 1 else "⚠️ 冗余 (Redundant)"
                    delta_pi = group.iloc[0]['ΔpI']
                    qc_status = "⚠️ 关注: ΔpI过大" if delta_pi and delta_pi > 2.0 else "✅ 正常"
                        
                    cluster_data.append({
                        '唯一性 (Unique)': unique_flag, '包含相同配对数': count, '代表分子名': rep_name,
                        '具体链组合': chain_details, '合并来源分子名': merged_names, '重链_pI': group.iloc[0]['重链_pI'],
                        '轻链_pI': group.iloc[0]['轻链_pI'], 'ΔpI': delta_pi, 'Fv质控状态': qc_status,
                        'PTM风险汇总': group.iloc[0]['PTM风险汇总'], '_Fv_Seq_Fingerprint': fp,
                        '_VH_ID': group.iloc[0]['_VH_ID'], '_VL_ID': group.iloc[0]['_VL_ID'],
                        '_VH_Seq': group.iloc[0]['_VH_Seq'], '_VL_Seq': group.iloc[0]['_VL_Seq']
                    })
                    
                df_paired_full = pd.DataFrame(cluster_data).sort_values(by=['包含相同配对数', '代表分子名'], ascending=[False, True])
                display_cols = ['唯一性 (Unique)', '包含相同配对数', '代表分子名', '具体链组合', '合并来源分子名', '重链_pI', '轻链_pI', 'ΔpI', 'Fv质控状态', 'PTM风险汇总']
                df_paired_final = df_paired_full[display_cols]
                
                def highlight_fv(row):
                    colors = [''] * len(row)
                    if '⚠️ 冗余' in str(row.get('唯一性 (Unique)', '')): colors[row.index.get_loc('唯一性 (Unique)')] = 'background-color: #ffebb5; color: black;'
                    if '⚠️ 关注' in str(row.get('Fv质控状态', '')): colors[row.index.get_loc('Fv质控状态')] = 'background-color: #e2e3e5; color: #383d41;'
                    if 'VH' in str(row.get('PTM风险汇总', '')) or 'VL' in str(row.get('PTM风险汇总', '')): colors[row.index.get_loc('PTM风险汇总')] = 'background-color: #f8d7da; color: #721c24;'
                    return colors
                st.dataframe(df_paired_final.style.apply(highlight_fv, axis=1), use_container_width=True)
                
                st.subheader("🌈 智能序列差异化推荐池 (Diversity Picker)")
                num_unique_clones = len(df_paired_full)
                if num_unique_clones > 2:
                    target_n = st.slider("请选择计划挑取进行下游验证的克隆数量：", min_value=2, max_value=min(200, num_unique_clones), value=min(5, num_unique_clones))
                    pool = df_paired_full.to_dict('records')
                    def calc_dist(seq1, seq2): return 1.0 - difflib.SequenceMatcher(None, seq1, seq2).ratio()
                    
                    dist_matrix = {}
                    for i in range(num_unique_clones):
                        for j in range(i+1, num_unique_clones):
                            d = calc_dist(pool[i]['_Fv_Seq_Fingerprint'], pool[j]['_Fv_Seq_Fingerprint'])
                            dist_matrix[(i, j)] = d
                            dist_matrix[(j, i)] = d
                            
                    max_d = -1; seed1, seed2 = 0, 1
                    for i in range(num_unique_clones):
                        for j in range(i+1, num_unique_clones):
                            if dist_matrix[(i, j)] > max_d: max_d = dist_matrix[(i, j)]; seed1, seed2 = i, j
                                
                    selected_indices = [seed1, seed2]
                    unselected_indices = [i for i in range(num_unique_clones) if i not in selected_indices]
                    
                    while len(selected_indices) < target_n and unselected_indices:
                        best_cand, max_of_mins = -1, -1
                        for cand in unselected_indices:
                            min_dist_to_sel = min([dist_matrix[(cand, sel)] for sel in selected_indices])
                            if min_dist_to_sel > max_of_mins: max_of_mins = min_dist_to_sel; best_cand = cand
                        selected_indices.append(best_cand)
                        unselected_indices.remove(best_cand)
                    
                    df_diverse = df_paired_full.iloc[selected_indices].copy()
                    df_diverse.insert(0, '建议纯化优先级', [f"🥇 Top {i+1}" for i in range(len(df_diverse))])
                    df_diverse_final = df_diverse
                    
                    st.success(f"已为您锁定 {target_n} 个在全局序列空间中发散度最高的种子克隆：")
                    st.dataframe(df_diverse_final[['建议纯化优先级'] + display_cols].style.apply(highlight_fv, axis=1), use_container_width=True)
                else:
                    st.info("独立克隆数量不足以触发多样性挑选。")

                # ==================== 新增 2：全量保存至云端知识库的 UI 面板 ====================
                st.markdown("---")
                st.markdown("### 💾 全量资产沉淀：写入云端知识库")
                st.info(f"💡 无论是否被推荐，点击下方按钮，将本次去重后的 **全部 {len(df_paired_full)} 条** 唯一双链 Fv 序列存入 Google Sheets。")
                
                with st.form(key="batch_save_form"):
                    col_tgt, col_nts = st.columns(2)
                    with col_tgt:
                        batch_target = st.text_input("统一设置靶点名称 (Target):", placeholder="例如: PD-1")
                    with col_nts:
                        batch_notes = st.text_input("统一添加备注 (Notes):", placeholder="全量解析入库")
                    
                    submit_btn = st.form_submit_button("🚀 一键全量入库 (Save ALL Unique Fvs to DB)")
                    
                    if submit_btn:
                        if not batch_target:
                            st.warning("⚠️ 建议输入一个靶点名称，方便在知识库中归类检索。")
                        with st.spinner(f"正在将 {len(df_paired_full)} 条全量无冗余序列推送到 Google Drive..."):
                            # 注意这里传入的是 df_paired_full (全量去重库)，而不是推荐池
                            success, msg = save_batch_to_database(
                                target=batch_target, 
                                notes=batch_notes, 
                                full_df=df_paired_full
                            )
                            if success:
                                st.success(f"🎉 {msg} 您可以前往左侧『3_📂序列知识库』查看最新大盘数据。")
                            else:
                                st.error(f"⚠️ 全量入库失败: {msg}")
                # ==============================================================================

            else:
                st.info("未成功配对任何有效 Fv 结构。")

            if len(df_all) > 0:
                output = io.BytesIO()
                with pd.ExcelWriter(output, engine='openpyxl') as writer:
                    df_all.to_excel(writer, index=False, sheet_name='完整单链数据')
                    if df_cdr3_final is not None and not df_cdr3_final.empty:
                        df_cdr3_final.to_excel(writer, index=False, sheet_name='CDR3克隆聚类')
                    if df_paired_final is not None and not df_paired_final.empty:
                        df_paired_final.to_excel(writer, index=False, sheet_name='Fv组装与排重')
                    if df_diverse_final is not None and not df_diverse_final.empty:
                        df_export_diverse = df_diverse_final[['建议纯化优先级'] + display_cols].copy()
                        df_export_diverse['重链序列'] = df_diverse_final['_VH_Seq']
                        df_export_diverse['轻链序列'] = df_diverse_final['_VL_Seq']
                        df_export_diverse.to_excel(writer, index=False, sheet_name='推荐纯化池')
                
                st.markdown("---")
                st.download_button(
                    label="💾 下载本地 Excel 多维度生信分析报告",
                    data=output.getvalue(),
                    file_name="工业级抗体大屏分析报告_V16_Final.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
