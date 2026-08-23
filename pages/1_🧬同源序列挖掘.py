import streamlit as st
import pandas as pd
import difflib
import io
import re
import gspread

# 页面配置
st.set_page_config(page_title="抗体同源挖掘 | 聚焦 CDR", page_icon="🧬", layout="wide")

# ================= 数据库连接模块 =================
@st.cache_resource(ttl=3600)
def init_gspread():
    """纯净版初始化连接"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client
    except Exception:
        return None

@st.cache_data(ttl=600)
def load_library_from_db():
    """从云端知识库直接提取所有重轻链作为候选池"""
    client = init_gspread()
    if not client: return None
    try:
        sheet = client.open("Antibody_Database").sheet1
        data = sheet.get_all_records()
        if not data: return {}
        
        lib_dict = {}
        for row in data:
            clone_id = str(row.get('克隆ID', '')).strip()
            vh = str(row.get('重链序列', '')).strip()
            vl = str(row.get('轻链序列', '')).strip()
            
            # 分离重轻链，以符合搜索引擎的字典格式
            if clone_id:
                if vh and vh.upper() != 'NONE' and vh.upper() != 'N/A' and vh != '':
                    lib_dict[f"{clone_id}_VH"] = vh.upper()
                if vl and vl.upper() != 'NONE' and vl.upper() != 'N/A' and vl != '':
                    lib_dict[f"{clone_id}_VL"] = vl.upper()
        return lib_dict
    except Exception:
        return None

# ================= 核心算法模块 =================
@st.cache_data
def parse_fasta(fasta_text):
    sequences = {}
    current_name = ""
    current_seq = []
    
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if current_name:
                sequences[current_name] = "".join(current_seq).replace(" ", "").upper()
            current_name = line[1:].strip()
            current_seq = []
        else:
            clean_line = re.sub(r'[^A-Za-z]', '', line)
            current_seq.append(clean_line)
            
    if current_name:
        sequences[current_name] = "".join(current_seq).replace(" ", "").upper()
    return sequences

def levenshtein_distance(s1, s2):
    """计算精准的氨基酸编辑距离（支持突变、插入、缺失）"""
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for index2, char2 in enumerate(s2):
        new_distances = [index2 + 1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                new_distances.append(distances[index1])
            else:
                new_distances.append(1 + min((distances[index1], distances[index1+1], new_distances[-1])))
        distances = new_distances
    return distances[-1]

def scan_motif_in_sequence(motif, full_seq):
    """滑动窗口局部雷达扫描"""
    m_len = len(motif)
    s_len = len(full_seq)
    
    if motif in full_seq:
        return 100.0, 0, motif
        
    best_dist = m_len
    best_match = ""
    
    for w in range(max(1, m_len - 2), min(s_len + 1, m_len + 3)):
        for i in range(s_len - w + 1):
            window = full_seq[i:i+w]
            dist = levenshtein_distance(motif, window)
            if dist < best_dist:
                best_dist = dist
                best_match = window
                
    sim = max(0.0, ((m_len - best_dist) / m_len) * 100)
    return sim, best_dist, best_match

def calculate_global_similarity(seq1, seq2):
    """传统的全局相似度（包含骨架）"""
    matcher = difflib.SequenceMatcher(None, seq1, seq2)
    identity = matcher.ratio()
    match_len = sum(triple.size for triple in matcher.get_matching_blocks())
    mutations = max(len(seq1), len(seq2)) - match_len
    return identity * 100, mutations

# ================= 侧边栏 UI =================
with st.sidebar:
    st.header("⚙️ 挖掘引擎设置")
    
    mode = st.radio(
        "选择比对策略：",
        ["🎯 CDR/局部基序扫描 (推荐)", "🌍 全局全长比对"]
    )
    
    st.markdown("---")
    if "CDR" in mode:
        st.info("当前模式：**无视骨架差异**，专门在库中搜寻与种子 CDR 高度相似的片段。")
        sim_threshold = st.slider("最低 CDR 相似度阈值 (%)", 50.0, 99.9, 80.0, 1.0)
    else:
        st.info("当前模式：对比整条序列（包含 FR 区），适合寻找完全属于同一个大克隆家族的兄弟分子。")
        sim_threshold = st.slider("最低全局相似度阈值 (%)", 50.0, 99.9, 85.0, 1.0)

    top_k = st.number_input("每个种子最多保留命中数", min_value=1, max_value=500, value=50)

# ================= 主界面 =================
st.title("🧬 抗体同源挖掘引擎 | 深度定制作战版")

col1, col2 = st.columns(2)
with col1:
    if "CDR" in mode:
        st.subheader("🥇 第一步：输入种子 CDR (Motif)")
        seed_input = st.text_area("请提取优选分子的 CDR 序列：", height=200, 
                                  placeholder=">Seed_1_CDRH3\nARQGYGMDVW\n>Seed_2_CDRH123\nGFNIKDTY-RIDPANGN-YYGMDY")
    else:
        st.subheader("🥇 第一步：输入优选分子 (全长)")
        seed_input = st.text_area("请在此粘贴优选分子的完整 FASTA 序列：", height=200, 
                                  placeholder=">Seed_Clone_1\nEVQLVQSGAEVKKPGASVKVSCKASGYTFT...")

with col2:
    st.subheader("📚 第二步：选择候选序列库")
    
    # 新增：让用户选择是从云端加载，还是手动输入
    lib_source = st.radio("选择底层大库来源：", 
                          ["☁️ 自动调取云端知识库 (Antibody_Database)", "✍️ 手动粘贴 FASTA"], horizontal=True)
    
    if lib_source == "✍️ 手动粘贴 FASTA":
        lib_input = st.text_area("请在此粘贴海量候选库的 FASTA 全序列：", height=155)
    else:
        st.info("💡 系统将自动连接云端，静默下载您积累的所有重链(VH)和轻链(VL)作为本次挖掘的庞大候选池。")
        lib_input = ""

# ================= 数据处理与挖掘 =================
if st.button("🚀 启动扫描引擎", type="primary", use_container_width=True):
    if not seed_input.strip():
        st.error("⚠️ 请确保优选种子已输入数据！")
    elif lib_source == "✍️ 手动粘贴 FASTA" and not lib_input.strip():
        st.error("⚠️ 请粘贴候选序列库！")
    else:
        with st.spinner("雷达正在穿透骨架扫描靶心序列，请稍候..."):
            seeds = parse_fasta(seed_input)
            
            # 根据用户的选择，决定从哪里读取 library
            if lib_source == "✍️ 手动粘贴 FASTA":
                library = parse_fasta(lib_input)
            else:
                library = load_library_from_db()
                if library is None:
                    st.error("⚠️ 无法连接云端数据库，请检查密钥或回退到手动粘贴模式。")
                    st.stop()
                elif not library:
                    st.warning("⚠️ 云端数据库目前为空，请先前往主程序入库序列。")
                    st.stop()
                else:
                    st.success(f"✅ 成功从云端知识库装载 **{len(library)}** 条重/轻链序列参战！")
            
            if not seeds or not library:
                st.error("⚠️ 序列解析失败，请检查输入格式。")
            else:
                results = []
                progress_bar = st.progress(0)
                
                for i, (seed_name, seed_seq) in enumerate(seeds.items()):
                    seed_matches = []
                    seed_seq_clean = seed_seq.replace("-", "") 
                    
                    for lib_name, lib_seq in library.items():
                        if seed_seq_clean == lib_seq: continue
                        
                        if "CDR" in mode:
                            identity, mutations, best_match = scan_motif_in_sequence(seed_seq_clean, lib_seq)
                        else:
                            identity, mutations = calculate_global_similarity(seed_seq_clean, lib_seq)
                            best_match = lib_seq
                        
                        if identity >= sim_threshold:
                            seed_matches.append({
                                '种子名称': seed_name,
                                '命中序列名称': lib_name,
                                '核心相似度 (%)': round(identity, 2),
                                '氨基酸突变数': mutations,
                                '局部命中片段': best_match if "CDR" in mode else "全局匹配"
                            })
                    
                    seed_matches = sorted(seed_matches, key=lambda x: x['核心相似度 (%)'], reverse=True)[:top_k]
                    results.extend(seed_matches)
                    progress_bar.progress((i + 1) / len(seeds))
                
                # ================= 结果展示 =================
                if results:
                    df_results = pd.DataFrame(results)
                    st.success(f"挖掘完成！共为您锁定 **{len(df_results)}** 条具有高同源价值的变体序列！")
                    
                    def highlight_sim(val):
                        if isinstance(val, float):
                            if val >= 95.0: return 'background-color: #d4edda; color: #155724; font-weight: bold;'
                            elif val >= 85.0: return 'background-color: #fff3cd; color: #856404;'
                        return ''
                    
                    st.dataframe(df_results.style.map(highlight_sim, subset=['核心相似度 (%)']), use_container_width=True)
                    
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='openpyxl') as writer:
                        df_results.to_excel(writer, index=False, sheet_name='同源挖掘结果')
                    output.seek(0)
                    
                    st.download_button(
                        label="📥 导出同源挖掘报告 (Excel)",
                        data=output,
                        file_name="CDR_同源分子挖掘报告.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        type="primary"
                    )
                else:
                    st.warning("🧐 在当前阈值下，未能找到相似变体。")
