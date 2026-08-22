import streamlit as st
import pandas as pd
import difflib
import io
import re

# 页面配置
st.set_page_config(page_title="抗体同源序列智能挖掘系统", page_icon="🧬", layout="wide")

# ================= 核心算法模块 =================
@st.cache_data
def parse_fasta(fasta_text):
    """解析 FASTA 格式数据"""
    sequences = {}
    current_name = ""
    current_seq = []
    
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name:
                sequences[current_name] = "".join(current_seq).replace(" ", "").upper()
            current_name = line[1:].strip()
            current_seq = []
        else:
            # 移除非字母字符（如数字、空格）
            clean_line = re.sub(r'[^A-Za-z]', '', line)
            current_seq.append(clean_line)
            
    if current_name:
        sequences[current_name] = "".join(current_seq).replace(" ", "").upper()
        
    return sequences

def calculate_similarity(seq1, seq2):
    """计算两条序列的相似度和氨基酸差异数"""
    matcher = difflib.SequenceMatcher(None, seq1, seq2)
    identity = matcher.ratio()
    
    # 粗略计算突变数（基于比对块）
    match_len = sum(triple.size for triple in matcher.get_matching_blocks())
    mutations = max(len(seq1), len(seq2)) - match_len
    
    return identity, mutations

# ================= 侧边栏 UI =================
with st.sidebar:
    st.header("⚙️ 挖掘参数设置")
    sim_threshold = st.slider("最低同源相似度阈值 (%)", min_value=50.0, max_value=99.9, value=85.0, step=0.5, 
                              help="低于此相似度的序列将被系统自动过滤")
    top_k = st.number_input("每个种子最多保留的相似克隆数", min_value=1, max_value=500, value=50)
    
    st.markdown("---")
    st.markdown("### 💡 挖掘策略建议")
    st.markdown("- **85% - 95%**：寻找可能有显著亲和力提升的家族变体。\n- **> 95%**：通常是单点或双点突变，常用于寻找去除 PTM 风险的天然备用分子。")

# ================= 主界面 =================
st.title("🧬 抗体同源序列智能挖掘引擎 (Homology Miner)")
st.markdown("输入您的**优选种子克隆**，系统将从海量**候选序列库**中，为您深挖具有高同源性的潜在优秀变体分子。")

col1, col2 = st.columns(2)
with col1:
    st.subheader("🥇 第一步：输入优选分子 (Seed Clones)")
    seed_input = st.text_area("请在此粘贴优选分子的 FASTA 序列（支持输入多条）：", height=200, 
                              placeholder=">Seed_Clone_1\nEVQLVQSGAEVKKPGASVKVSCKASGYTFT...\n>Seed_Clone_2\nDIQMTQSPSSLSASVGDRVTITCRASQ...")

with col2:
    st.subheader("📚 第二步：输入候选序列库 (Database)")
    lib_input = st.text_area("请在此粘贴海量候选库的 FASTA 序列（将被用于遍历比对）：", height=200,
                             placeholder=">Library_Clone_001\nEVQLVQSGAEVKKPGASVKVSCKASGYTFT...\n...")

# ================= 数据处理与挖掘 =================
if st.button("🚀 开始极速同源挖掘", type="primary", use_container_width=True):
    if not seed_input.strip() or not lib_input.strip():
        st.error("⚠️ 请确保优选分子和候选序列库均已输入数据！")
    else:
        with st.spinner("正在构建序列空间并进行全局动态比对，请稍候..."):
            seeds = parse_fasta(seed_input)
            library = parse_fasta(lib_input)
            
            if not seeds or not library:
                st.error("⚠️ 序列解析失败，请检查 FASTA 格式是否规范。")
            else:
                st.success(f"解析成功！已加载 **{len(seeds)}** 个种子分子，**{len(library)}** 条候选序列库。")
                
                results = []
                progress_bar = st.progress(0)
                
                # 双重循环比对挖掘
                for i, (seed_name, seed_seq) in enumerate(seeds.items()):
                    seed_matches = []
                    for lib_name, lib_seq in library.items():
                        # 跳过完全相同的克隆（名字和序列都一样认为是同一个）
                        if seed_seq == lib_seq:
                            continue
                            
                        identity, mutations = calculate_similarity(seed_seq, lib_seq)
                        
                        if identity * 100 >= sim_threshold:
                            seed_matches.append({
                                '种子分子 (Seed)': seed_name,
                                '挖掘命中分子 (Hit)': lib_name,
                                '序列相似度 (%)': round(identity * 100, 2),
                                '氨基酸突变数': mutations,
                                '命中序列长度': len(lib_seq),
                                '命中序列 (Hit Sequence)': lib_seq
                            })
                    
                    # 按相似度降序排序，并截取 Top K
                    seed_matches = sorted(seed_matches, key=lambda x: x['序列相似度 (%)'], reverse=True)[:top_k]
                    results.extend(seed_matches)
                    
                    # 更新进度条
                    progress_bar.progress((i + 1) / len(seeds))
                
                # ================= 结果展示 =================
                if results:
                    df_results = pd.DataFrame(results)
                    
                    st.markdown("### 🎯 挖掘结果总览")
                    st.info(f"在 {sim_threshold}% 的相似度阈值下，共为您挖掘到 **{len(df_results)}** 条高潜力同源序列！")
                    
                    # 样式高亮函数
                    def highlight_sim(val):
                        if isinstance(val, float):
                            if val >= 95.0:
                                return 'background-color: #d4edda; color: #155724; font-weight: bold;'
                            elif val >= 90.0:
                                return 'background-color: #fff3cd; color: #856404;'
                        return ''
                    
                    st.dataframe(
                        df_results.style.map(highlight_sim, subset=['序列相似度 (%)']),
                        use_container_width=True,
                        height=400
                    )
                    
                    # ================= 导出功能 =================
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='openpyxl') as writer:
                        df_results.to_excel(writer, index=False, sheet_name='同源挖掘结果')
                    output.seek(0)
                    
                    st.download_button(
                        label="📥 导出同源挖掘报告 (Excel)",
                        data=output,
                        file_name="同源分子挖掘报告.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        type="primary"
                    )
                else:
                    st.warning("🧐 在当前阈值下，序列库中未能找到与种子分子足够相似的序列。您可以尝试在左侧侧边栏调低【最低同源相似度阈值】。")
