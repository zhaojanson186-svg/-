import streamlit as st
import pandas as pd
import re

# 页面配置
st.set_page_config(page_title="CMC 成药性预警 | 疏水聚集扫描", page_icon="🧪", layout="wide")

# ================= 科学参数与算法基座 =================
# Kyte-Doolittle 氨基酸疏水性指数 (正值越大约疏水，负值越大约亲水)
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

@st.cache_data
def parse_fasta(fasta_text):
    sequences = {}
    current_name = ""
    current_seq = []
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if current_name: sequences[current_name] = "".join(current_seq).upper()
            current_name = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(re.sub(r'[^A-Za-z]', '', line))
    if current_name: sequences[current_name] = "".join(current_seq).upper()
    return sequences

def calculate_gravy(seq):
    """计算总平均疏水性 (GRAVY)"""
    valid_aa = [aa for aa in seq if aa in KD_SCALE]
    if not valid_aa: return 0.0
    return sum(KD_SCALE[aa] for aa in valid_aa) / len(valid_aa)

def scan_aggregation_patches(seq, window_size, threshold):
    """滑动窗口扫描局部聚集热点 (Aggregation-Prone Regions, APRs)"""
    patches = []
    profile = []
    
    # 构建全序列疏水性图谱
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        score = sum(KD_SCALE.get(aa, 0) for aa in window) / window_size
        profile.append({"起始位点": i+1, "局部疏水得分": score, "窗口序列": window})
        
        # 抓取超过高危阈值的片段
        if score >= threshold:
            patches.append({
                "起始位点": i+1,
                "结束位点": i+window_size,
                "高危序列": window,
                "平均疏水得分": round(score, 2)
            })
            
    # 合并相邻或重叠的高危片段，避免重复报警
    merged_patches = []
    if patches:
        current_patch = patches[0].copy()
        for next_patch in patches[1:]:
            if next_patch["起始位点"] <= current_patch["结束位点"] + 1:
                # 发生重叠，扩展当前片段
                current_patch["结束位点"] = next_patch["结束位点"]
                current_patch["高危序列"] = seq[current_patch["起始位点"]-1 : current_patch["结束位点"]]
                current_patch["平均疏水得分"] = max(current_patch["平均疏水得分"], next_patch["平均疏水得分"])
            else:
                merged_patches.append(current_patch)
                current_patch = next_patch.copy()
        merged_patches.append(current_patch)
        
    return merged_patches, pd.DataFrame(profile)

# ================= 侧边栏与系统参数 =================
with st.sidebar:
    st.header("⚙️ 疏水雷达参数设置")
    st.info("基于 **Kyte-Doolittle** 物理化学量表进行预测。")
    window_size = st.slider("滑动窗口大小 (Window Size)", min_value=5, max_value=15, value=7, step=1,
                            help="建议 7-9。窗口过小噪音多，过大可能漏掉短距疏水斑块。")
    apr_threshold = st.slider("聚集热点报警阈值", min_value=0.5, max_value=3.0, value=1.5, step=0.1,
                              help="KD得分超过此值的片段将被标记为高危易聚集区 (APR)。")

# ================= 主界面 UI =================
st.title("🧪 CMC 成药性预警 | 疏水聚集扫描 (Developability)")
st.markdown("通过计算整体亲水/疏水平衡 (**GRAVY**)，并利用滑动窗口揪出导致多聚体和高粘度的 **局部聚集热点 (APRs)**。")

seq_input = st.text_area("请在此粘贴待进行 CMC 评估的抗体单链序列 (FASTA格式)：", height=150, 
                         placeholder=">S0728-CC2-B08_VH\nEVQLQESGPELVKPGTSVKISCKASGYPFTDYYINWVKQRPGQGLEWIGRIFPGSGSTYYNAKFMVKATLTVDKSSSTAYMLLSRLTSEDSAVYFCARIYDGYYPSDYWGQGTTLTVSS")

if st.button("🔍 开始成药性风险扫描", type="primary", use_container_width=True):
    if not seq_input.strip():
        st.error("⚠️ 请先输入抗体序列！")
    else:
        sequences = parse_fasta(seq_input)
        
        for name, seq in sequences.items():
            st.markdown("---")
            st.subheader(f"🧬 分子: `{name}`")
            
            # 计算指标
            gravy = calculate_gravy(seq)
            aprs, profile_df = scan_aggregation_patches(seq, window_size, apr_threshold)
            
            # 1. 宏观指标展示
            col1, col2, col3 = st.columns(3)
            with col1:
                # GRAVY 评估逻辑：抗体通常是水溶性蛋白，GRAVY应为负值。正值极大增加沉淀风险。
                if gravy < -0.2:
                    st.metric("总平均疏水性 (GRAVY)", f"{gravy:.3f}", "🟢 亲水性良好", delta_color="normal")
                elif gravy < 0:
                    st.metric("总平均疏水性 (GRAVY)", f"{gravy:.3f}", "🟡 轻微疏水倾向", delta_color="off")
                else:
                    st.metric("总平均疏水性 (GRAVY)", f"{gravy:.3f}", "🔴 极高沉淀风险", delta_color="inverse")
            
            with col2:
                st.metric("序列总长度", f"{len(seq)} AA")
                
            with col3:
                apr_count = len(aprs)
                if apr_count == 0:
                    st.metric("局部聚集热点 (APRs)", "0", "🟢 未发现高危斑块", delta_color="normal")
                else:
                    st.metric("局部聚集热点 (APRs)", f"{apr_count} 处", "🔴 存在聚集风险", delta_color="inverse")

            # 2. 聚集斑块具体报警信息
            if aprs:
                st.error(f"**高危预警**：雷达扫描发现 **{len(aprs)}** 处极度疏水的连续斑块，下游极易诱发多聚体或不溶！建议重点排查是否位于表面或 CDR 区。")
                apr_df = pd.DataFrame(aprs)
                st.table(apr_df)
            else:
                st.success("**CMC 质控通过**：未在当前阈值下扫描到高危连续疏水斑块。")

            # 3. 疏水性动态图谱
            st.markdown("#### 🌊 动态疏水性全景图谱")
            st.markdown("曲线向上（正值区）代表疏水，曲线向下（负值区）代表亲水。横线为设置的高危预警线。")
            
            # 使用 Streamlit 自带的轻量级折线图 (无需安装新依赖)
            # 整理画图数据
            chart_data = profile_df.set_index("起始位点")[["局部疏水得分"]]
            chart_data["高危预警线"] = apr_threshold
            
            st.line_chart(chart_data, color=["#1f77b4", "#d62728"], height=300)
