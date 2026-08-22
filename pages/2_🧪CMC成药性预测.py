import streamlit as st
import pandas as pd
import re

st.set_page_config(page_title="CMC 成药性预警 | HIC与聚集扫描", page_icon="🧪", layout="wide")

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
    valid_aa = [aa for aa in seq if aa in KD_SCALE]
    return sum(KD_SCALE[aa] for aa in valid_aa) / len(valid_aa) if valid_aa else 0.0

def scan_aggregation_patches(seq, window_size, threshold):
    patches = []
    profile = []
    max_score = -99.0
    
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        score = sum(KD_SCALE.get(aa, 0) for aa in window) / window_size
        profile.append({"起始位点": i+1, "局部疏水得分": score, "窗口序列": window})
        
        if score > max_score:
            max_score = score
            
        if score >= threshold:
            patches.append({
                "起始位点": i+1, "结束位点": i+window_size,
                "高危序列": window, "平均疏水得分": round(score, 2)
            })
            
    merged_patches = []
    if patches:
        current_patch = patches[0].copy()
        for next_patch in patches[1:]:
            if next_patch["起始位点"] <= current_patch["结束位点"] + 1:
                current_patch["结束位点"] = next_patch["结束位点"]
                current_patch["高危序列"] = seq[current_patch["起始位点"]-1 : current_patch["结束位点"]]
                current_patch["平均疏水得分"] = max(current_patch["平均疏水得分"], next_patch["平均疏水得分"])
            else:
                merged_patches.append(current_patch)
                current_patch = next_patch.copy()
        merged_patches.append(current_patch)
        
    return merged_patches, pd.DataFrame(profile), max_score

with st.sidebar:
    st.header("⚙️ 雷达参数设置")
    window_size = st.slider("滑动窗口大小", 5, 15, 7, 1)
    apr_threshold = st.slider("聚集热点报警阈值", 0.5, 3.0, 1.5, 0.1)

st.title("🧪 CMC 成药性预警 | HIC与聚集扫描")
st.markdown("不仅评估整体溶解度 (GRAVY)，更通过提取 **表面高危疏水斑块 (Max Patch Score)**，精准预测抗体在 **HIC (疏水色谱)** 中的异常滞留与高粘度风险。")

seq_input = st.text_area("请在此粘贴序列 (FASTA格式)：", height=150)

if st.button("🔍 开始成药性风险扫描", type="primary", use_container_width=True):
    if not seq_input.strip():
        st.error("⚠️ 请先输入抗体序列！")
    else:
        sequences = parse_fasta(seq_input)
        for name, seq in sequences.items():
            st.markdown("---")
            st.subheader(f"🧬 分子: `{name}`")
            
            gravy = calculate_gravy(seq)
            aprs, profile_df, max_score = scan_aggregation_patches(seq, window_size, apr_threshold)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                if gravy < -0.2:
                    st.metric("整体 GRAVY (预测溶解度)", f"{gravy:.3f}", "🟢 整体亲水", delta_color="normal")
                else:
                    st.metric("整体 GRAVY (预测溶解度)", f"{gravy:.3f}", "🔴 极高沉淀风险", delta_color="inverse")
            
            with col2:
                # ================= 新增：专属 HIC 预测指标 =================
                if max_score > 1.8:
                    st.metric("最高局部疏水峰 (HIC 预测)", f"{max_score:.2f}", "🔴 极端 HIC 滞留/高粘度", delta_color="inverse")
                elif max_score > 1.5:
                    st.metric("最高局部疏水峰 (HIC 预测)", f"{max_score:.2f}", "🟡 中高 HIC 风险", delta_color="off")
                else:
                    st.metric("最高局部疏水峰 (HIC 预测)", f"{max_score:.2f}", "🟢 HIC 表现良好", delta_color="normal")
                
            with col3:
                apr_count = len(aprs)
                if apr_count == 0:
                    st.metric("高危斑块数量 (APRs)", "0", "🟢 未见聚集热点", delta_color="normal")
                else:
                    st.metric("高危斑块数量 (APRs)", f"{apr_count} 处", "🔴 存在局部聚集", delta_color="inverse")

            if aprs:
                st.error("⚠️ **HIC/粘度高危预警**：发现以下极度疏水的连续斑块。这些区域如果暴露在 3D 表面，将导致 HIC 晚出峰（如 >20min）以及极高的制剂粘度！")
                st.table(pd.DataFrame(aprs))
            
            st.markdown("#### 🌊 动态疏水性全景图谱")
            chart_data = profile_df.set_index("起始位点")[["局部疏水得分"]]
            chart_data["高危预警线"] = apr_threshold
            st.line_chart(chart_data, color=["#1f77b4", "#d62728"], height=300)
