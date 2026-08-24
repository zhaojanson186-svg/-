import streamlit as st
import pandas as pd
import re
import io

# 页面配置
st.set_page_config(page_title="蛋白理化性质分析小工具", page_icon="🧮", layout="wide")

# 1. 氨基酸单体平均分子量 (扣除水分子后的残基重量)
AA_MW = {
    'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
    'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
    'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
    'S': 87.0773, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
}
WATER_MW = 18.01524 # 游离 N 端和 C 端加起来多出的一个水分子

# 2. pKa 常数 (Henderson-Hasselbalch 方程使用)
PI_DICT = {
    'D': -3.90, 'E': -4.07, 'C': -8.18, 'Y': -10.46,
    'H': 6.04, 'K': 10.54, 'R': 12.48,
    'N_term': 8.0, 'C_term': -3.1
}

# 3. Kyte-Doolittle 亲水疏水性指数 (GRAVY)
HYDROPATHY_INDEX = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

@st.cache_data
def calculate_mw(seq):
    """计算分子量 (Da)"""
    if not seq: return 0.0
    mw = sum(AA_MW.get(aa, 0) for aa in seq) + WATER_MW
    return mw

@st.cache_data
def calculate_pi(seq):
    """计算理论等电点 (pI)"""
    if not seq: return 0.0
    counts = {aa: seq.count(aa) for aa in PI_DICT.keys() if len(aa) == 1}
    
    def net_charge(pH):
        charge = 0.0
        # 碱性带正电
        charge += counts.get('H', 0) / (1.0 + 10**(pH - PI_DICT['H']))
        charge += counts.get('K', 0) / (1.0 + 10**(pH - PI_DICT['K']))
        charge += counts.get('R', 0) / (1.0 + 10**(pH - PI_DICT['R']))
        charge += 1.0 / (1.0 + 10**(pH - PI_DICT['N_term']))
        # 酸性带负电
        charge -= counts.get('D', 0) / (1.0 + 10**(PI_DICT['D'] * -1 - pH))
        charge -= counts.get('E', 0) / (1.0 + 10**(PI_DICT['E'] * -1 - pH))
        charge -= counts.get('C', 0) / (1.0 + 10**(PI_DICT['C'] * -1 - pH))
        charge -= counts.get('Y', 0) / (1.0 + 10**(PI_DICT['Y'] * -1 - pH))
        charge -= 1.0 / (1.0 + 10**(PI_DICT['C_term'] * -1 - pH))
        return charge
        
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2.0
        if net_charge(mid) > 0: low = mid
        else: high = mid
    return (low + high) / 2.0

@st.cache_data
def calculate_extinction_coefficient(seq):
    """计算 280nm 处的消光系数 (M-1 cm-1) 和 0.1% 吸光度"""
    num_W = seq.count('W')
    num_Y = seq.count('Y')
    num_C = seq.count('C')
    
    # 公式一：所有 Cys 形成二硫键 (Oxidized)
    ext_ox = (num_W * 5500) + (num_Y * 1490) + (num_C * 125)
    # 公式二：所有 Cys 均处于还原态 (Reduced)
    ext_red = (num_W * 5500) + (num_Y * 1490)
    
    mw = calculate_mw(seq)
    
    # 0.1% (= 1 mg/ml) 吸光度计算: Abs = Extinction / MW
    abs_ox = ext_ox / mw if mw > 0 else 0
    abs_red = ext_red / mw if mw > 0 else 0
    
    return ext_ox, ext_red, abs_ox, abs_red

@st.cache_data
def calculate_gravy(seq):
    """计算总平均亲水性 (GRAVY)"""
    if not seq: return 0.0
    total_hydropathy = sum(HYDROPATHY_INDEX.get(aa, 0) for aa in seq)
    return total_hydropathy / len(seq)

@st.cache_data
def analyze_composition(seq):
    """分析氨基酸组分"""
    length = len(seq)
    composition = {}
    for aa in HYDROPATHY_INDEX.keys():
        count = seq.count(aa)
        composition[aa] = {
            "数量": count,
            "百分比 (%)": round((count / length) * 100, 2) if length > 0 else 0.0
        }
    df = pd.DataFrame.from_dict(composition, orient='index')
    df.index.name = "氨基酸"
    return df.sort_values(by="百分比 (%)", ascending=False)

def parse_input(raw_input):
    """解析输入，兼容 FASTA 格式和纯文本格式"""
    raw_input = raw_input.strip()
    sequences = {}
    if raw_input.startswith(">"):
        # FASTA 格式处理
        current_name = ""
        current_seq = []
        for line in raw_input.splitlines():
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(re.sub(r'[^A-Za-z]', '', line).upper())
        if current_name:
            sequences[current_name] = "".join(current_seq)
    else:
        # 纯文本处理
        clean_seq = re.sub(r'[^A-Za-z]', '', raw_input).upper()
        if clean_seq:
            sequences["Seq_1"] = clean_seq
            
    return sequences

st.title("🧮 蛋白理化性质分析计算器 (Protein Calculator)")
st.markdown("快速计算蛋白质或抗体片段的**分子量、等电点(pI)、消光系数、浓度换算(Abs 0.1%)** 及 **GRAVY 指数**。算法严格对标 ExPASy ProtParam 标准。")

with st.sidebar:
    st.header("💡 工具说明")
    st.markdown("""
    **计算公式基准：**
    *   **分子量 (MW)**: 平均同位素质量 (Average mass)。
    *   **等电点 (pI)**: ExPASy pKa 常数 & Henderson-Hasselbalch 迭代法。
    *   **消光系数 (ε)**: 以 280nm 为基准。
        *   氧化态：假定所有 Cys 参与形成二硫键。
        *   还原态：无二硫键存在。
    *   **Abs 0.1%**: 1 mg/ml 浓度下，1cm 光程的 280nm 吸光度值 (常用于 Nanodrop 定量)。
    """)
    st.divider()
    st.info("💡 **提示：** 支持直接粘贴纯氨基酸序列，也支持粘贴一条或多条标准的 FASTA 格式数据。")

input_text = st.text_area("✍️ 请在此粘贴蛋白质序列 (FASTA 或纯氨基酸序列)：", height=200, 
                          placeholder="例如:\n>My_Protein\nEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYWMHWVRQAPGKGLEWVSAIN...")

if st.button("🚀 开始极速分析", type="primary", use_container_width=True):
    if not input_text.strip():
        st.warning("⚠️ 请先输入氨基酸序列！")
    else:
        seqs = parse_input(input_text)
        if not seqs:
            st.error("⚠️ 无法识别到有效的氨基酸序列，请检查输入格式。")
        else:
            st.success(f"✅ 成功解析 {len(seqs)} 条序列，分析结果如下：")
            
            all_results = []
            
            for name, seq in seqs.items():
                length = len(seq)
                mw = calculate_mw(seq)
                pi = calculate_pi(seq)
                ext_ox, ext_red, abs_ox, abs_red = calculate_extinction_coefficient(seq)
                gravy = calculate_gravy(seq)
                
                # 保存至汇总列表用于导出
                all_results.append({
                    "分子名称": name,
                    "长度 (AA)": length,
                    "分子量 (Da)": round(mw, 2),
                    "分子量 (kDa)": round(mw / 1000, 2),
                    "理论 pI": round(pi, 2),
                    "GRAVY": round(gravy, 3),
                    "消光系数_氧化态 (M-1 cm-1)": ext_ox,
                    "Abs 0.1% (1 mg/ml)_氧化态": round(abs_ox, 3),
                    "消光系数_还原态 (M-1 cm-1)": ext_red,
                    "Abs 0.1% (1 mg/ml)_还原态": round(abs_red, 3),
                    "全序列": seq
                })
                
                with st.expander(f"🧬 分子报告: **{name}** (点击展开/折叠)", expanded=(len(seqs) == 1)):
                    # 核心指标卡片
                    cols = st.columns(5)
                    cols[0].metric("分子长度", f"{length} AA")
                    cols[1].metric("分子量 (MW)", f"{mw/1000:.2f} kDa")
                    cols[2].metric("理论等电点 (pI)", f"{pi:.2f}")
                    cols[3].metric("GRAVY 指数", f"{gravy:.3f}")
                    
                    # Cys 检测警告
                    num_cys = seq.count('C')
                    if num_cys % 2 != 0:
                        cols[4].metric("半胱氨酸 (Cys)", f"{num_cys} 个", "奇数风险 🚨", delta_color="inverse")
                    else:
                        cols[4].metric("半胱氨酸 (Cys)", f"{num_cys} 个", "偶数匹配 ✅", delta_color="normal")
                    
                    st.divider()
                    
                    # 浓度测定面板 (消光系数)
                    st.markdown("#### 🧪 光谱定量参数 (A280, 1 cm 光程)")
                    col_ox, col_red = st.columns(2)
                    with col_ox:
                        st.info(f"""**氧化态 (Oxidized / 含二硫键)** 
                        - 消光系数 (Extinction Coefficient): **{ext_ox}** M⁻¹ cm⁻¹
                        - **Abs 0.1% (= 1 mg/ml)**: **{abs_ox:.3f}**""")
                    with col_red:
                        st.warning(f"""**还原态 (Reduced / 游离巯基)**
                        - 消光系数 (Extinction Coefficient): **{ext_red}** M⁻¹ cm⁻¹
                        - **Abs 0.1% (= 1 mg/ml)**: **{abs_red:.3f}**""")
                    
                    # 氨基酸组分柱状图
                    st.markdown("#### 📊 氨基酸组成分析 (Amino Acid Composition)")
                    comp_df = analyze_composition(seq)
                    st.bar_chart(comp_df["百分比 (%)"], use_container_width=True)

            st.markdown("---")
            df_export = pd.DataFrame(all_results)
            st.subheader("📥 批量结果导出")
            st.dataframe(df_export.drop(columns=['全序列']), use_container_width=True)
            
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                df_export.to_excel(writer, index=False, sheet_name='理化性质分析')
            
            st.download_button(
                label="💾 下载 Excel 完整分析报告",
                data=output.getvalue(),
                file_name="蛋白理化性质分析报告.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary",
                use_container_width=True
            )
