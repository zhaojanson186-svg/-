import streamlit as st
import pandas as pd
import gspread

st.set_page_config(page_title="全景数据追踪 | 干湿实验闭环", page_icon="📊", layout="wide")

# ================= 数据库连接与加载模块 =================
@st.cache_resource(ttl=3600)
def init_gspread():
    """初始化并连接到 Google Sheets"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None, "未检测到 GCP 服务账号配置"
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client, "OK"
    except Exception as e:
        return None, str(e)

@st.cache_data(ttl=300) # 缓存5分钟，保证数据相对新鲜
def load_all_data():
    client, msg = init_gspread()
    if not client:
        return None, None, f"数据库连接失败: {msg}"
        
    try:
        # 1. 读取序列知识库 (Dry Lab)
        try:
            sheet_seq = client.open("Antibody_Database").sheet1
            data_seq = sheet_seq.get_all_records()
            df_seq = pd.DataFrame(data_seq) if data_seq else pd.DataFrame()
        except Exception:
            df_seq = pd.DataFrame()

        # 2. 读取样品信息库 (Wet Lab)
        try:
            sheet_sample = client.open("Protein_Sample_DB").sheet1
            data_sample = sheet_sample.get_all_records()
            df_sample = pd.DataFrame(data_sample) if data_sample else pd.DataFrame()
        except Exception:
            df_sample = pd.DataFrame()
            
        return df_seq, df_sample, "OK"
    except Exception as e:
        return None, None, f"读取数据时发生错误: {str(e)}"

# ================= 数据清洗与融合引警 =================
import re

def merge_databases(df_seq, df_sample):
    if df_seq.empty or df_sample.empty:
        return pd.DataFrame()
        
    # 为了保证匹配的鲁棒性，去除空格并统一大写进行准备
    df_seq['Clean_ID'] = df_seq['克隆ID'].astype(str).str.strip().str.upper()
    df_sample['Clean_Name'] = df_sample['Protein Name'].astype(str).str.strip().str.upper()
    
    merged_records = []
    
    # 执行智能模糊交叉比对
    for _, seq_row in df_seq.iterrows():
        seq_name = seq_row['Clean_ID']
        if not seq_name or seq_name in ['NAN', 'NONE', '']: continue
        
        for _, sample_row in df_sample.iterrows():
            sample_name = sample_row['Clean_Name']
            if not sample_name or sample_name in ['NAN', 'NONE', '']: continue
            
            # 智能模糊匹配：核心名称互相包含即可
            # 找出短名字和长名字
            short_n, long_n = (seq_name, sample_name) if len(seq_name) < len(sample_name) else (sample_name, seq_name)
            
            # 💡 核心防御正则：短名称的前后不能紧挨着其他字母或数字（必须是边界或如 - _ 等连接符）
            # 这能完美防止 "M1" 错误匹配到 "M12"
            pattern = r'(?:^|[^a-zA-Z0-9])' + re.escape(short_n) + r'(?:[^a-zA-Z0-9]|$)'
            
            if re.search(pattern, long_n):
                # 匹配成功，将两张表的这一行数据合并为一条大记录
                combined = {**seq_row.to_dict(), **sample_row.to_dict()}
                # 设立一个统一的 Match_Key 用于前端展示和去重计数
                combined['Match_Key'] = seq_row['克隆ID']
                merged_records.append(combined)
                
    df_merged = pd.DataFrame(merged_records)
    
    # 按照表达时间倒序排列，优先看最新批次
    if not df_merged.empty and 'Record Time' in df_merged.columns:
        df_merged.sort_values(by='Record Time', ascending=False, inplace=True)
        
    return df_merged

# ================= 页面 UI 与交互 =================
st.title("📊 研发全景数据追踪 (Dry-Wet Loop)")
st.markdown("这里是整个中台的神经中枢。系统会自动将 **序列生信库 (Dry Lab)** 与 **QC 样品库 (Wet Lab)** 打通，为您提供同一分子的多维全景视图，实现成药性数据的干湿闭环验证。")

with st.spinner("🔗 正在从云端并发拉取两大数据库，并进行智能关联融合..."):
    df_seq, df_sample, status = load_all_data()

if status != "OK":
    st.error(status)
    st.stop()

# 顶部全局大盘数据
col1, col2, col3 = st.columns(3)
col1.metric("🧬 序列库资产总量 (Dry Lab)", f"{len(df_seq)} 条" if not df_seq.empty else "0 条")
col2.metric("🧪 样品批次总量 (Wet Lab)", f"{len(df_sample)} 批" if not df_sample.empty else "0 批")

df_merged = merge_databases(df_seq, df_sample)
unique_merged_clones = df_merged['Match_Key'].nunique() if not df_merged.empty else 0
col3.metric("🤝 成功闭环分子数 (已印证)", f"{unique_merged_clones} 个", "高价值数据资产")

st.divider()

if df_merged.empty:
    st.info("💡 暂无成功关联的闭环数据。请确保在《序列库》中存在的克隆，在《样品库》中也录入了名称完全一致 (无多余空格) 的表达记录。")
    st.stop()

st.subheader("🔍 分子全景档案检索")
search_query = st.selectbox(
    "请选择或输入要查看的闭环分子档案：", 
    options=df_merged['克隆ID'].unique(),
    help="下拉列表中的分子，均是既有序列记录、又有表达QC记录的闭环资产。"
)

if search_query:
    # 抽取选定分子的所有数据 (可能有一个序列对应多个表达批次的情况)
    target_data = df_merged[df_merged['克隆ID'] == search_query]
    seq_info = target_data.iloc[0] # 序列信息取第一条即可，因为是同一个克隆
    
    st.markdown(f"### 🪪 档案：{search_query}")
    st.caption(f"🎯 对应靶点：**{seq_info.get('靶点', '未登记')}** | 首次入库时间：{seq_info.get('入库时间', '未知')}")
    
    # 核心干湿对比视图
    col_dry, col_wet = st.columns(2)
    
    # ============ 👈 左侧：生信预测卡片 ============
    with col_dry:
        st.markdown("#### 🧬 In Silico (生信预测)")
        st.info(f"**生信质控状态：**\n{seq_info.get('质控状态', 'N/A')}")
        
        # 指标展示
        m1, m2 = st.columns(2)
        m1.metric("理论 GRAVY (疏水性)", seq_info.get('GRAVY', 'N/A'))
        m2.metric("Max HIC Score (成药风险)", seq_info.get('Max_HIC_Score', 'N/A'))
        
        with st.expander("👁️ 查看底层序列"):
            st.code(f">VH\n{seq_info.get('重链序列', 'N/A')}\n>VL\n{seq_info.get('轻链序列', 'N/A')}")
            st.caption(f"生信备注：{seq_info.get('备注_x', '无')}") # merge后左表备注为备注_x
            
    # ============ 👉 右侧：湿实验真相卡片 ============
    with col_wet:
        st.markdown("#### 🧪 In Vitro (湿实验真相)")
        batch_count = len(target_data)
        st.success(f"共匹配到 **{batch_count}** 个实际表达与纯化批次。")
        
        # 提取历史批次的最佳数据进行展示
        best_yield = target_data['Yield(mg/L)'].max()
        best_sec = target_data['Purity by SEC-HPLC(%)'].max()
        
        m3, m4 = st.columns(2)
        m3.metric("最高表达量 (Yield)", f"{best_yield} mg/L")
        m4.metric("最高 SEC 纯度", f"{best_sec}")
        
        with st.expander("📦 查看各批次 QC 详情", expanded=True):
            # 过滤出展示列，避免过于臃肿
            display_cols = ['Lot No', 'Yield(mg/L)', 'Purity by SEC-HPLC(%)', 'Concentration(ug/ul)', 'PI']
            display_cols = [c for c in display_cols if c in target_data.columns]
            
            st.dataframe(target_data[display_cols], hide_index=True, use_container_width=True)

st.markdown("---")
st.subheader("🌐 全景闭环数据底座 (Global Joined View)")
st.markdown("这里汇聚了系统内所有达成干湿闭环的分子。您可以导出该表，用于**训练 AI 序列预测大模型**或进行宏观的大规模成药性相关分析。")

# 整理最终呈现列
cols_to_show = [
    '克隆ID', '靶点', 'Lot No', 'Yield(mg/L)', 'Purity by SEC-HPLC(%)', 
    '质控状态', 'GRAVY', 'PI', '重链序列', '轻链序列'
]
cols_to_show = [c for c in cols_to_show if c in df_merged.columns]

st.dataframe(df_merged[cols_to_show], use_container_width=True)
