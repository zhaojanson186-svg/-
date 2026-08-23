import streamlit as st
import pandas as pd
from streamlit_gsheets import GSheetsConnection

st.set_page_config(page_title="抗体序列云端知识库", page_icon="📂", layout="wide")

st.title("📂 抗体序列云端知识库 (Knowledge Base)")
st.markdown("连接至您的私人 Google Drive。您可以在这里检索历史入库的优质分子，进行资产沉淀与追踪。")

# ================= 1. 使用官方连接器极简读取 =================
try:
    # 直接初始化官方连接器
    conn = st.connection("gsheets", type=GSheetsConnection)
    
    # ⚠️ 请在这里填入您那个 100% 准确的表格公开链接！
    # （只需要填下面这个 URL，不需要再拆分 ID 了）
    sheet_url = "https://docs.google.com/spreadsheets/d/1V8x9vN5mp0l2dw2FsCeUSMPdGB-zPDhTGHoE5OrEx82/edit"
    
    # 直接读取数据
    db_df = conn.read(spreadsheet=sheet_url, ttl="10m")
    
    if db_df is None or db_df.empty:
         # 如果表格存在但是空的，提供标准表头
         db_df = pd.DataFrame(columns=["克隆ID", "靶点", "入库时间", "重链序列", "轻链序列", "GRAVY", "Max_HIC_Score", "质控状态", "备注"])

    st.success(f"✅ 成功连接至云端数据库！当前总入库序列数：**{len(db_df)}**")
    
    # ======== 2. 数据面板 ========
    st.markdown("### 📊 数据库全景")
    
    search_query = st.text_input("🔍 快速全局搜索 (支持克隆名 / 靶点 / 备注等)：")
    
    if search_query:
        mask = db_df.apply(lambda row: row.astype(str).str.contains(search_query, case=False).any(), axis=1)
        display_df = db_df[mask]
    else:
        display_df = db_df
        
    edited_df = st.data_editor(
        display_df,
        use_container_width=True,
        num_rows="dynamic",
        hide_index=True
    )

except Exception as e:
    st.error(f"⚠️ 读取表格失败，底层报错: {e}")
