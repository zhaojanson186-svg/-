import streamlit as st
import pandas as pd
import json
import gspread
from google.oauth2.credentials import Credentials
from google.auth.exceptions import RefreshError

st.set_page_config(page_title="抗体序列云端知识库", page_icon="📂", layout="wide")

# ================= 1. 核心通讯引擎 =================
@st.cache_resource(ttl=3600)  # 缓存连接，避免频繁请求断线
def init_gspread():
    """初始化并连接到 Google Sheets"""
    try:
        # 从 Streamlit Secrets 中安全读取密钥
        if "gcp_token" not in st.secrets:
            st.error("⚠️ 未在后台检测到 GCP 密钥配置 (Secrets)。")
            return None
            
        token_info = json.loads(st.secrets["gcp_token"]["oauth_json"])
        
        # 使用 OAuth 凭证构建通讯卡
        creds = Credentials.from_authorized_user_info(token_info)
        client = gspread.authorize(creds)
        return client
    except Exception as e:
        st.error(f"连接云端数据库失败: {e}")
        return None

def load_database():
    """读取指定的 Google Sheet 数据"""
    client = init_gspread()
    if client:
        try:
            # 尝试打开名为 Antibody_Database 的表格 (请确保您在网盘里新建了这个名字的表)
            sheet = client.open_by_key("1V8x9vN5mp0l2dw2FsCeUSMPdGB-zPDhTGHoE5OrEx82").sheet1
            data = sheet.get_all_records()
            if not data:
                return pd.DataFrame(columns=["克隆ID", "靶点", "入库时间", "重链序列", "轻链序列", "GRAVY", "Max_HIC_Score", "质控状态", "备注"])
            return pd.DataFrame(data)
        except gspread.exceptions.SpreadsheetNotFound:
            st.error("⚠️ 未在您的 Google Drive 中找到名为 `Antibody_Database` 的表格，请先创建。")
            return None
        except RefreshError:
            st.error("⚠️ 您的授权 Token 已过期，请重新生成并更新 Secrets。")
            return None
    return None

# ================= 2. 主界面 UI =================
st.title("📂 抗体序列云端知识库 (Knowledge Base)")
st.markdown("连接至您的私人 Google Drive。您可以在这里检索历史入库的优质分子，进行资产沉淀与追踪。")

db_df = load_database()

if db_df is not None:
    st.success(f"✅ 成功连接至云端数据库！当前总入库序列数：**{len(db_df)}**")
    
    # ======== 数据面板 ========
    st.markdown("### 📊 数据库全景")
    
    # 提供基础的交互式检索
    search_query = st.text_input("🔍 快速全局搜索 (支持克隆名 / 靶点 / 备注等)：")
    
    if search_query:
        # 模糊搜索所有列
        mask = db_df.apply(lambda row: row.astype(str).str.contains(search_query, case=False).any(), axis=1)
        display_df = db_df[mask]
    else:
        display_df = db_df
        
    # 渲染数据表 (加入勾选框，方便后续做批量操作)
    edited_df = st.data_editor(
        display_df,
        use_container_width=True,
        num_rows="dynamic",
        hide_index=True,
        disabled=["克隆ID", "入库时间", "重链序列", "轻链序列", "GRAVY", "Max_HIC_Score"] # 锁定核心数据不可篡改，但允许改备注
    )
    
    # ======== 拓展功能区 ========
    st.markdown("### 🛠️ 资产流转引擎")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.button("🧬 发送选中分子至同源挖掘引擎", disabled=True, help="该功能正在研发中...")
    with col2:
        st.button("🧪 发送选中分子重新进行 CMC 评估", disabled=True, help="该功能正在研发中...")
    with col3:
        if st.button("💾 同步修改到云端"):
            st.info("编辑同步功能暂未实装。当前请直接在 Google Sheets 中修改。")

else:
    st.warning("🔄 数据库离线或未初始化。")
