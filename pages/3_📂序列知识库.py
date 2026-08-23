import streamlit as st
import pandas as pd
import gspread

st.set_page_config(page_title="抗体序列云端知识库", page_icon="📂", layout="wide")

st.title("📂 抗体序列云端知识库 (Knowledge Base)")
st.markdown("连接至您的私人 Google Drive。您可以在这里检索历史入库的优质分子，进行资产沉淀与追踪。")

@st.cache_resource(ttl=3600)
def init_gspread():
    """纯净版初始化连接"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None, "未在后台 Secrets 检测到配置"
        
        # 极简官方读取法
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client, "OK"
    except Exception as e:
        return None, str(e)

# ================= 核心读取逻辑 =================
client, error_msg = init_gspread()

if not client:
    st.error(f"🤖 机器人启动失败，原因: {error_msg}")
else:
    try:
        # 🚀 重新回归按名称读取！
        # 机器人会在它自己的视野里，搜索名字精确匹配 "Antibody_Database" 的表格
        sheet = client.open("Antibody_Database").sheet1
        
        data = sheet.get_all_records()
        
        if not data:
             db_df = pd.DataFrame(columns=["克隆ID", "靶点", "入库时间", "重链序列", "轻链序列", "GRAVY", "Max_HIC_Score", "质控状态", "备注"])
        else:
             db_df = pd.DataFrame(data)

        st.success(f"✅ 成功连接至云端数据库！当前总入库序列数：**{len(db_df)}**")
        
        # ======== 数据面板 ========
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

    except gspread.exceptions.SpreadsheetNotFound:
        st.warning("⚠️ 机器人已经成功上线，但在它的视野里找不到名叫 `Antibody_Database` 的表格！")
        st.info("💡 终极排查建议：\n1. 检查表格名字是否多敲了空格。\n2. **必须将表格共享给这个机器人**：`drive-bot@paper-downloader-491323.iam.gserviceaccount.com`")
    except Exception as e:
        st.error(f"⚠️ 读取过程中发生未知错误: {e}")
