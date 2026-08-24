import streamlit as st
import pandas as pd
import gspread
import io
from datetime import datetime

st.set_page_config(page_title="蛋白样品信息库", page_icon="🧪", layout="wide")

# 严格对齐图片中的 10 个字段名
HEADERS = [
    "Protein Name", 
    "Lot No", 
    "Buffer", 
    "Concentration(ug/ul)", 
    "Total(mg)", 
    "Yield(mg/L)", 
    "Purity by SEC-HPLC(%)", 
    "Purity by SDS-PAGE under NR(%)", 
    "M.W.(KDa)", 
    "PI",
    "Record Time" # 系统自动附加的录入时间
]

@st.cache_resource(ttl=3600)
def init_gspread():
    """纯净版初始化连接"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None, "未在后台 Secrets 检测到配置"
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client, "OK"
    except Exception as e:
        return None, str(e)

def load_sample_db():
    """读取云端样品库"""
    client, error_msg = init_gspread()
    if not client:
        return None, error_msg
    
    try:
        # 我们约定一个独立的新表格名字
        sheet = client.open("Protein_Sample_DB").sheet1
        data = sheet.get_all_records()
        if not data:
             return pd.DataFrame(columns=HEADERS), "OK"
        return pd.DataFrame(data), "OK"
    except gspread.exceptions.SpreadsheetNotFound:
        return None, "NOT_FOUND"
    except Exception as e:
        return None, str(e)

def save_to_sample_db(df_to_save):
    """批量追加数据到云端样品库"""
    client, error_msg = init_gspread()
    if not client:
        return False, error_msg
    try:
        sheet = client.open("Protein_Sample_DB").sheet1
        # 将 DataFrame 转换为二维列表并填充 NaN 为空字符串
        df_clean = df_to_save.fillna("")
        rows = df_clean.values.tolist()
        sheet.append_rows(rows)
        return True, "OK"
    except Exception as e:
        return False, str(e)

st.title("🧪 蛋白样品交付与质控信息库 (Protein Sample DB)")
st.markdown("结构化记录每一批次纯化交付的蛋白样品参数。支持批量从 Excel 复制粘贴，实现 QC 数据一键沉淀。")

# 加载数据
db_df, status = load_sample_db()

# 异常处理：如果没有找到表格
if status == "NOT_FOUND":
    st.error("🚨 无法连接到名为 `Protein_Sample_DB` 的云端表格！")
    st.info("""
    **请按照以下步骤初始化该数据库（仅需一次）：**
    1. 去您的 Google Drive 创建一个全新的空白表格，命名为 **`Protein_Sample_DB`** （一字不差）。
    2. 将第一行 (A1 到 K1) 依次填入以下表头：
       `Protein Name` | `Lot No` | `Buffer` | `Concentration(ug/ul)` | `Total(mg)` | `Yield(mg/L)` | `Purity by SEC-HPLC(%)` | `Purity by SDS-PAGE under NR(%)` | `M.W.(KDa)` | `PI` | `Record Time`
    3. 点击右上角“共享”，赋予机器人 `drive-bot@paper-downloader-491323.iam.gserviceaccount.com` 编辑者权限。
    4. 刷新本页面即可！
    """)
    st.stop()
elif status != "OK":
    st.error(f"⚠️ 数据库连接出错: {status}")
    st.stop()

tab_view, tab_add = st.tabs(["📊 样品数据总览", "➕ 批量/单条录入样品"])

with tab_view:
    # 状态概览卡片
    st.markdown("### 🗂️ 已入库样品概览")
    cols = st.columns(4)
    cols[0].metric("累计收录批次", f"{len(db_df)} 批")
    if not db_df.empty and 'Total(mg)' in db_df.columns:
        try:
            total_mg = pd.to_numeric(db_df['Total(mg)'], errors='coerce').sum()
            cols[1].metric("历史累计交付总量", f"{total_mg:.2f} mg")
        except:
            cols[1].metric("历史累计交付总量", "N/A")
            
        unique_proteins = db_df['Protein Name'].nunique() if 'Protein Name' in db_df.columns else 0
        cols[2].metric("独立蛋白种类", f"{unique_proteins} 种")
    
    st.divider()
    
    # 全局搜索
    search_query = st.text_input("🔍 搜索样品 (支持按 Protein Name, Lot No 等快速检索)：")
    
    if search_query and not db_df.empty:
        mask = db_df.apply(lambda row: row.astype(str).str.contains(search_query, case=False).any(), axis=1)
        display_df = db_df[mask]
    else:
        display_df = db_df
        
    st.dataframe(
        display_df,
        use_container_width=True,
        hide_index=True,
        height=500
    )
    
    if not display_df.empty:
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            display_df.to_excel(writer, index=False, sheet_name='Sample_DB')
        st.download_button(
            label="📥 导出当前数据为 Excel",
            data=output.getvalue(),
            file_name=f"蛋白样品库导出_{datetime.now().strftime('%Y%m%d')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

with tab_add:
    st.markdown("### 📝 从 Excel 批量粘贴入库")
    st.info("💡 **操作指南**：请在 Excel 中按照下方的字段顺序整理好数据（**不需要复制表头，只复制数据行**），然后直接粘贴到下方文本框中。")
    
    # 显示模板列以便用户对照
    st.code("数据列顺序要求：\n" + " | ".join(HEADERS[:-1]), language="text")
    
    paste_data = st.text_area("📋 在此粘贴 TSV/Excel 数据：", height=200, 
                              placeholder="F0630-m132\t20260818-20-1\tPBS pH 7.4\t2\t7.8\t780\t98.28%\t>90%\t146.19\t6.76\nF0630-m134\t20260818-20-2\tPBS pH 7.4\t2\t8\t800\t99.86%\t>90%\t144.1\t6.53")
    
    if st.button("👁️ 预览数据并校验", type="primary"):
        if paste_data.strip():
            try:
                # 尝试解析粘贴的数据 (Excel 复制默认是 Tab 分隔符)
                df_parsed = pd.read_csv(io.StringIO(paste_data), sep='\t', header=None)
                
                # 校验列数是否匹配 (减1是因为 Record Time 是系统自动生成的)
                expected_cols = len(HEADERS) - 1
                if len(df_parsed.columns) != expected_cols:
                    st.error(f"⚠️ 列数不匹配！系统期待 {expected_cols} 列，但您粘贴的数据解析出了 {len(df_parsed.columns)} 列。请检查是否有合并单元格或缺失列。")
                else:
                    # 分配标准表头
                    df_parsed.columns = HEADERS[:-1]
                    # 追加当前时间
                    df_parsed["Record Time"] = datetime.now().strftime("%Y-%m-%d %H:%M")
                    
                    st.success(f"✅ 成功解析 {len(df_parsed)} 条数据！请核对下方预览表：")
                    st.dataframe(df_parsed, use_container_width=True, hide_index=True)
                    
                    # 存入 session_state 以备提交
                    st.session_state.pending_upload = df_parsed
            except Exception as e:
                st.error(f"⚠️ 数据解析失败，请确保格式是纯文本制表符分隔: {e}")
        else:
            st.warning("请先粘贴数据。")

    if 'pending_upload' in st.session_state and not st.session_state.pending_upload.empty:
        if st.button("🚀 确认无误，批量写入云端数据库", use_container_width=True):
            with st.spinner("正在将样品数据同步至 Google Drive..."):
                success, msg = save_to_sample_db(st.session_state.pending_upload)
                if success:
                    st.success("🎉 数据已成功入库！您可以切换到『📊 样品数据总览』查看。")
                    # 清空暂存
                    st.session_state.pending_upload = pd.DataFrame() 
                else:
                    st.error(f"写入失败: {msg}")
