import streamlit as st
import pandas as pd
import gspread
from datetime import datetime
import io

st.set_page_config(page_title="蛋白样品智能入库系统", page_icon="🧪", layout="wide")

# ================= 数据库连接模块 =================
@st.cache_resource(ttl=3600)
def init_gspread():
    """初始化并连接到 Google Sheets"""
    try:
        if "gcp_service_account" not in st.secrets:
            return None, "未在后台检测到 GCP 服务账号配置"
        client = gspread.service_account_from_dict(st.secrets["gcp_service_account"])
        return client, "OK"
    except Exception as e:
        return None, str(e)

def save_to_db(df_to_save):
    """批量写入数据库"""
    client, error_msg = init_gspread()
    if not client:
        return False, f"数据库连接失败: {error_msg}"
    
    try:
        sheet = client.open("Protein_Sample_DB").sheet1
        
        # 将 DataFrame 转换为二维列表，处理 NaN 为 空字符串
        df_clean = df_to_save.fillna("")
        
        # 加上打卡时间戳
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        df_clean['Record Time'] = current_time
        
        # 确保列顺序绝对正确
        target_columns = [
            "Protein Name", "Lot No", "Buffer", "Concentration(ug/ul)", 
            "Total(mg)", "Yield(mg/L)", "Purity by SEC-HPLC(%)", 
            "Purity by SDS-PAGE under NR(%)", "M.W.(KDa)", "PI", "Record Time"
        ]
        
        # 补齐可能缺失的列
        for col in target_columns:
            if col not in df_clean.columns:
                df_clean[col] = ""
                
        df_final = df_clean[target_columns]
        data_to_insert = df_final.values.tolist()
        
        # 批量追加
        sheet.append_rows(data_to_insert)
        return True, f"成功将 {len(data_to_insert)} 条样品记录写入云端数据库！"
        
    except gspread.exceptions.SpreadsheetNotFound:
        return False, "找不到名为 Protein_Sample_DB 的表格，请检查名称和共享权限。"
    except Exception as e:
        return False, str(e)

# ================= 智能字段映射引擎 =================
def smart_extract_columns(raw_df):
    """
    智能分析原始 Excel 的表头，自动映射到我们的 10 个目标字段
    """
    mapping = {}
    raw_columns = raw_df.columns.tolist()
    
    for col in raw_columns:
        # 清理表头名称中的换行符和多余空格，统一转大写
        col_str = str(col).replace('\n', ' ').strip().upper()
        
        if "PROTEIN NAME" in col_str:
            mapping[col] = "Protein Name"
        elif "LOT NO" in col_str:
            mapping[col] = "Lot No"
        elif "BUFFER" in col_str:
            mapping[col] = "Buffer"
        elif "CONCENTRATION" in col_str:
            mapping[col] = "Concentration(ug/ul)"
        elif "TOTAL" in col_str and "MG" in col_str:
            mapping[col] = "Total(mg)"
        elif "YIELD" in col_str:
            mapping[col] = "Yield(mg/L)"
        elif "SEC-HPLC" in col_str:
            mapping[col] = "Purity by SEC-HPLC(%)"
        elif "SDS-PAGE" in col_str and "CE" not in col_str:
            mapping[col] = "Purity by SDS-PAGE under NR(%)"
        elif "M.W" in col_str and "REDUCING" not in col_str:
            mapping[col] = "M.W.(KDa)"
        elif col_str == "PI":
            mapping[col] = "PI"
            
    # 根据映射关系，截取子表并重命名列
    extracted_df = raw_df[list(mapping.keys())].copy()
    extracted_df.rename(columns=mapping, inplace=True)
    
    return extracted_df, mapping

# ================= UI 界面 =================
st.title("🧪 蛋白样品自动入库系统 (Smart Paste)")
st.markdown("突破加密文件限制！只需从原始交付大表里 **全选并复制 (包含表头)**，粘贴到下方。系统将自动过滤无关列，精准提取 10 项核心 QC 记录并安全入库。")

pasted_data = st.text_area("📋 在此粘贴从 Excel 复制的数据 (请务必包含表头行)：", height=300, 
                           placeholder="请在您的 Excel 文件中，选中包含 Order ID, Protein Name, Lot No... 等所有内容的区域（包括表头），按 Ctrl+C，然后在此处 Ctrl+V。")

if pasted_data.strip():
    with st.spinner("🤖 正在智能解析粘贴的数据..."):
        try:
            # 模拟 Excel 的复制逻辑：以制表符 (\t) 作为列的分隔符
            raw_df = pd.read_csv(io.StringIO(pasted_data), sep='\t')
            
            if raw_df.empty:
                st.warning("⚠️ 解析到的数据为空，请确保复制了数据行（不只是空行或纯文本）。")
            else:
                # 启动智能提取引擎
                clean_df, matched_cols = smart_extract_columns(raw_df)
                
                # 检查最核心的两个标识符是否找到
                if "Protein Name" not in clean_df.columns or "Lot No" not in clean_df.columns:
                    st.warning("⚠️ 未能识别到核心标识 (Protein Name 或 Lot No)，请检查复制的区域是否完整包含了表头行。")
                else:
                    st.success(f"✅ 解析成功！从 {len(raw_df.columns)} 列杂乱数据中，精准提取了 {len(matched_cols)} 列核心指标，共计 {len(clean_df)} 个样品。")
                    
                    st.markdown("### 👁️ 提取结果预览验证")
                    st.dataframe(clean_df, use_container_width=True)
                    
                    # 提交按钮
                    st.markdown("---")
                    if st.button("🚀 确认无误，一键批量入库 (Save to DB)", type="primary", use_container_width=True):
                        with st.spinner("正在安全同步至 Google Drive 知识库..."):
                            success, msg = save_to_db(clean_df)
                            if success:
                                st.balloons()
                                st.success(msg)
                            else:
                                st.error(msg)
                                
        except Exception as e:
            st.error(f"解析粘贴的内容时发生错误，请确保是从 Excel 表格直接复制的结构化数据。底层错误: {e}")
