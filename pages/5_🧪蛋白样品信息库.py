import streamlit as st
import pandas as pd
import gspread
from datetime import datetime
import io
import csv

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
    """批量写入数据库（带有防覆盖绝对坐标写入逻辑）"""
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
        
        # ====== 修复覆盖问题的核心防线 ======
        # 获取当前表格真实存在的行数（包含表头）
        existing_rows = sheet.get_all_values()
        next_empty_row = len(existing_rows) + 1
        
        # 强制从 A 列的真实空白行开始定点写入，100% 杜绝覆盖
        sheet.update(f"A{next_empty_row}", data_to_insert)
        # ====================================
        
        return True, f"成功将 {len(data_to_insert)} 条样品记录安全写入云端数据库的第 {next_empty_row} 行！"
        
    except gspread.exceptions.SpreadsheetNotFound:
        return False, "找不到名为 Protein_Sample_DB 的表格，请检查名称和共享权限。"
    except Exception as e:
        return False, str(e)

# ================= 智能字段映射引擎 =================
def smart_extract_columns(raw_text):
    """
    全新 V3 提取引擎 (终极地狱级表格克星)：
    1. 纯手工构建矩阵：100% 杜绝因格式错乱导致的读取崩溃。
    2. 垂直压缩前5行：无视多级合并表头产生的换行和跨列。
    3. 右到左 (Right-to-Left) 扫描：绝对优先锁定最右侧的 Final 最终交付数据。
    4. 主键过滤法：精准剔除多级表头中因合并单元格产生的碎片化空行。
    """
    # 1. 完美处理 Excel 单元格内部的 Alt+Enter 换行符
    f = io.StringIO(raw_text.strip())
    reader = csv.reader(f, delimiter='\t', quotechar='"')
    data_matrix = list(reader)
    
    max_cols = max(len(row) for row in data_matrix) if data_matrix else 0
    for row in data_matrix:
        while len(row) < max_cols:
            row.append("")
            
    raw_df = pd.DataFrame(data_matrix, dtype=str)
    
    mapping = {}
    mapped_targets = set()
    
    num_cols = len(raw_df.columns)
    scan_depth = min(6, len(raw_df))
    
    # 2. 从右向左遍历列，确保右侧的 Final 列被优先记录
    for j in range(num_cols - 1, -1, -1):
        # 将前几行的文本垂直合并为一个超级字符串 blob
        blob = " ".join([str(raw_df.iat[i, j]).upper() for i in range(scan_depth)])
        
        target = None
        # 👇 在这里加入了 "分子名称" 和 "分子名" 等扩展词汇
        if "PROTEIN NAME" in blob or "项目名称" in blob or "蛋白名称" in blob or "分子名称" in blob or "分子名" in blob:
            target = "Protein Name"
        elif "LOT" in blob or "批号" in blob or "批次" in blob:
            target = "Lot No"
        elif "BUFFER" in blob or "缓冲液" in blob or "纯化步骤" in blob or "纯化方式" in blob:
            target = "Buffer"
        elif "CONCENTRATION" in blob or "CONC" in blob or "浓度" in blob:
            target = "Concentration(ug/ul)"
        elif ("TOTAL" in blob and "MG" in blob) or "AMOUNT" in blob or "总量" in blob or "总计" in blob:
            target = "Total(mg)"
        elif "YIELD" in blob or "TITER" in blob or "产量" in blob or "表达量" in blob:
            target = "Yield(mg/L)"
        elif "SEC-HPLC" in blob or ("PURITY" in blob and "SEC" in blob) or "PURITY" in blob or "纯度" in blob:
            target = "Purity by SEC-HPLC(%)"
        elif ("SDS-PAGE" in blob and "CE" not in blob):
            target = "Purity by SDS-PAGE under NR(%)"
        elif ("M.W" in blob and "REDUCING" not in blob) or "分子量" in blob or "KD" in blob or "MW" in blob:
            target = "M.W.(KDa)"
        elif "PI" in blob or "等电点" in blob:
            target = "PI"
            
        # 只有当该目标指标还没被记录过时才映射，防止被左侧的 AC 过程数据篡改
        if target and target not in mapped_targets:
            mapping[j] = target
            mapped_targets.add(target)
            
    if not mapping:
        return pd.DataFrame(), {}
        
    # 恢复正常的从左到右列序
    sorted_mapping = dict(sorted(mapping.items()))
    extracted_df = raw_df[list(sorted_mapping.keys())].copy()
    extracted_df.rename(columns=sorted_mapping, inplace=True)
    
    # 3. 寻找主键列来判断这一行是不是真实的蛋白数据
    check_col = "Protein Name" if "Protein Name" in extracted_df.columns else (
                "Lot No" if "Lot No" in extracted_df.columns else list(sorted_mapping.values())[0])
                
    def is_real_data(val):
        val_str = str(val).strip().upper()
        # 空白行或者由于合并单元格产生的 NaN 绝对不是数据
        if val_str in ["", "NAN", "NONE", "NAT"]: return False
        # 👇 同样在这里补充 "分子名称" 等过滤词，防止把表头当成数据写入
        header_keywords = ["项目名称", "PROTEIN NAME", "蛋白名称", "分子名称", "分子名", "LOT号", "LOT NO", "LOT", "批号"]
        if val_str in header_keywords: return False
        return True
        
    # 精准切除表头碎片行，留下纯净的数据
    clean_df = extracted_df[extracted_df[check_col].apply(is_real_data)].copy()
    clean_df.reset_index(drop=True, inplace=True)
    
    return clean_df, sorted_mapping

# ================= UI 界面 =================
st.title("🧪 蛋白样品自动入库系统 (Smart Paste)")
st.markdown("突破加密文件限制！只需从原始交付大表里 **全选并复制 (包含表头)**，粘贴到下方。系统将自动过滤无关列，精准提取 10 项核心 QC 记录并安全入库。")

pasted_data = st.text_area("📋 在此粘贴从 Excel 复制的数据 (请务必包含表头行)：", height=300, 
                           placeholder="请在您的 Excel 文件中，选中包含 Order ID, Protein Name, Lot No... 等所有内容的区域（包括表头），按 Ctrl+C，然后在此处 Ctrl+V。")

if pasted_data.strip():
    with st.spinner("🤖 正在穿透多级表头进行深度解析..."):
        try:
            # 启动 V3 智能提取引擎
            clean_df, matched_cols = smart_extract_columns(pasted_data)
            
            if clean_df.empty:
                st.warning("⚠️ 解析到的有效数据为空。请确保复制了真实的数据行，且包含核心指标列。")
            else:
                # 检查最核心的两个标识符是否找到
                if "Protein Name" not in clean_df.columns or "Lot No" not in clean_df.columns:
                    st.warning("⚠️ 未能识别到核心标识 (Protein Name 或 Lot No)，请检查复制的区域是否完整包含了表头行。")
                else:
                    st.success(f"✅ 解析成功！从极度错综复杂的格式中精准穿透，滤出 {len(matched_cols)} 列核心指标，共计 {len(clean_df)} 个样品。")
                    
                    st.markdown("### 👁️ 提取结果预览验证")
                    st.dataframe(clean_df, use_container_width=True)
                    
                    # 提交按钮
                    st.markdown("---")
                    if st.button("🚀 确认无误，一键批量追加至云端数据库 (Save to DB)", type="primary", use_container_width=True):
                        with st.spinner("正在将数据安全写入云端，请稍候..."):
                            success, msg = save_to_db(clean_df)
                            if success:
                                st.success(msg)
                                st.balloons()
                            else:
                                st.error(f"写入失败: {msg}")

        except Exception as e:
            st.error(f"解析或处理数据时发生意外错误: {str(e)}")
