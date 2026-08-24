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
def smart_extract_columns(raw_df):
    """
    智能分析原始 Excel 的表头，自带中英双语词典，并解决重名列冲突（优先取靠左的最终数据）
    """
    mapping = {}
    mapped_targets = set() # 用于防止重复映射，确保重名列（如两个 Conc）只取第一个
    raw_columns = raw_df.columns.tolist()
    
    for col in raw_columns:
        # 清理表头名称中的换行符和多余空格，统一转大写
        col_str = str(col).replace('\n', ' ').strip().upper()
        target = None
        
        # 1. 蛋白名称
        if "PROTEIN NAME" in col_str or "项目名称" in col_str or "蛋白名称" in col_str:
            target = "Protein Name"
        # 2. 批号
        elif "LOT" in col_str or "批号" in col_str:
            target = "Lot No"
        # 3. 缓冲液
        elif "BUFFER" in col_str or "缓冲液" in col_str:
            target = "Buffer"
        # 4. 浓度
        elif "CONCENTRATION" in col_str or "CONC" in col_str or "浓度" in col_str:
            target = "Concentration(ug/ul)"
        # 5. 总量
        elif ("TOTAL" in col_str and "MG" in col_str) or "AMOUNT" in col_str or "总量" in col_str:
            target = "Total(mg)"
        # 6. 产量/滴度
        elif "YIELD" in col_str or "TITER" in col_str or "产量" in col_str or "表达量" in col_str:
            target = "Yield(mg/L)"
        # 7. SEC-HPLC 纯度 (如果只写了 Purity，默认先往 SEC 归类)
        elif "SEC-HPLC" in col_str or ("PURITY" in col_str and "SEC" in col_str) or "PURITY" in col_str or "纯度" in col_str:
            target = "Purity by SEC-HPLC(%)"
        # 8. SDS-PAGE 纯度
        elif ("SDS-PAGE" in col_str and "CE" not in col_str):
            target = "Purity by SDS-PAGE under NR(%)"
        # 9. 分子量
        elif ("M.W" in col_str and "REDUCING" not in col_str) or "分子量" in col_str or "KD" in col_str or "MW" in col_str:
            target = "M.W.(KDa)"
        # 10. PI
        elif col_str == "PI" or "等电点" in col_str:
            target = "PI"
            
        # 核心防撞车逻辑：只有当该目标字段还没有被映射过时，才进行映射。
        # 这样当遇到 Final 和 AC 下面两个同名的 Conc 时，左边的 Final 会被优先截获。
        if target and target not in mapped_targets:
            mapping[col] = target
            mapped_targets.add(target)
            
    # 根据映射关系，截取子表并重命名列
    extracted_df = raw_df[list(mapping.keys())].copy()
    extracted_df.rename(columns=mapping, inplace=True)
    
    return extracted_df, mapping

# ================= UI 界面 =================
st.title("🧪 蛋白样品自动入库系统 (Smart Paste)")
st.markdown("突破加密文件限制！只需从原始交付大表里 **全选并复制 (包含表头)**，粘贴到下方。系统内置双语视觉引擎，**无视多级合并表头与中文列名**，精准提取 10 项核心 QC 记录。")

pasted_data = st.text_area("📋 在此粘贴从 Excel 复制的数据 (连同乱七八糟的表头一起复制即可)：", height=300, 
                           placeholder="支持中英双语！选中包含 项目名称、Conc、purity 等区域，Ctrl+C 然后在此处 Ctrl+V。")

if pasted_data.strip():
    with st.spinner("🤖 正在智能解析粘贴的数据..."):
        try:
            # 🚀 极其强壮的“不规则抗干扰”解析逻辑（完美解决 Expected x fields saw y 报错）
            # 按换行符切分每一行
            lines = pasted_data.strip("\r\n").split('\n')
            # 按制表符切分每个单元格
            data_matrix = [line.split('\t') for line in lines]
            
            # 找到全表最长的一行（容错 Excel 复制时末尾多出的空白格子）
            max_cols = max(len(row) for row in data_matrix)
            
            # 补齐所有数据行的空白，使矩阵完美对齐，先不设列名
            for row in data_matrix:
                while len(row) < max_cols:
                    row.append("")
                    
            raw_df = pd.DataFrame(data_matrix)
            
            # 🚀 新增：雷达自动扫描“真实表头”所在行
            best_row_idx = 0
            max_matches = 0
            # 设置中英高频关键词雷达
            keywords = ['PROTEIN', 'LOT', 'CONC', 'AMOUNT', 'YIELD', 'PURITY', '项目', '分子', 'PI', '浓度', '纯度']
            
            # 扫描前 5 行，谁包含的关键词最多，谁就是真正的列名行
            for idx, row in raw_df.head(5).iterrows():
                row_str = " ".join([str(x).upper() for x in row.tolist()])
                matches = sum(1 for kw in keywords if kw in row_str)
                if matches > max_matches:
                    max_matches = matches
                    best_row_idx = idx
            
            # 设定真正的表头，并直接砍掉它上面没用的多级合并垃圾行
            raw_df.columns = raw_df.iloc[best_row_idx].astype(str)
            raw_df = raw_df.iloc[best_row_idx+1:].reset_index(drop=True)
            
            if raw_df.empty:
                st.warning("⚠️ 解析到的数据为空，请确保复制了真实的数据行。")
            else:
                # 启动智能提取引擎
                clean_df, matched_cols = smart_extract_columns(raw_df)
                
                # 检查最核心的两个标识符是否找到
                if "Protein Name" not in clean_df.columns or "Lot No" not in clean_df.columns:
                    st.warning("⚠️ 未能识别到核心标识 (Protein Name/项目名称 或 Lot No/批号)，请检查复制区域。")
                else:
                    st.success(f"✅ 解析成功！从极度杂乱的格式中精准过滤出 {len(matched_cols)} 列核心指标，共计 {len(clean_df)} 个样品。")
                    
                    st.markdown("### 👁️ 提取结果预览验证")
                    st.dataframe(clean_df, use_container_width=True)
                    
                    # 提交按钮
                    st.markdown("---")
                    if st.button("🚀 确认无误，一键批量追加至云端数据库 (Save to DB)", type="primary", use_container_width=True):
