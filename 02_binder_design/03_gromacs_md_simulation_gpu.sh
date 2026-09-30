#!/bin/bash
# -----------------------------------------------------------------------------
# GROMACS 分子动力学模拟自动化脚本 (GPU 加速版)
# 针对 project_root/protein_name/{input,workspace,output} 目录结构优化
# -----------------------------------------------------------------------------

# 任何命令失败时立即退出
set -e

echo "--- GROMACS 模拟自动化脚本 (GPU版) 开始执行 ---"
echo "启动时间: $(date)"
echo ""

# -----------------------------------------------------------------------------
# 1. 变量定义和环境设置
# -----------------------------------------------------------------------------
# ---------- 用户环境配置 ----------
WORK_ROOT="/home/你的工作目录"
CONDA_BASE="/home/你的Conda安装目录"
GMX_RC="/opt/gromacsgpu/bin/GMXRC"

# 项目和蛋白名称
export PROJECT_ROOT_DIR="$WORK_ROOT/binder260906/MD"
export PROTEIN_NAME="001"  # 001 为示例编号，具体以实际项目为准
# ----------------------------------

# 定义子目录路径
export PROTEIN_DIR="${PROJECT_ROOT_DIR}/${PROTEIN_NAME}"
export INPUT_DIR="${PROTEIN_DIR}/input"
export WORKSPACE_DIR="${PROTEIN_DIR}/workspace"
export OUTPUT_DIR="${PROTEIN_DIR}/output"

# 输入 PDB 文件路径
export PDB_INPUT="${INPUT_DIR}/001_624_dldesign_0_cycle1_af2pred.pdb"

# MDP 文件目录
export MDP_DIR="${PROJECT_ROOT_DIR}/gromacs_mdp"

# 确保所有目录存在
mkdir -p "$INPUT_DIR" "$WORKSPACE_DIR" "$OUTPUT_DIR"

# 进入 workspace 目录
cd "$WORKSPACE_DIR" || { echo "错误: 无法进入工作目录 $WORKSPACE_DIR"; exit 1; }
echo "--- 当前工作目录: $(pwd) ---"

# --- GPU 加速设置 开始 ---
echo "--- 加载 Conda 环境以支持 GROMACS GPU 版本 ---"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate gromacs
echo "--- Conda 环境 'gromacs' 已激活 ---"
# --- GPU 加速设置 结束 ---

# 加载 GROMACS 环境
# 如果 Conda 环境已自动配置 GROMACS，可将下一行改为: GMX_RC=""
if [[ -n "$GMX_RC" ]]; then
    source "$GMX_RC" || { echo "错误: 无法 source GMXRC 文件，请检查 $GMX_RC 或 Conda 环境配置"; exit 1; }
fi

# 检查GROMACS版本
echo "GROMACS 版本信息:"
gmx -version | head -n 3
echo "--- 环境设置完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 2. 蛋白质预处理 (pdb2gmx)
# -----------------------------------------------------------------------------
echo "--- 步骤 2: 蛋白质预处理 (pdb2gmx) ---"
gmx pdb2gmx -f "$PDB_INPUT" \
 -ff "charmm27" \
 -water tip3p \
 -ignh \
 -p "${PROTEIN_NAME}.top" \
 -o "${PROTEIN_NAME}_processed.gro" \
 || { echo "错误: pdb2gmx 步骤失败"; exit 1; }
echo "--- pdb2gmx 完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 3. 定义模拟盒子 (editconf)
# -----------------------------------------------------------------------------
echo "--- 步骤 3: 定义模拟盒子 (editconf) ---"
gmx editconf -f "${PROTEIN_NAME}_processed.gro" \
 -o "${PROTEIN_NAME}_box.gro" \
 -bt dodecahedron \
 -d 3.0 \
 || { echo "错误: editconf 步骤失败"; exit 1; }
echo "--- editconf 完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 4. 溶剂化 (solvate)
# -----------------------------------------------------------------------------
echo "--- 步骤 4: 溶剂化 (solvate) ---"
gmx solvate -cp "${PROTEIN_NAME}_box.gro" \
 -cs spc216.gro \
 -p "${PROTEIN_NAME}.top" \
 -o "${PROTEIN_NAME}_solv.gro" \
 || { echo "错误: solvate 步骤失败"; exit 1; }
echo "--- solvate 完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 5. 添加离子 (grompp for ions & genion)
# -----------------------------------------------------------------------------
echo "--- 步骤 5: 添加离子 ---"
gmx grompp -f "${MDP_DIR}/ions.mdp" \
 -c "${PROTEIN_NAME}_solv.gro" \
 -p "${PROTEIN_NAME}.top" \
 -o ions.tpr \
 -maxwarn 1 \
 || { echo "错误: ions grompp 步骤失败"; exit 1; }

echo "SOL" | \
gmx genion -s ions.tpr \
 -p "${PROTEIN_NAME}.top" \
 -o "${PROTEIN_NAME}_ions.gro" \
 -neutral \
 -pname NA \
 -nname CL \
 || { echo "错误: genion 步骤失败"; exit 1; }
echo "--- 添加离子完成 ---"
echo ""

# OpenMP 线程数；可根据 CPU/GPU 配置调整。
export NT_THREADS="8"

# -----------------------------------------------------------------------------
# 6. 能量最小化 (em)
# -----------------------------------------------------------------------------
echo "--- 步骤 6: 能量最小化 (GPU) ---"
gmx grompp -f "${MDP_DIR}/em.mdp" \
 -c "${PROTEIN_NAME}_ions.gro" \
 -p "${PROTEIN_NAME}.top" \
 -o em.tpr \
 -maxwarn 1 \
 || { echo "错误: em grompp 步骤失败"; exit 1; }

gmx mdrun -v -deffnm em \
 -nt "${NT_THREADS}" \
 -nb gpu \
 || { echo "错误: em mdrun 步骤失败"; exit 1; }
echo "--- 能量最小化完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 7. NVT 平衡
# -----------------------------------------------------------------------------
echo "--- 步骤 7: NVT 平衡 (GPU) ---"
gmx grompp -f "${MDP_DIR}/nvt.mdp" \
 -c em.gro \
 -r em.gro \
 -p "${PROTEIN_NAME}.top" \
 -o nvt.tpr \
 -maxwarn 1 \
 || { echo "错误: nvt grompp 步骤失败"; exit 1; }

gmx mdrun -v -deffnm nvt \
 -nt "${NT_THREADS}" \
 -nb gpu \
 -pme gpu \
 || { echo "错误: nvt mdrun 步骤失败"; exit 1; }
echo "--- NVT 平衡完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 8. NPT 平衡
# -----------------------------------------------------------------------------
echo "--- 步骤 8: NPT 平衡 (GPU) ---"
gmx grompp -f "${MDP_DIR}/npt.mdp" \
 -c nvt.gro \
 -r nvt.gro \
 -t nvt.cpt \
 -p "${PROTEIN_NAME}.top" \
 -o npt.tpr \
 -maxwarn 1 \
 || { echo "错误: npt grompp 步骤失败"; exit 1; }

gmx mdrun -v -deffnm npt \
 -nt "${NT_THREADS}" \
 -nb gpu \
 -pme gpu \
 || { echo "错误: npt mdrun 步骤失败"; exit 1; }
echo "--- NPT 平衡完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 9. 正式模拟 (md)
# -----------------------------------------------------------------------------
echo "--- 步骤 9: 正式模拟 (GPU) ---"
echo "--- 创建 index 文件 ---"
if [ ! -f npt.gro ]; then
    echo "错误: npt.gro 文件未找到，无法创建 index 文件。"
    exit 1
fi
gmx make_ndx -f npt.gro -o index.ndx <<EOF
0
q
EOF
if [ ! -f index.ndx ]; then
    echo "错误: index.ndx 文件创建失败。"
    exit 1
fi
echo "--- index 文件创建完成 ---"

echo "--- 准备 md.tpr 文件 ---"
gmx grompp -f "${MDP_DIR}/md.mdp" \
 -c npt.gro \
 -t npt.cpt \
 -p "${PROTEIN_NAME}.top" \
 -pp processed.top \
 -n index.ndx \
 -o md.tpr \
 -maxwarn 1 \
 || { echo "错误: md grompp 步骤失败"; exit 1; }
echo "--- md.tpr 文件准备完成 ---"

echo "--- 运行正式模拟 (mdrun GPU) ---"
gmx mdrun -v -deffnm md \
 -nt "${NT_THREADS}" \
 -nb gpu \
 -dlb yes \
 || { echo "错误: md mdrun 步骤失败"; exit 1; }
echo "--- 正式模拟完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 10. 清理和结果归档
# -----------------------------------------------------------------------------
echo "--- 步骤 10: 清理和结果归档 ---"
echo "将最终模拟结果移动到 $OUTPUT_DIR..."
[ -f em.gro ] && mv em.gro "${OUTPUT_DIR}/"
[ -f nvt.gro ] && mv nvt.gro "${OUTPUT_DIR}/"
[ -f npt.gro ] && mv npt.gro "${OUTPUT_DIR}/"
[ -f md.gro ] && mv md.gro "${OUTPUT_DIR}/"
[ -f processed.top ] && mv processed.top "${OUTPUT_DIR}/"
[ -f "${PROTEIN_NAME}.top" ] && mv "${PROTEIN_NAME}.top" "${OUTPUT_DIR}/"
[ -f md.log ] && mv md.log "${OUTPUT_DIR}/"
[ -f md.edr ] && mv md.edr "${OUTPUT_DIR}/"
[ -f md.xtc ] && mv md.xtc "${OUTPUT_DIR}/"
[ -f md.cpt ] && mv md.cpt "${OUTPUT_DIR}/"
echo "--- 结果归档完成 ---"
echo ""

# -----------------------------------------------------------------------------
# 11. 报告体系原子数
# -----------------------------------------------------------------------------
echo "--- 步骤 11: 报告体系原子数 ---"
ATOM_COUNT_FILE="${OUTPUT_DIR}/atom_count.out"
TPR_FOR_ATOM_COUNT="${WORKSPACE_DIR}/md.tpr" # 优先使用md.tpr

if [ ! -f "${TPR_FOR_ATOM_COUNT}" ]; then
    # 如果md.tpr不存在，尝试使用em.tpr (可能是流程早期出错)
    TPR_FOR_ATOM_COUNT="${WORKSPACE_DIR}/em.tpr"
fi

if [ -f "${TPR_FOR_ATOM_COUNT}" ]; then
    ATOM_COUNT=$(gmx dump -s "${TPR_FOR_ATOM_COUNT}" 2>/dev/null | grep "System has" | awk '{print $3}')
    if [[ "$ATOM_COUNT" =~ ^[0-9]+$ ]]; then
        echo "$ATOM_COUNT" > "${ATOM_COUNT_FILE}"
        echo "体系原子数 ($ATOM_COUNT) 已写入 ${ATOM_COUNT_FILE} 文件 (使用 ${TPR_FOR_ATOM_COUNT} 计算)。"
    else
        echo "错误: 从 ${TPR_FOR_ATOM_COUNT} 中提取到的原子数格式不正确或未能提取。" > "${ATOM_COUNT_FILE}"
        echo "请检查 ${TPR_FOR_ATOM_COUNT} 文件或 gmx dump 输出。"
    fi
else
    echo "错误: TPR 文件 (尝试了 ${WORKSPACE_DIR}/md.tpr 和 ${WORKSPACE_DIR}/em.tpr) 未在工作目录中找到，无法计算原子数。" > "${ATOM_COUNT_FILE}"
fi
echo "--- 原子数报告完成 ---"
echo ""

# 在原子数统计完成后归档 md.tpr，供轨迹后处理和分析使用。
[ -f "${WORKSPACE_DIR}/md.tpr" ] && mv "${WORKSPACE_DIR}/md.tpr" "${OUTPUT_DIR}/"

echo "--- GROMACS 模拟流程 (GPU版) 全部完成！ ---"
echo "所有中间文件在 $WORKSPACE_DIR 目录下。"
echo "最终结果文件在 $OUTPUT_DIR 目录下。"
echo "结束时间: $(date)"
