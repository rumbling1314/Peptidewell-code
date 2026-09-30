#!/bin/bash
# ============================================================
# 通用 GROMACS 轨迹分析脚本 (自动解析链分组)
# 用法:
#   1. 修改下面“用户环境配置”区的参数
#   2. 在配置完成后运行: bash 05_gromacs_md_analysis.sh
# 产出: output/xvg/ 下 8 个 xvg 文件
# ============================================================
set -euo pipefail

# ---------- 用户环境配置 ----------
WORK_ROOT="/home/你的工作目录"
GMX_RC="/opt/gromacsgpu/bin/GMXRC"

PROTEIN_NAME="001"                                                                      # 001 为示例编号，具体以实际项目为准
WORK_DIR="$WORK_ROOT/binder260906/MD/001/output"                                        # 本任务 output 目录
PDB_INPUT="$WORK_ROOT/binder260906/MD/001/input/001_624_dldesign_0_cycle1_af2pred.pdb"  # 原始 PDB(自动解析链范围用)
BINDER_CHAIN="A"                                                                        # binder 多肽所在链 ID (A 或 B)
# ------------------------------------

OUT_DIR="${WORK_DIR}/xvg"

cd "$WORK_DIR" || { echo "错误: 无法进入 $WORK_DIR"; exit 1; }
mkdir -p "$OUT_DIR"

set +u
source "$GMX_RC"
set -u

# 检查必需输入文件
for f in md.tpr md_fit.xtc md.edr md_fit.gro; do
    [[ -f "$f" ]] || { echo "错误: 缺少必需文件 $f"; exit 1; }
done

echo "==== [0/4] 从 PDB 自动解析链残基范围 ===="
if [[ ! -f "$PDB_INPUT" ]]; then
    echo "错误: PDB 文件不存在: $PDB_INPUT"
    exit 1
fi
# 解析两条链的残基范围(写入临时文件, 便于查看)
grep "^ATOM" "$PDB_INPUT" | awk -v bc="$BINDER_CHAIN" '{
    chain = substr($0,22,1)
    res   = substr($0,23,4) + 0
    if (!(chain in min) || res < min[chain]) min[chain] = res
    if (!(chain in max) || res > max[chain]) max[chain] = res
} END {
    if (!(bc in min)) { print "0 0 0 0"; exit }
    oc = ""
    for (c in min) if (c != bc) { oc = c; break }
    if (oc == "") { print "0 0 0 0"; exit }
    printf "%d %d %d %d", min[bc], max[bc], min[oc], max[oc]
}' > chain_ranges.txt || { echo "PDB 解析失败"; exit 1; }

echo "解析结果: $(cat chain_ranges.txt)"
read PEP_MIN PEP_MAX PROT_MIN PROT_MAX < chain_ranges.txt || true
if [[ -z "${PEP_MIN:-}" || "${PEP_MIN:-0}" == "0" || -z "${PROT_MIN:-}" || "${PROT_MIN:-0}" == "0" ]]; then
    echo "错误: 未能从 PDB 解析两条链, 请检查 PDB_INPUT 路径和 BINDER_CHAIN 设置"
    echo "PDB 文件: $PDB_INPUT"
    exit 1
fi
echo "binder(链$BINDER_CHAIN): 残基 $PEP_MIN-$PEP_MAX | protein: 残基 $PROT_MIN-$PROT_MAX"

# ============ 1. 热力学数据 ============
echo "==== [1/4] 热力学数据 ===="
echo Temperature | gmx energy -f md.edr -o "$OUT_DIR/temp.xvg" -xvg none
echo Potential   | gmx energy -f md.edr -o "$OUT_DIR/potential.xvg" -xvg none
echo Pressure    | gmx energy -f md.edr -o "$OUT_DIR/pressure.xvg" -xvg none

# ============ 2. 结构稳定性 ============
echo "==== [2/4] 结构稳定性 ===="
printf "4\n4\n" | gmx rms   -s md.tpr -f md_fit.xtc -o "$OUT_DIR/rmsd.xvg" -tu ns -xvg none
printf "1\n"   | gmx gyrate -s md.tpr -f md_fit.xtc -o "$OUT_DIR/rg.xvg" -tu ns -xvg none
printf "4\n"   | gmx rmsf   -s md.tpr -f md_fit.xtc -o "$OUT_DIR/rmsf.xvg" -res -xvg none

# ============ 3. 蛋白-多肽结合分析 ============
echo "==== [3/4] 蛋白-多肽结合分析 ===="
rm -f index_chain.ndx
printf "r ${PEP_MIN}-${PEP_MAX}\nr ${PROT_MIN}-${PROT_MAX}\nq\n" | gmx make_ndx -f md.tpr -o index_chain.ndx

get_group_id() {
    awk -v tgt="$1" '
        /^\[/ {
            id++
            name = $0
            sub(/^\[[[:space:]]*/, "", name)
            sub(/[[:space:]]*\]$/, "", name)
        }
        name == tgt { print id - 1; exit }
    ' "$2"
}
PEP_GROUP="r_${PEP_MIN}-${PEP_MAX}"
PROT_GROUP="r_${PROT_MIN}-${PROT_MAX}"
PEP_ID=$(get_group_id "$PEP_GROUP" index_chain.ndx)
PROT_ID=$(get_group_id "$PROT_GROUP" index_chain.ndx)
if [[ -z "$PEP_ID" || -z "$PROT_ID" ]]; then
    echo "错误: 索引中找不到 $PEP_GROUP / $PROT_GROUP 组"
    grep '^\[' index_chain.ndx
    exit 1
fi
echo "组编号: $PEP_GROUP=$PEP_ID  $PROT_GROUP=$PROT_ID"

# 蛋白-多肽最小距离 (0.4 nm = 4Å 接触截断)
printf "${PEP_ID}\n${PROT_ID}\n" | gmx mindist -s md.tpr -f md_fit.xtc -n index_chain.ndx \
    -od "$OUT_DIR/mindist.xvg" -d 0.4 -tu ns -xvg none

# 蛋白-多肽氢键数目
printf "${PEP_ID}\n${PROT_ID}\n" | gmx hbond -s md.tpr -f md_fit.xtc -n index_chain.ndx \
    -num "$OUT_DIR/hbond.xvg" -tu ns -xvg none

# ============ 4. 汇总 ============
echo "==== [4/4] 全部完成 ===="
ls -lh "$OUT_DIR"/*.xvg
echo "分析完成!"
