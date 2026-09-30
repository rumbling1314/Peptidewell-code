#!/bin/bash
# -----------------------------------------------------------------------------
# AF2 候选 binder 提取与 Rosetta 界面评估
#
# 工作目录应包含 AF2 initial guess 生成的 out.silent。
# 脚本会:
#   1. 根据 counts.py 的统计结果生成 select.txt
#   2. 从 out.silent 中提取指定 PDB
#   3. 使用 Rosetta score_jd2 打分
#   4. 使用 InterfaceAnalyzer 评估界面指标
#
# 运行前请修改下方“用户环境配置”中的路径。
# -----------------------------------------------------------------------------

ulimit -v unlimited || true

# ---------- 用户环境配置 ----------
CONDA_BASE="/home/你的Conda安装目录"
SOFTWARE_ROOT="/home/你的软件目录"

HELPER_SCRIPT_DIR="$SOFTWARE_ROOT/helper_scripts"
SILENT_TOOLS_DIR="$SOFTWARE_ROOT/silent_tools"
ROSETTA_BIN="$SOFTWARE_ROOT/rosetta/main/source/bin"

COUNTS_PY="$HELPER_SCRIPT_DIR/counts.py"
PACK_OPTIONS="$HELPER_SCRIPT_DIR/pack_input_options.txt"
# ----------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

for required_file in out.silent "$COUNTS_PY"; do
    if [[ ! -f "$required_file" ]]; then
        echo "ERROR: missing required file: $required_file" >&2
        exit 1
    fi
done

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate proteinmpnn_binder_design

cp "$COUNTS_PY" ./
python counts.py > select.txt

mkdir -p select
mkdir -p last

while IFS= read -r line
do
    if [[ "$line" =~ [0-9]+_[0-9]+_dldesign_0_cycle1_af2pred ]]
    then
        x="${BASH_REMATCH[0]}"
        "$SILENT_TOOLS_DIR/silentextractspecific" \
            out.silent "$x"
    fi
done < select.txt

shopt -s nullglob
pdb_files=(./*.pdb)

if (( ${#pdb_files[@]} == 0 ))
then
    echo "没有从 out.silent 中提取到 PDB 文件"
    exit 1
fi

mv "${pdb_files[@]}" select/

cd last || exit 1

if [[ ! -f "$PACK_OPTIONS" ]]; then
    echo "ERROR: missing required file: $PACK_OPTIONS" >&2
    exit 1
fi
cp "$PACK_OPTIONS" ./

conda activate af2_binder_design

"$ROSETTA_BIN/score_jd2.linuxgccrelease" \
    -s ../select/*.pdb \
    -no_optH false \
    -ignore_unrecognized_res

"$ROSETTA_BIN/InterfaceAnalyzer.linuxgccrelease" \
    -s ../select/*.pdb \
    @pack_input_options.txt
