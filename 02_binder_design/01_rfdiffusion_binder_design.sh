#!/bin/bash
# -----------------------------------------------------------------------------
# RFdiffusion + ProteinMPNN + AlphaFold2 binder 设计流程
#
# 输入:
#   - 目标蛋白 PDB: ${WORK_ROOT}/binder260906/001/input.pdb
#
# 输出:
#   - project_output/: RFdiffusion 生成的 PDB
#   - project_output/test/design.silent: ProteinMPNN 输入 silent 文件
#   - project_output/test/design2.silent: ProteinMPNN 序列设计结果
#   - AlphaFold2 initial guess 预测结果
#
# 运行前请修改下方“用户环境配置”中的路径。
# -----------------------------------------------------------------------------

# Avoid locale warnings from the login environment.
export LANG=C
export LC_ALL=C
export PYTHONUNBUFFERED=1

# This removes only the shell's virtual-memory limit if the system permits it.
ulimit -v unlimited || true

# ---------- 用户环境配置 ----------
WORK_ROOT="/home/你的工作目录"
CONDA_BASE="/home/你的Conda安装目录"
SOFTWARE_ROOT="/home/你的软件目录"

# 任务配置
this_dir="$WORK_ROOT/binder260906/001"  # 001 为示例编号，具体以实际项目为准
name_protein=001
output_dir="$this_dir/project_output"
test_dir="$output_dir/test"

# 软件路径
RFDIFFUSION_DIR="$SOFTWARE_ROOT/RFdiffusion"
DL_BINDER_DESIGN_DIR="$SOFTWARE_ROOT/dl_binder_design"
SILENT_TOOLS_DIR="$SOFTWARE_ROOT/silent_tools"
MPNN_WEIGHTS="$DL_BINDER_DESIGN_DIR/mpnn_fr/ProteinMPNN/soluble_model_weights/v_48_020.pt"
# ----------------------------------

source "$CONDA_BASE/etc/profile.d/conda.sh"

activate_env() {
    conda activate "$1"
    echo "[INFO] environment: $(basename "$CONDA_PREFIX")"
    echo "[INFO] python: $(which python)"
}

check_torch_gpu() {
    python - <<'PY'
import torch

if not torch.cuda.is_available():
    raise SystemExit("ERROR: PyTorch cannot see a CUDA GPU")
if torch.cuda.device_count() < 1:
    raise SystemExit("ERROR: PyTorch sees zero CUDA devices")
print(f"[INFO] torch={torch.__version__}")
print(f"[INFO] cuda_available={torch.cuda.is_available()}")
print(f"[INFO] cuda_device_count={torch.cuda.device_count()}")
print(f"[INFO] cuda_device_0={torch.cuda.get_device_name(0)}")
PY
}

if [[ ! -f "$this_dir/input.pdb" ]]; then
    echo "ERROR: missing input file: $this_dir/input.pdb" >&2
    exit 1
fi

# Do not silently mix a rerun with old RFdiffusion results.
if compgen -G "$output_dir/${name_protein}_*.pdb" > /dev/null; then
    echo "ERROR: existing design PDB files were found in $output_dir" >&2
    echo "Please move the old project_output directory aside before rerunning." >&2
    exit 1
fi

mkdir -p "$output_dir" "$test_dir"
cd "$this_dir"

# RFdiffusion secondary-structure and adjacency features
activate_env SE3nv
check_torch_gpu

python "$RFDIFFUSION_DIR/helper_scripts/make_secstruc_adj.py" \
    --input_pdb "$this_dir/input.pdb" \
    --out_dir "$this_dir/target_secstruct"

# RFdiffusion
python "$RFDIFFUSION_DIR/scripts/run_inference.py" \
    'ppi.hotspot_res=[A182,A155,A169,A175,A189,A204]' \
    scaffoldguided.target_path="$this_dir/input.pdb" \
    scaffoldguided.target_ss="$this_dir/target_secstruct/input_ss.pt" \
    scaffoldguided.target_adj="$this_dir/target_secstruct/input_adj.pt" \
    inference.output_prefix="$output_dir/$name_protein" \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.target_pdb=True \
    scaffoldguided.scaffold_dir="$RFDIFFUSION_DIR/examples/ppi_scaffolds" \
    inference.num_designs=2000 \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0

# 备注：上方 ppi.hotspot_res 为实例热点，具体以实际为准。

shopt -s nullglob
design_pdbs=("$output_dir/${name_protein}_"*.pdb)
if [[ ${#design_pdbs[@]} -ne 2000 ]]; then
    echo "ERROR: RFdiffusion produced ${#design_pdbs[@]} design PDB files; expected 2000." >&2
    printf '%s\n' "${design_pdbs[@]}" >&2
    exit 1
fi
echo "[INFO] RFdiffusion produced ${#design_pdbs[@]} design PDB files"

# ProteinMPNN
activate_env proteinmpnn_binder_design
check_torch_gpu

python "$DL_BINDER_DESIGN_DIR/helper_scripts/addFIXEDlabels.py" \
    --pdbdir "$output_dir" \
    --trbdir "$output_dir" \
    --verbose

"$SILENT_TOOLS_DIR/silentfrompdbs" \
    "${design_pdbs[@]}" \
    > "$test_dir/design.silent"

cd "$test_dir"
rm -f check.point out.silent design2.silent design2.silent.idx

python "$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py" \
    -silent design.silent \
    -checkpoint_path "$MPNN_WEIGHTS"

if [[ ! -s out.silent ]]; then
    echo "ERROR: ProteinMPNN did not create a non-empty out.silent" >&2
    exit 1
fi
mv out.silent design2.silent

# AF2 initial-guess prediction
activate_env af2_binder_design

python - <<'PY'
import jax

devices = jax.devices()
print(f"[INFO] JAX devices: {devices}")
if not any(device.platform == "gpu" for device in devices):
    raise SystemExit("ERROR: JAX cannot see a GPU")
PY

python "$DL_BINDER_DESIGN_DIR/af2_initial_guess/predict.py" \
    -silent "$test_dir/design2.silent"

echo "[INFO] RFdiffusion, ProteinMPNN and AF2 stages completed."
