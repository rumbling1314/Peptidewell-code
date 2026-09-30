#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 GROMACS 输出的 .xvg 结果批量转换为 PNG 图片。

用法：
  1. cd 到 xvg 所在目录后运行：
       python /path/to/06_plot_xvg_results.py
  2. 或指定输入和输出目录：
       python 06_plot_xvg_results.py <xvg目录> <图片输出目录>

依赖：Python 3、NumPy、Matplotlib。
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SPECS = {
    "temp.xvg": {
        "kind": "series",
        "time_unit": "ps",
        "ylabel": "Temperature (K)",
        "title": "System Temperature",
        "output": "temperature.png",
        "color": "#dc2626",
        "reference": 300.0,
        "reference_label": "300 K",
    },
    "potential.xvg": {
        "kind": "series",
        "time_unit": "ps",
        "ylabel": "Potential Energy (kJ/mol)",
        "title": "Total Potential Energy",
        "output": "potential.png",
        "color": "#334155",
    },
    "pressure.xvg": {
        "kind": "series",
        "time_unit": "ps",
        "ylabel": "Pressure (bar)",
        "title": "System Pressure",
        "output": "pressure.png",
        "color": "#ea580c",
        "reference": 1.0,
        "reference_label": "1 bar",
    },
    "rmsd.xvg": {
        "kind": "series",
        "time_unit": "ns",
        "ylabel": "Backbone RMSD (nm)",
        "title": "Backbone RMSD",
        "output": "rmsd.png",
        "color": "#2563eb",
    },
    "rg.xvg": {
        "kind": "series",
        "time_unit": "ns",
        "ylabel": "Radius of Gyration (nm)",
        "title": "Radius of Gyration",
        "output": "rg.png",
        "color": "#059669",
    },
    "rmsf.xvg": {
        "kind": "rmsf",
        "xlabel": "Residue Number",
        "ylabel": "RMSF (nm)",
        "title": "Per-Residue RMSF",
        "output": "rmsf.png",
        "color": "#059669",
    },
    "mindist.xvg": {
        "kind": "series",
        "time_unit": "ns",
        "ylabel": "Minimum Distance (nm)",
        "title": "Protein-Binder Minimum Distance",
        "output": "mindist.png",
        "color": "#0891b2",
        "reference": 0.4,
        "reference_label": "0.4 nm cutoff",
    },
    "hbond.xvg": {
        "kind": "series",
        "time_unit": "ns",
        "ylabel": "Number of Hydrogen Bonds",
        "title": "Protein-Binder Hydrogen Bonds",
        "output": "hbond.png",
        "color": "#0891b2",
    },
}


def load_xvg(path):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            try:
                values = [float(item) for item in line.split()]
            except ValueError:
                continue
            if len(values) >= 2:
                rows.append(values[:2])

    if not rows:
        raise ValueError("no numeric data rows found")
    return np.asarray(rows, dtype=float)


def plot_time_series(data, spec, output_dir):
    time = data[:, 0]
    value = data[:, 1]

    # gmx energy 默认输出 ps；结构分析脚本已通过 -tu ns 指定纳秒。
    if spec.get("time_unit") == "ps" and np.nanmax(time) > 200:
        time = time / 1000.0

    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    ax.plot(time, value, color=spec["color"], linewidth=1.2)

    reference = spec.get("reference")
    if reference is not None:
        ax.axhline(
            reference,
            color="#dc2626",
            linestyle="--",
            linewidth=1.0,
            label=spec.get("reference_label"),
        )
        ax.legend(fontsize=10, loc="best")

    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel(spec["ylabel"], fontsize=12)
    ax.set_title(spec["title"], fontsize=13, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--")
    ax.set_xlim(np.nanmin(time), np.nanmax(time))

    fig.tight_layout()
    output_path = output_dir / spec["output"]
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_rmsf(data, spec, output_dir):
    residue = data[:, 0]
    rmsf = data[:, 1]

    fig, ax = plt.subplots(figsize=(11, 4), dpi=300)
    ax.plot(residue, rmsf, color=spec["color"], linewidth=1.0)
    ax.set_xlabel(spec["xlabel"], fontsize=12)
    ax.set_ylabel(spec["ylabel"], fontsize=12)
    ax.set_title(spec["title"], fontsize=13, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--")
    ax.set_xlim(np.nanmin(residue), np.nanmax(residue))

    fig.tight_layout()
    output_path = output_dir / spec["output"]
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert GROMACS .xvg result files to PNG figures."
    )
    parser.add_argument(
        "xvg_dir",
        nargs="?",
        default=".",
        help="directory containing .xvg files (default: current directory)",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help="directory for PNG files (default: <xvg_dir>/plots)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    xvg_dir = Path(args.xvg_dir).expanduser().resolve()
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else xvg_dir / "plots"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not xvg_dir.is_dir():
        print(f"ERROR: input directory does not exist: {xvg_dir}", file=sys.stderr)
        return 1

    generated = []
    missing = []

    for name, spec in SPECS.items():
        xvg_path = xvg_dir / name
        if not xvg_path.is_file():
            missing.append(name)
            continue

        try:
            data = load_xvg(xvg_path)
            if spec["kind"] == "rmsf":
                output_path = plot_rmsf(data, spec, output_dir)
            else:
                output_path = plot_time_series(data, spec, output_dir)
        except Exception as exc:
            print(f"SKIP {name}: {exc}", file=sys.stderr)
            continue

        generated.append(output_path)
        print(f"Generated: {output_path}")

    if missing:
        print("Missing expected files: " + ", ".join(missing), file=sys.stderr)

    if not generated:
        print("No figures generated.", file=sys.stderr)
        return 1

    print(f"Completed. {len(generated)} figures saved to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
