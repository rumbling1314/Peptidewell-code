#!/bin/bash
# ==============================================================
# 01_download_proteomes.sh
# 功能：批量下载 32 种人类肠道菌参考基因组的蛋白序列（NCBI RefSeq）
# 依赖：NCBI Datasets CLI（datasets，v14+）
#       安装：https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
# 用法：bash 01_download_proteomes.sh [输出目录]
#       默认输出目录：./proteomes
# 说明：使用 RefSeq 组装编号（GCF_xxx）精确下载目标菌株，
#       每条基因组输出一个 <菌名>.faa（蛋白 FASTA）
# ==============================================================

set -u
OUTDIR="${1:-proteomes}"
mkdir -p "$OUTDIR"

# ---------- 32 菌株：菌名 <TAB> RefSeq 组装编号 ----------
STRAINS="
Faecalibacterium_prausnitzii	GCF_000154385.1
Bacteroides_vulgatus	GCF_000012825.1
Agathobacter_rectalis	GCF_022453685.1
Ruminococcus_bromii	GCF_928721825.1
Blautia_wexlerae	GCF_025148125.1
Fusicatenibacter_saccharivorans	GCF_027667425.1
Anaerostipes_hadrus	GCF_030864025.1
Phascolarctobacterium_faecium	GCF_003945365.1
Roseburia_faecis	GCF_045061065.1
Coprococcus_eutactus	GCF_025149915.1
Butyrivibrio_crossotus	GCF_025148445.1
Akkermansia_muciniphila	GCF_009731575.1
Alistipes_shahii	GCF_025145845.1
Parabacteroides_distasonis	GCF_018279895.1
Bacteroides_ovatus	GCF_001314995.1
Anaerobutyricum_hallii	GCF_018140635.1
Bacteroides_fragilis	GCF_000025985.1
Mitsuokella_jalaludinii	GCF_928722185.1
Bacteroides_thetaiotaomicron	GCF_014131755.1
Ruminococcus_lactaris	GCF_025152405.1
Ruminococcus_callidus	GCF_049562855.1
Roseburia_intestinalis	GCF_900537995.1
Eubacterium_siraeum	GCF_001404415.1
Dysosmobacter_welbionis	GCF_965585535.1
Odoribacter_laneus	GCF_000243215.1
Flavonifractor_plautii	GCF_010508875.1
Eubacterium_ramulus	GCF_001406295.1
Bacillus_subtilis	GCF_000009045.1
Clostridium_scindens	GCF_004295125.1
Lactobacillus_gasseri	GCF_000014425.1
Enterocloster_clostridioformis	GCF_012273195.1
Christensenella_minuta	GCF_003628755.1
"

TOTAL=$(echo "$STRAINS" | grep -c "GCF_")
COUNT=0

echo "=========================================="
echo "下载 $TOTAL 个菌株蛋白组 → $OUTDIR/"
echo "=========================================="

echo "$STRAINS" | while IFS=$'\t' read -r name acc; do
    [ -z "$name" ] && continue
    COUNT=$((COUNT + 1))
    faa="$OUTDIR/${name}.faa"

    if [ -s "$faa" ]; then
        echo "[$COUNT/$TOTAL] 已存在，跳过: $name"
        continue
    fi

    echo "[$COUNT/$TOTAL] 下载 $name ($acc) ..."
    zip="$OUTDIR/${name}.zip"

    datasets download genome accession "$acc" --include protein --filename "$zip" 2>/dev/null

    if [ ! -s "$zip" ]; then
        echo "  ⚠️  下载失败: $name ($acc)"
        continue
    fi

    # 解压并取出 protein.faa
    tmpdir="$OUTDIR/${name}_tmp"
    mkdir -p "$tmpdir"
    unzip -q "$zip" -d "$tmpdir"
    src=$(find "$tmpdir" -name "*protein.faa" | head -1)
    if [ -n "$src" ]; then
        mv "$src" "$faa"
        n=$(grep -c "^>" "$faa")
        echo "  ✓ 完成: $name（$n 条蛋白）"
    else
        echo "  ⚠️  未找到 protein.faa: $name"
    fi
    rm -rf "$tmpdir" "$zip"
done

echo ""
echo "=========================================="
echo "下载结束。蛋白组文件：$OUTDIR/*.faa"
ls -lh "$OUTDIR"/*.faa 2>/dev/null | awk '{print $9, $5}'
