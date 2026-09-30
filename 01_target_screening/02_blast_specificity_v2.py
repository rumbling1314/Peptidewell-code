#!/usr/bin/env python3
"""
菌种特有蛋白批量筛选脚本（支持断点续跑）
用法: python3 blast_specificity_v2.py <蛋白组目录> <输出目录>
已有结果文件的菌种自动跳过。
"""

import os
import sys
import subprocess
from pathlib import Path


def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"命令失败: {cmd}")
        print(result.stderr)
        sys.exit(1)
    return result


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    proteome_dir = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    faa_files = sorted(proteome_dir.glob("*.faa"))
    strains = [f.stem for f in faa_files]
    print(f"共 {len(strains)} 个蛋白组\n")

    done = skip = 0

    for i, target in enumerate(strains):
        out_file = output_dir / f"{target}_specificity_ranked.tsv"

        # 断点续跑：已有结果则跳过
        if out_file.exists() and out_file.stat().st_size > 0:
            skip += 1
            print(f"[{i+1}/{len(strains)}] 跳过: {target}")
            continue

        print(f"\n[{i+1}/{len(strains)}] 处理: {target}")

        target_faa = proteome_dir / f"{target}.faa"
        others = [s for s in strains if s != target]
        combined_faa = output_dir / f"_tmp_db_{target}.faa"

        # 合并排除数据库
        with open(combined_faa, 'w') as outf:
            for s in others:
                with open(proteome_dir / f"{s}.faa") as inf:
                    outf.write(inf.read())

        db_prefix = output_dir / f"_tmp_db_{target}"
        run_cmd(f"makeblastdb -in {combined_faa} -dbtype prot -out {db_prefix} -logfile /dev/null")

        blast_out = output_dir / f"_tmp_blast_{target}.tsv"
        cmd = (f"blastp -query {target_faa} -db {db_prefix} "
               f"-evalue 1e-3 -outfmt 6 -max_target_seqs 5 "
               f"-num_threads 4 -out {blast_out}")
        run_cmd(cmd)

        # 解析最佳匹配
        best_hits = {}
        with open(blast_out) as f:
            for line in f:
                p = line.strip().split('\t')
                q, s, pid = p[0], p[1], float(p[2])
                aln_len = int(p[3])
                ev = float(p[10])
                bs = float(p[11])
                if q not in best_hits or bs > best_hits[q][2]:
                    best_hits[q] = (pid, ev, bs, s, aln_len)

        # 所有蛋白ID
        all_prots = []
        with open(target_faa) as f:
            for line in f:
                if line.startswith('>'):
                    all_prots.append(line[1:].split()[0])

        # 排序输出
        ranked = []
        for pid in all_prots:
            if pid in best_hits:
                ident, ev, bs, sid, alen = best_hits[pid]
                strain = "unknown"
                for s in others:
                    if s.replace('_', '') in sid.replace('_', '').lower():
                        strain = s
                        break
                ranked.append((pid, ident, ev, bs, strain, sid, alen))
            else:
                ranked.append((pid, 0.0, 1.0, 0.0, "NONE", "NO_HIT", 0))

        ranked.sort(key=lambda x: (x[1] > 0, x[1], -x[2]))

        with open(out_file, 'w') as f:
            f.write("protein_id\tidentity_%\tevalue\tbitscore\tmatched_strain\tsubject_id\taln_length\n")
            for row in ranked:
                f.write("\t".join(str(v) for v in row) + "\n")

        n_nohit = sum(1 for r in ranked if r[4] == "NONE")
        n_low = sum(1 for r in ranked if 0 < r[1] < 25)
        print(f"  完成: 总{len(ranked)} | 无同源{n_nohit} | <25%{n_low}")

        # 清理临时文件
        combined_faa.unlink(missing_ok=True)
        blast_out.unlink(missing_ok=True)
        for ext in ['.phr', '.pin', '.psq']:
            (output_dir / f"_tmp_db_{target}{ext}").unlink(missing_ok=True)

        done += 1

    print(f"\n完成！新增 {done} 个，跳过 {skip} 个")
    print(f"结果目录: {output_dir}/")


if __name__ == '__main__':
    main()
