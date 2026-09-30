#!/usr/bin/env python3
"""
按同源性分级筛选特有蛋白
S级: 完全无同源 (matched_strain == NONE)
A级: 一致性 < 30%
B级: 一致性 30-40%
C级: 一致性 > 40%
"""

import os
import sys
import csv

def grade_protein(identity_str, matched_strain):
    if matched_strain.strip() == 'NONE':
        return 'S'
    try:
        identity = float(identity_str)
    except:
        return 'C'
    if identity < 30:
        return 'A'
    elif identity < 40:
        return 'B'
    else:
        return 'C'

def load_fasta(fasta_path):
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq)
    return seqs

def process_strain(blast_file, proteome_file, output_dir, grade_output_dir):
    strain_name = os.path.basename(blast_file).replace('_specificity_ranked.tsv', '')
    
    proteins = []
    with open(blast_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pid = row['protein_id']
            identity = row.get('identity_%', '100')
            matched = row.get('matched_strain', 'UNKNOWN')
            grade = grade_protein(identity, matched)
            proteins.append({
                'protein_id': pid,
                'identity_%': identity,
                'evalue': row.get('evalue', ''),
                'bitscore': row.get('bitscore', ''),
                'matched_strain': matched,
                'grade': grade
            })
    
    os.makedirs(grade_output_dir, exist_ok=True)
    out_file = os.path.join(grade_output_dir, f'{strain_name}_graded.tsv')
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['grade', 'protein_id', 'identity_%', 'evalue', 'bitscore', 'matched_strain'], delimiter='\t')
        writer.writeheader()
        for p in proteins:
            writer.writerow(p)
    
    sa_ids = set(p['protein_id'] for p in proteins if p['grade'] in ('S', 'A'))
    seqs = load_fasta(proteome_file)
    
    sa_fasta = os.path.join(output_dir, f'{strain_name}_SA_grade.faa')
    count = 0
    with open(sa_fasta, 'w') as f:
        for pid, seq in seqs.items():
            if pid in sa_ids:
                f.write(f'>{pid}\n{seq}\n')
                count += 1
    
    s_count = sum(1 for p in proteins if p['grade'] == 'S')
    a_count = sum(1 for p in proteins if p['grade'] == 'A')
    b_count = sum(1 for p in proteins if p['grade'] == 'B')
    c_count = sum(1 for p in proteins if p['grade'] == 'C')
    
    print(f'{strain_name}: 总{len(proteins)}个 | S:{s_count} | A:{a_count} | B:{b_count} | C:{c_count} | S+A提取:{count}个')
    return count

def main():
    if len(sys.argv) < 4:
        print('用法: python3 rank_specificity_graded.py <blast_dir> <proteome_dir> <output_dir>')
        sys.exit(1)
    
    blast_dir = sys.argv[1]
    proteome_dir = sys.argv[2]
    output_dir = sys.argv[3]
    grade_output_dir = 'graded_specificity'
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(grade_output_dir, exist_ok=True)
    
    total = 0
    for blast_file in sorted(os.listdir(blast_dir)):
        if not blast_file.endswith('_specificity_ranked.tsv'):
            continue
        strain_name = blast_file.replace('_specificity_ranked.tsv', '')
        proteome_file = os.path.join(proteome_dir, f'{strain_name}.faa')
        
        if not os.path.exists(proteome_file):
            print(f'跳过 {strain_name}: 找不到蛋白组文件')
            continue
        
        count = process_strain(
            os.path.join(blast_dir, blast_file),
            proteome_file,
            output_dir,
            grade_output_dir
        )
        total += count
    
    print(f'\n完成！S+A级序列输出到 {output_dir}/ 目录')
    print(f'分级详情输出到 {grade_output_dir}/ 目录')
    print(f'总计提取 {total} 个 S+A 级蛋白')

if __name__ == '__main__':
    main()
