#!/usr/bin/env python3
"""
通用版：批量提交 NCBI BLASTp 任务
支持任意菌种，只需修改输入 FASTA 文件和输出目录
支持断点续跑

用法：
  python3 batch_blast_ncbi.py \
    --input /path/to/proteins.fasta \
    --output /path/to/output_dir \
    --database nr \
    --evalue 0.001 \
    --hitlist 500

特点：
- 逐条提交，自动保存 RID，支持断点续跑
- 轮询等待结果，完成后自动下载 XML
- 每个蛋白一个 XML 文件
- 失败自动重试
"""

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request


def parse_fasta(fasta_file):
    """简单的FASTA解析器，不依赖BioPython"""
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # 保存上一条序列
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                # 开始新序列
                # 提取ID（第一个空格前的内容，去掉>）
                header = line[1:]
                current_id = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    
    # 保存最后一条序列
    if current_id is not None:
        sequences.append((current_id, ''.join(current_seq)))
    
    return sequences


# ============ NCBI BLAST API 配置 ============
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
POLL_INTERVAL = 30  # 轮询间隔（秒）
MAX_RETRIES = 3     # 最大重试次数
SUBMIT_DELAY = 3    # 每次提交间隔（秒），避免触发NCBI限速
# ============================================


def parse_args():
    parser = argparse.ArgumentParser(description='批量提交 NCBI BLASTp 任务')
    parser.add_argument('--input', '-i', required=True, help='输入 FASTA 文件')
    parser.add_argument('--output', '-o', required=True, help='输出目录（存放XML结果）')
    parser.add_argument('--database', '-db', default='nr', help='BLAST数据库，默认nr')
    parser.add_argument('--program', '-p', default='blastp', help='BLAST程序，默认blastp')
    parser.add_argument('--evalue', '-e', type=float, default=0.001, help='E-value阈值，默认0.001')
    parser.add_argument('--hitlist', '-n', type=int, default=500, help='保留的hit数量，默认500')
    parser.add_argument('--matrix', default='BLOSUM62', help='打分矩阵，默认BLOSUM62')
    parser.add_argument('--gapopen', type=int, default=11, help='Gap打开罚分，默认11')
    parser.add_argument('--gapextend', type=int, default=1, help='Gap延伸罚分，默认1')
    parser.add_argument('--filter', default='F', help='低复杂度过滤，默认F(关闭)')
    return parser.parse_args()


def submit_blast(sequence, program, database, evalue, hitlist, matrix, gapopen, gapextend, filter_str):
    """提交单个BLAST任务，返回RID"""
    params = {
        'CMD': 'Put',
        'PROGRAM': program,
        'DATABASE': database,
        'QUERY': sequence,
        'EXPECT': str(evalue),
        'HITLIST_SIZE': str(hitlist),
        'MATRIX_NAME': matrix,
        'GAPCOSTS': f'{gapopen} {gapextend}',
        'FILTER': filter_str,
        'FORMAT_TYPE': 'XML',
    }
    # key 来源优先级：环境变量 > ~/.ncbi_api_key 文件
    api_key = os.environ.get('NCBI_API_KEY', '')
    if not api_key:
        _kf = os.path.expanduser('~/.ncbi_api_key')
        if os.path.exists(_kf):
            api_key = open(_kf).read().strip()
    if api_key:
        params['API_KEY'] = api_key
    
    data = urllib.parse.urlencode(params).encode('utf-8')
    
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(BLAST_URL, data=data)
            with urllib.request.urlopen(req, timeout=30) as response:
                result = response.read().decode('utf-8')
            
            # 提取 RID
            rid = None
            for line in result.split('\n'):
                if 'RID =' in line:
                    rid = line.split('RID =')[1].strip()
                    break
            
            if rid:
                return rid
            else:
                print(f'    警告：未获取到RID，重试 {attempt+1}/{MAX_RETRIES}')
                time.sleep(5)
                
        except Exception as e:
            print(f'    提交失败: {e}，重试 {attempt+1}/{MAX_RETRIES}')
            time.sleep(5)
    
    return None


def check_blast_status(rid):
    """检查BLAST任务状态"""
    params = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_OBJECT': 'SearchInfo',
    }
    
    url = f"{BLAST_URL}?{urllib.parse.urlencode(params)}"
    
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            result = response.read().decode('utf-8')
        
        # 提取状态
        status = None
        for line in result.split('\n'):
            if 'Status=' in line:
                status = line.split('Status=')[1].strip()
                break
        
        return status
        
    except Exception as e:
        print(f'    状态检查失败: {e}')
        return None


def get_blast_results(rid):
    """获取BLAST结果（XML格式）"""
    params = {
        'CMD': 'Get',
        'RID': rid,
        'FORMAT_TYPE': 'XML',
    }
    
    url = f"{BLAST_URL}?{urllib.parse.urlencode(params)}"
    
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                result = response.read().decode('utf-8')
            
            if 'BlastOutput' in result:
                return result
            else:
                print(f'    结果获取不完整，重试 {attempt+1}/{MAX_RETRIES}')
                time.sleep(5)
                
        except Exception as e:
            print(f'    结果下载失败: {e}，重试 {attempt+1}/{MAX_RETRIES}')
            time.sleep(5)
    
    return None


def load_rids(rids_file):
    """加载已提交的RID记录（断点续跑用）"""
    rids = {}
    if os.path.exists(rids_file):
        with open(rids_file, 'r') as f:
            for line in f:
                if line.strip() and '\t' in line:
                    protein_id, rid, status = line.strip().split('\t', 2)
                    rids[protein_id] = {'rid': rid, 'status': status}
    return rids


def save_rids(rids, rids_file):
    """保存RID记录"""
    with open(rids_file, 'w') as f:
        for protein_id, info in rids.items():
            f.write(f"{protein_id}\t{info['rid']}\t{info['status']}\n")


def main():
    args = parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    rids_file = os.path.join(args.output, 'rids.tsv')
    
    print(f'输入文件: {args.input}')
    print(f'输出目录: {args.output}')
    print(f'数据库: {args.database}')
    print(f'程序: {args.program}')
    print(f'E-value: {args.evalue}')
    print(f'Hitlist大小: {args.hitlist}')
    print()
    
    # 读取FASTA文件
    sequences = parse_fasta(args.input)
    
    print(f'共读取 {len(sequences)} 条蛋白序列')
    
    # 加载已提交的RID（断点续跑）
    rids = load_rids(rids_file)
    print(f'已提交记录: {len(rids)} 个')
    
    # 统计已完成的
    completed = 0
    for protein_id, info in rids.items():
        xml_file = os.path.join(args.output, f'{protein_id}.xml')
        if os.path.exists(xml_file):
            completed += 1
    
    print(f'已完成: {completed} 个')
    print()
    
    # ===== 第一阶段：提交所有任务 =====
    print('===== 阶段1：提交BLAST任务 =====')
    
    submitted_count = 0
    for i, (protein_id, seq) in enumerate(sequences):
        # 跳过已提交的
        if protein_id in rids:
            continue
        
        print(f'[{i+1}/{len(sequences)}] 提交 {protein_id}...', end=' ', flush=True)
        
        rid = submit_blast(
            seq, args.program, args.database, args.evalue, args.hitlist,
            args.matrix, args.gapopen, args.gapextend, args.filter
        )
        
        if rid:
            rids[protein_id] = {'rid': rid, 'status': 'submitted'}
            save_rids(rids, rids_file)
            submitted_count += 1
            print(f'RID={rid} ✓')
        else:
            print('✗ 提交失败')
            rids[protein_id] = {'rid': 'FAILED', 'status': 'failed'}
            save_rids(rids, rids_file)
        
        # 限速
        time.sleep(SUBMIT_DELAY)
    
    print(f'\n提交完成：新提交 {submitted_count} 个，总计 {len(rids)} 个')
    
    # ===== 第二阶段：等待并下载结果 =====
    print('\n===== 阶段2：等待结果并下载 =====')
    
    while True:
        pending = []
        done_count = 0
        failed_count = 0
        
        for protein_id, info in rids.items():
            xml_file = os.path.join(args.output, f'{protein_id}.xml')
            
            if os.path.exists(xml_file):
                done_count += 1
                continue
            
            if info['status'] == 'failed':
                failed_count += 1
                continue
            
            if info['rid'] == 'FAILED':
                failed_count += 1
                continue
            
            pending.append(protein_id)
        
        total = len(rids)
        print(f'\n进度: {done_count}/{total} 完成, {len(pending)} 等待中, {failed_count} 失败')
        
        if not pending:
            print('\n所有任务完成！')
            break
        
        # 检查每个pending任务的状态
        new_completed = 0
        
        for protein_id in pending:
            info = rids[protein_id]
            rid = info['rid']
            
            print(f'  检查 {protein_id} (RID: {rid})...', end=' ', flush=True)
            
            status = check_blast_status(rid)
            
            if status == 'READY':
                print('完成！下载中...', end=' ', flush=True)
                result = get_blast_results(rid)
                
                if result:
                    xml_file = os.path.join(args.output, f'{protein_id}.xml')
                    with open(xml_file, 'w') as f:
                        f.write(result)
                    
                    info['status'] = 'completed'
                    save_rids(rids, rids_file)
                    new_completed += 1
                    print('✓')
                else:
                    print('✗ 下载失败')
                    
            elif status == 'WAITING' or status == 'RUNNING':
                print(f'{status}')
            else:
                print(f'状态未知: {status}')
            
            # 限速
            time.sleep(1)
        
        if new_completed == 0:
            print(f'\n没有新完成的任务，等待 {POLL_INTERVAL} 秒后重试...')
            time.sleep(POLL_INTERVAL)
    
    # ===== 最终统计 =====
    print('\n===== 最终统计 =====')
    success = 0
    failed = 0
    
    for protein_id, info in rids.items():
        xml_file = os.path.join(args.output, f'{protein_id}.xml')
        if os.path.exists(xml_file):
            success += 1
        else:
            failed += 1
            print(f'  失败: {protein_id}')
    
    print(f'\n成功: {success}')
    print(f'失败: {failed}')
    print(f'结果保存在: {args.output}')


if __name__ == '__main__':
    main()
