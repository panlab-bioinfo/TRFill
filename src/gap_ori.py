#!/usr/bin/env python
import sys
from collections import defaultdict

def parse_config(config_file):
    """解析配置文件，返回gap列表"""
    config_vars = {}
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip('();')
                if key in ['chrs', 'gap_starts', 'gap_ends']:
                    items = [item.strip(' "\'') for item in value.split(',')]
                    if key in ['gap_starts', 'gap_ends']:
                        items = list(map(int, items))
                    config_vars[key] = items
    gaps = list(zip(config_vars['chrs'], config_vars['gap_starts'], config_vars['gap_ends']))
    return gaps

def process_paf(paf_file):
    """读取并过滤PAF文件，保留高质量主要比对"""
    paf_records = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:  # 确保包含基本字段
                continue

            # 解析基础字段
            try:
                mapq = int(parts[11])  # 第12列是MAPQ
            except ValueError:
                continue
            
            # 质量控制过滤
            if mapq <= 30:
                continue

            # 解析比对类型标签
            is_primary = False
            for opt in parts[12:]:
                if opt.startswith('tp:A:'):
                    if opt == 'tp:A:P':
                        is_primary = True
                    break
            
            if not is_primary:
                continue
            if parts[0] == parts[5]:
                continue
            # 存储有效记录
            record = {
                'query_length': int(parts[1]),
                'query_start': int(parts[2]),
                'query_end': int(parts[3]),
                'target_name': parts[5],
                'target_start': int(parts[7]),
                'target_end': int(parts[8]),
            }
            paf_records[parts[0]].append(record)
    return paf_records

def determine_orientation(gap, paf_records):
    """确定单个gap的填回方向"""
    chr_name, gap_start, gap_end = gap
    gap_id = f"{chr_name}_gap"
    records = paf_records.get(gap_id, [])
    if not records:
        return None

    # 获取gap序列长度
    L = records[0]['query_length']
    # 确定左右部分
    if L < 200000:
        left_part_end = L // 2
        right_part_start = L // 2
    else:
        left_part_end = 100000
        right_part_start = L - 100000

    # 岸区范围
    left_flank_start = max(0, gap_start - 20000)
    left_flank_end = gap_start
    right_flank_start = gap_end
    right_flank_end = gap_end + 20000

    counters = {'left_left':0, 'left_right':0, 'right_left':0, 'right_right':0}

    for record in records:
        if record['target_name'] != chr_name:
            continue

        # 判断query部分
        mid = (record['query_start'] + record['query_end']) // 2
        if mid < left_part_end:
            part = 'left'
        elif mid >= right_part_start:
            part = 'right'
        else:
            continue

        # 判断target区域
        ts, te = record['target_start'], record['target_end']
        overlap_left = (ts < left_flank_end) and (te > left_flank_start)
        overlap_right = (ts < right_flank_end) and (te > right_flank_start)
        if overlap_left and overlap_right:
            continue  # 忽略重叠两岸的情况
        elif overlap_left:
            flank = 'left'
        elif overlap_right:
            flank = 'right'
        else:
            continue

        # 统计
        key = f"{part}_{flank}"
        if key in counters:
            counters[key] += 1

    correct = counters['left_left'] + counters['right_right']
    incorrect = counters['left_right'] + counters['right_left']
    return '+' if correct > incorrect else '-'

def main(config_path, paf_path):
    gaps = parse_config(config_path)
    paf_data = process_paf(paf_path)
    for gap in gaps:
        orientation = determine_orientation(gap, paf_data)
        print(f"{gap}\t{orientation}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py config.txt input.paf")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])