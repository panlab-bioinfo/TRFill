#!/usr/bin/env python3
import sys
import argparse

def read_and_compare_values(log_file):
    with open(log_file, 'r') as file:
        lines = file.readlines()
        if len(lines) < 2:
            print("Error: Log file does not contain enough lines.")
            return
        

        last_line = lines[-1].strip()
        second_last_line = lines[-2].strip()
        
        try:
            val1 = float(last_line.split()[1])
            val2 = float(second_last_line.split()[1])
        except (IndexError, ValueError) as e:
            print(f"Error processing lines: {e}")
            return
        
        # 比较两个浮点数值
        print(val1, val2)
        if val1 >= val2:
            print("0")
        elif val1 < val2:
            print("1")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare the last two lines' second column values in a log file.")
    parser.add_argument('log_file', type=str, help='Path to the log file')
    
    args = parser.parse_args()
    read_and_compare_values(args.log_file)