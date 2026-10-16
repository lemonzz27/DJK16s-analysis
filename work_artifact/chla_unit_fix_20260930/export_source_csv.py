import os
from pathlib import Path
import csv
import openpyxl
root = Path(__file__).resolve().parents[2]
src = root / '01_分析项目' / '环境变化曲线' / '环境因子_按月份平均_17个月.xlsx'
out = Path(__file__).resolve().parent / 'monthly_environment_source.csv'
wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
ws = wb['按月平均']
with out.open('w', newline='', encoding='utf-8-sig') as handle:
    csv.writer(handle).writerows(ws.values)
print(out)
