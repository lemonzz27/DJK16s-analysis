import os
from pathlib import Path
import csv
import openpyxl
root = Path(os.environ.get('DJK_PROJECT_ROOT', '.'))
src = root / '01_分析项目' / '环境变化曲线' / '环境因子_按月份平均_17个月.xlsx'
out = root / 'work_artifact/chla_unit_fix_20260930/monthly_environment_source.csv'
out.parent.mkdir(parents=True, exist_ok=True)
wb = openpyxl.load_workbook(src, read_only=True, data_only=True)
ws = wb['按月平均']
with out.open('w', newline='', encoding='utf-8-sig') as handle:
    csv.writer(handle).writerows(ws.values)
print(out)
