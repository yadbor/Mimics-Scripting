import os
import openpyxl
import pandas
import mimics

def process_row(r):
    mimics.file.open_project(r['file'], read_only_mode=True)
    mask_name = [m.name for m in mimics.data.masks]
    print(f"masks")

root = 'D:\Projects & Research\Enophthalmos Study'
sheet = 'files_to_analyse.xlsx'

fullname = os.path.join(root, sheet)

tbl = pandas.read_excel(fullname)
# Drop rows with no file and empty columns
tbl = tbl.dropna(axis=0, subset='file').dropna(axis=1, how='all')

# Now have the clean table with only the used columns and rows, procees each row
nrows = tbl.shape[0]
for i, r in tbl.iterrows():
    print(f"process file {i+1}/{nrows} {r['file']}")
    process_row(r)


# Clean masks
masks_to_keep = [m.name for m in mimics.data.masks.filter("Manual|Haematoma|Volume", regex=True)]
for m in mimics.data.masks:
    if not m in masks_to_keep:
      mimics.data.masks.delete(m)

