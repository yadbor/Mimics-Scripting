import os
import openpyxl
import pandas


root = 'D:\Projects & Research\Enophthalmos Study'
sheet = 'files_to_analyse.xlsx'

fullname = os.path.join(root, sheet)

tbl = pandas.read_excel(fullname)
# Drop rows with no file and empty columns
tbl = tbl.dropna(axis=0, subset='file').dropna(axis=1, how='all')

# Now have the clean table with only the used columns and rows
