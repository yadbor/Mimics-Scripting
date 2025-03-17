import os
import datetime
import re
import csv

import mimics

# def process_row(r):
#     mimics.file.open_project(r['file'], read_only_mode=True)
#     mask_name = [m.name for m in mimics.data.masks]
#     print(f"masks")

# root = 'D:\Projects & Research\Enophthalmos Study'
# sheet = 'files_to_analyse.xlsx'

# fullname = os.path.join(root, sheet)

# tbl = pandas.read_excel(fullname)
# # Drop rows with no file and empty columns
# tbl = tbl.dropna(axis=0, subset='file').dropna(axis=1, how='all')

# # Now have the clean table with only the used columns and rows, procees each row
# nrows = tbl.shape[0]
# for i, r in tbl.iterrows():
#     print(f"process file {i+1}/{nrows} {r['file']}")
#     process_row(r)

def directory_spider(input_dir, path_pattern="", file_pattern="", maxResults=500):
    file_paths = []
    if not os.path.exists(input_dir):
        raise FileNotFoundError("Could not find path: %s"%(input_dir))
    for dirpath, dirnames, filenames in os.walk(input_dir):
        if re.search(path_pattern, dirpath):
            file_list = [item for item in filenames if re.search(file_pattern,item)]
            file_path_list = [os.path.join(dirpath, item) for item in file_list]
            file_paths += file_path_list
            if len(file_paths) > maxResults:
                break
    return file_paths[0:maxResults]

def file_mod_time(f):
    return datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime('%Y-%m-%d_%H:%M:%S')

def parse_project_name(project_name):
  if re.match(pattern=".*_orbit.mcs", string = project_name):
    (stem, date, series, extra) = project_name.split('_')
    return (stem, date, series, "orbit", "")

  try: 
    (stem, date, series, extra) = re.match(pattern="(.*) ([0-9\\.]+) ([A-Za-z]+) SS 01(.*).mcs", string=project_name).groups()
  except AttributeError:
        return ("", "", "", "", "")

  state = ""
  if re.match(pattern="(?i)processed", string=extra):
    state = "processed"
  operator = ""
  if re.match(pattern="(?i)rob", string=extra):
    operator = "rob"
  if re.match(pattern="(?i)ryan", string=extra):  
    operator = "ryan"

  series = series.casefold()
  if series[0:5] == "merge":
    series = "merged"

  return (stem, date, series, state, operator)

# Clean masks
masks_to_keep = [m.name for m in mimics.data.masks.filter("(?i)Manual|Haematoma|Volume|globe", regex=True)]
for m in mimics.data.masks:
    if not m in masks_to_keep:
      mimics.data.masks.delete(m)

# Clean objects
for p in mimics.data.planes:
    mimics.data.planes.delete(p)
# Delete points, keeping the apex points
for i, p in enumerate(mimics.data.points):
    if not re.match(pattern="^(left|right)", string=p.name):
        mimics.data.points.delete(p)

if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  #root = r'D:\Projects & Research\Enophthalmos Study\re-do_DICOM'
  root = r'D:\Projects & Research\Enophthalmos Study'

  csv_file_name = os.path.join(root, 'all_masks.csv')
  
  # Get a list of all the .mcs files in root
#   projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]
#   num_projects = len(projects)
  
  # Get a list of all the .mcs files in all folder under root
  projects = directory_spider(root, file_pattern = ".*mcs", maxResults=1000)
  num_projects = len(projects)
  
  with open(csv_file_name, 'w', newline='') as csvfile:
    maskwriter = csv.writer(csvfile)
    for i, p in enumerate(projects):
        print(f"process {i+1}/{num_projects}\t{p}")

        project_name = os.path.basename(p)
        stem, date, series, state, operator = parse_project_name(project_name)
        if re.match(pattern="Ryan Processing", string=p):
            operator = "ryan"
        if re.match(pattern="Rob Processing", string=p):
            operator = "rob"

        if stem == "":
            print(f"skipped {project_name} as not in study.")
            continue

        mimics.file.open_project(filename=p, read_only_mode=True)
        masks = [m.name for m in mimics.data.masks]
        
        op_row = [p, file_mod_time(p), project_name, stem, date, series, state, operator, *masks]
        maskwriter.writerow(op_row)

print("Finished.")