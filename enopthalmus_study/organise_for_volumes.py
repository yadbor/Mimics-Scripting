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

def parse_project_path(project_path):
  """Extract project information from the file name. Expects the full path of the current project as a string."""
  operator = "na"
  project_file = os.path.basename(project_path)

  # Look for the operator in the path if stored there
  if re.search(pattern="Ryan", string=project_path, flags=re.IGNORECASE):
            operator = "ryan"
  if re.search(pattern="Rob", string=project_path, flags=re.IGNORECASE):
            operator = "rob"

  # Special case for orbit names that have already been cleaned
  if re.search(pattern="orbit.mcs", string = project_file):
    (patient, scan_date, series, extra) = project_file.split('_')
    return (project_file, patient, scan_date, series, "orbit", operator)

  # Otherwise try to extract the information from the file name
  try: 
    (patient, scan_date, series, extra) = re.match(pattern="(.*) ([0-9\\.]+) ([A-Za-z]+) SS 01(.*).mcs", string=project_file).groups()
  except AttributeError:
    # Failed to match the pattern, so not a project we are interested in
    return (project_file, "", "", "", "", "")

  # If the pattern matched, try to deduce the state 
  state = ""
  if re.search(pattern="(?i)processed", string=extra):
    state = "processed"
  operator = ""
  if re.search(pattern="(?i)rob", string=extra):
    operator = "rob"
  if re.search(pattern="(?i)ryan", string=extra):  
    operator = "ryan"

  # Standardise the series name
  series = series.casefold()
  if series[0:5] == "merge":
    series = "merged"

  return (project_file, patient, scan_date, series, state, operator)

def clean_objects():
  # Clean masks
  masks_to_keep = [m for m in mimics.data.masks.filter("(?i)Manual|Haematoma|Orbital Volume|globe", regex=True)]
  for m in mimics.data.masks:
    if not m in masks_to_keep:
      mimics.data.masks.delete(m)

  # Clean objects
  # Delete all planes
  mimics.data.planes.delete(mimics.data.planes)
  # Delete points, keeping the apex points
  for p in mimics.data.points:
    if not re.match(pattern="^(left|right)", string=p.name):
        mimics.data.points.delete(p)

if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  #root = r'D:\Projects & Research\Enophthalmos Study\re-do_DICOM'
  root = r'D:\Projects & Research\Enophthalmos Study'

  csv_file_name = os.path.join(root, 'organised.csv')

  cleaned_folder = os.path.join(root, "to_measure")
  
  # list of all manual compare projects
  manual_compare_path = os.path.join(root, 're-do_DICOM')
  manual = [f.path for f in os.scandir(manual_compare_path) if re.match('.*.mcs', f.name)]
  # list of all the projects didgitised manually by Ryan
  ryan = [f.path for f in os.scandir(root) if re.match('.*_orbit.mcs', f.name)]
  # Combine for all that we wil measure
  projects = manual + ryan
  num_projects = len(projects)
  
  with open(csv_file_name, 'w', newline='') as csvfile:
    maskwriter = csv.writer(csvfile)
    for i, p in enumerate(projects):
      print(f"process {i+1}/{num_projects}\t{p}")

      # Extract information from the project file name
      project_name, patient, date, series, state, operator = parse_project_path(p)
      
      if patient == "":
        print(f"skipped {project_name} as not in study.")
        continue

      mimics.file.open_project(filename=p, read_only_mode=True)
      # Remove any masks or othe robject not used in the final analysis
      clean_objects()
      # list the masks that are left
      masks = [m.name for m in mimics.data.masks]
      # Are there manual measurement masks? If so, note that in the state
      if re.search(pattern= "Haematoma|Manual", string="_".join(masks), flags=re.IGNORECASE):
        state = "manual"

      # Record information about this project
      op_row = [p, file_mod_time(p), project_name, patient, date, series, state, operator, *masks]
      maskwriter.writerow(op_row)

      # Create a nice regular name and save the clean project to a designated folder
      clean_name = f"{patient}_{date}_{series}_{state}_{operator}.mcs"
      mimics.file.save_project(os.path.join(cleaned_folder, clean_name))
      mimics.file.close_project()
        
  print("Finished.")