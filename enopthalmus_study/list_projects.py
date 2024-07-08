# Measure the intersection volumes for various tissues, 
# given folders of projects with already segmented orbits.
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching
import numpy as np # array functions and vector maths

# CONSTANT definitions for this project (* is safe as only consts)
from const import *

import utils # Utility functions to simplify code

import mimics # API to access mimics

# Write measurements etc. to a .CSV file.
from result_logger import Path, log_to_file

def mimics_info_as_dict(info):
    """Return all informative fields in a mimics.ImageInformation object."""
    return dict([(a, getattr(info, a)) for a in dir(info) if a[:2] != "__"])

def image_get_info(img):
    info = img.get_image_information()
    info_dict = {
        "algorithm":   info.algorithm,
        "obliqueness":  info.obliqueness,
        "orientation":  info.orientation,
        "gantry_tilt":  info.gantry_tilt,
        "slice_thickness":  info.slice_thickness,
        "slice_increment":  info.slice_increment,
        "pixel_size": info.pixel_size
        }
    return info_dict

if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
  root = r'D:\Projects & Research\Enophthalmos Study'
 
  # Put the combined results in the root, rather than one file per user.
  results_file = Path(os.path.join(root, 'project_info.csv'))

  # List the mimics project files to be processed (all in this case).
  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]

  num_projects = len(projects)
  for i, p in enumerate(projects): 
    try:
        # Process a single project file.
        print(f"reading project {i} of {num_projects}: {os.path.basename(p)}")

        mimics.file.open_project(filename=p, read_only_mode=True)
        
        project_info = mimics.file.get_project_information()
        filename = os.path.basename(project_info.project_path)
        study_info = {
            'filename': filename,
            'num_series': len(mimics.data.images)
            }

        # process all image sets
        for img in mimics.data.images:
            img_info = {'name': img.name, **image_set_info(img)}
            basis = (i, j, k) = mimics_basis_vectors(img)
            img_info['i'] = i
            img_info['j'] = j
            img_info['k'] = k
            
            combined = {**study_info, **img_info} 
            log_to_file(results_file, headers = list(combined.keys()), results = list(combined.values()))

        mimics.file.close_project()
        
    except (IndexError, ValueError):
        # If something went wrong processing this project close it and move to the next one.
        mimics.file.close_project()
        continue

  print('********** Finished all files **********')

