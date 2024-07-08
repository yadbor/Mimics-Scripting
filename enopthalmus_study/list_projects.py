# Measure the intersection volumes for various tissues, 
# given folders of projects with already segmented orbits.
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching

import numpy as np  # Was only used for np.array() in bounding box construction; now for several more

# CONSTANT definitions for this project (* is safe as only consts)
from const import *
import utils  # Utility functions to simplify code

import mimics  # API to access mimics

# Write measurements to a .CSV file
from result_logger import Path, log_to_file


def mimcs_info_as_dict(info):
    """Return all informative fields in a mimics.ImageInformation object."""
    return dict([(a, getattr(info, a)) for a in dir(info) if a[:2] != "__"])

def image_set_info(img):
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

def v_hat(v):
  '''Return a normalised unit vector.''' 
  mag = sum([i**2 for i in v])**0.5
  return([i/mag for i in v])

def mimics_image_vectors(img):
  '''return the three unit vectors that descrivbe this image volume.'''
  active_img = [i for i in mimics.data.images if i.active][0]
  p0 = img.get_voxel_center([0, 0, 0])
  d = img.get_voxel_buffer().shape
  x = img.get_voxel_center([d[0]-1, 0, 0])
  y = img.get_voxel_center([0, d[1]-1, 0])
  z = img.get_voxel_center([0, 0, d[2]-1])
  span = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
  (i, j, k) = [b_i / np.linalg.norm(b_i, ord=1) for b_i in span]
  
  return (i,j,k)


if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  # This version has one folder with the segmenting person int he file name
  
  root = r'D:\Projects & Research\Enophthalmos Study'
  # Put the combined results in the root, rather than one file per user
  results_file = Path(os.path.join(root, 'project_info.csv'))

  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]

  num_projects = len(projects)
  for i, p in enumerate(projects): 
    try:
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
            basis = (i, j, k) = mimics_image_vectors(img)
            img_info['i'] = i
            img_info['j'] = j
            img_info['k'] = k
            
            combined = {**study_info, **img_info} 
            log_to_file(results_file, headers = list(combined.keys()), results = list(combined.values()))

        mimics.file.close_project()
        
    except (IndexError, ValueError):
        mimics.file.close_project()  # close the currently open project
        # Move to the next project
        continue

  print('********** Finished all files **********')

