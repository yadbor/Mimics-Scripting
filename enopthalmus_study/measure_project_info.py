import os
import re
import csv
import numpy as np

import mimics

import result_logger

from const import * # Safe to use * as only contains CONST variables
import utils # My utility functions

## Extracting info from the project
# Useful (?) information in the mimics project
project_info_fields = ("height, width, "
                       "slice_increment, slice_thickness, number_of_slices, "
                       "obliqueness, orientation, algorithm, gantry_tilt, "
                       "pixel_size, pixel_units",

                       "patient_id, patient_name, project_path, study_date"
                       )
info_fields = dict(zip(('image', 'study'), [re.split(
    r',\s*', name) for name in project_info_fields]))

def extract_info(info, fields):
    return (dict([(att, getattr(info, att)) for att in fields]))

def pop_dict_level(dict):
    return {k: v for d in dict.values() for k, v in d.items()}

def gather_project_info():
  # Extract project information
  project_info = mimics.file.get_project_information()
  info_dict = {k: extract_info(project_info, v)
                for k, v in info_fields.items()}

  # Add the DoB, which is only in the DICOM tags
  # DICOM tags come as a dict already, but the parts need decoding
  # and different studies may have very different tags
  dicom_tags = mimics.get_dicom_tags()
  t = dicom_tags[0x0010,0x0030].value
  # Add to the study info as an ISO date YYYY-MM-DD
  info_dict['study']['DOB'] = '-'.join((t[0:4], t[4:6], t[6:8]))    

  return info_dict, dicom_tags

def measure_project():
  '''Process a open mimics project, analysing as many eyes as exist within it.'''

  study_info, dicom_tags = gather_project_info()
 
  try:
               
    # Create dicts to hold measured volumes & inputs (spline bounds, globe & point location)
    # for each eye
    input_info = {}
    
    globes = [mimics.data.spheres[f"{side}_globe"] for side in ["left", "right"]]
    rims   = [mimics.data.splines[f"{side}_rim"]   for side in ["left", "right"]]

    for side in ["left", "right"]:
        rim = mimics.data.splines[f"{side}_rim"]
        globe = mimics.data.spheres[f"{side}_globe"]
        point = mimics.data.points[f"{side}_apex"]

        # Record the input info for this eye
        bbox_rim = mimics.measure.get_bounding_box([rim])
        p1, p2 = utils.bbox_to_points(bbox_rim)
        side_label = f"{side}_"
        input_info[side] = {
            **utils.labelled_point(prefix=side_label, name='rim_p1', point=p1),
            **utils.labelled_point(prefix=side_label, name='rim_p2', point=p2),
            **utils.labelled_point(prefix=side_label, name='globe', point=globe.center),
            side_label + 'radius': globe.radius,
            **utils.labelled_point(prefix=side_label, name='apex', point=point),
            side_label + 'rim.w': p2[X] - p1[X],
            side_label + 'rim.d': p2[Y] - p1[Y],
            side_label + 'rim.h': p2[Z] - p1[Z]
        }
    # Having processed as many eyes as exist, return the results for logging
    return study_info, input_info

  except (IndexError, ValueError):
    # Huston, we have a problem. Bail without returning results
    return


if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.

  root = r'D:\Projects & Research\Enophthalmos Study\to_measure'

  results_file = result_logger.Path(os.path.join(root, 'project_info.csv'))
  # Start with a clean empty file so that the headers are written correctly
  if os.path.exists(results_file):
    os.remove(results_file)

  # Get a list of all the .mcs files in root
  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]
  num_projects = len(projects) 
 
  for i, p in enumerate(projects):
      project_name = os.path.basename(p)
      print(f"process {i+1}/{num_projects}\t{project_name}")

      if re.match(".*_measured.mcs", p):
        print("skip already measured file")
        continue

      mimics.file.open_project(filename=p, read_only_mode=True)
      try:
        study_info, input_info = measure_project() 
        to_write = {'project': project_name[:-4], **pop_dict_level(study_info), **pop_dict_level(input_info)}

      except KeyboardInterrupt:
        exit
      except Exception as e:
        print("Error measuring distances: ", e)
      
      mimics.file.close_project()
      result_logger.log_to_file(results_file, headers = to_write.keys(), results = to_write.values())
  
  print("Finished.")