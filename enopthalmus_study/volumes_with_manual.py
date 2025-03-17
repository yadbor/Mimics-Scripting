import os
import re
import csv

import mimics

import result_logger

from const import * # Safe to use * as only contains CONST variables
import utils # My utility functions
import materials # Contains definitions of all materials for analysis
# Segment orbital contents into these Materials & measure volumes
# Define a new material in materials.py and add here to include in the analysis.
orbit_materials = {
  'air'    : materials.MATL_AIR, 
  'fat'    : materials.MATL_FAT,
  'muscle' : materials.MATL_MUSCLE
}

# ## Extracting info from the project
# # Useful (?) information in the mimics project
# project_info_fields = ("height, width, "
#                        "slice_increment, slice_thickness, number_of_slices, "
#                        "obliqueness, orientation, algorithm, gantry_tilt, "
#                        "pixel_size, pixel_units",

#                        "patient_id, patient_name, project_path, study_date"
#                        )
# info_fields = dict(zip(('image', 'study'), [re.split(
#     r',\s*', name) for name in project_info_fields]))

# def extract_info(info, fields):
#     return (dict([(att, getattr(info, att)) for att in fields]))

# def gather_project_info():
#   # Extract project information
#   project_info = mimics.file.get_project_information()
#   info_dict = {k: extract_info(project_info, v)
#                 for k, v in info_fields.items()}

#   # Add the DoB, which is only in the DICOM tags
#   # DICOM tags come as a dict already, but the parts need decoding
#   # and different studies may have very different tags
#   dicom_tags = mimics.get_dicom_tags()
#   t = dicom_tags[0x0010,0x0030].value
#   # Add to the study info as an ISO date YYYY-MM-DD
#   info_dict['study']['DOB'] = '-'.join((t[0:4], t[4:6], t[6:8]))    

#   return info_dict, dicom_tags

# def write_results(study_info, input_info, volumes, results_file):
#   # Collapse the study info to a single dict for logging, as don't care whether they are image or subject info
#   collapsed_study = {k: v for d in study_info.values() for k, v in d.items()}
#   # Collapse the input information into a single dict for logging, as side is encoded in the label
#   collapsed_inputs = {k: v for d in input_info.values() for k, v in d.items()}
#   # Collapse the volume data into a single dict for logging, encoding side (from the dict) in the label
#   collapsed_volumes = {side + "|" + k: v for side, d in volumes.items() for k, v in d.items()}
  
#   # Combine the two sets for logging
#   combined_for_log = {**collapsed_study, **collapsed_inputs, **collapsed_volumes}
#   # The header will be all the keys in this dict, and the data will be all the values
#   headers = list(combined_for_log.keys())
#   results = list(combined_for_log.values())

#   # On first call this will create the file and write the headers, then the results.
#   # Subsequent calls wil only write the results.
#   result_logger.log_to_file(results_file, headers, results)

def mask_volume(mask):
  if mask:
    vol = mask.volume
  else:
    vol = 0
  return vol

def measure_project():
  sides = ['left', 'right']
  # First ensure that the globe has been subtracted from each orbital volume
  for side in sides:
    # Select the orbital volume to analyse
    orbit_name = f"{side}_Orbital Volume"
    temp = mimics.data.masks.find(orbit_name)
    # Make sure the globe has been removed from the orbital volume
    globe = utils.sphere_to_mask(mimics.data.spheres[f"{side}_globe"])
    orbit = mimics.segment.boolean_operations(temp, globe, operation="Minus")
    orbit.name = orbit_name
    mimics.data.masks.delete([temp, globe])
    volumes[side] = {'orbital' : orbit.volume}

  # Find the bounding box for both orbits 
  orbits = mimics.data.masks.filter("(?i)Orbital Volume$", regex=True)
  bbox_both = mimics.measure.get_bounding_box(orbits)
  
  # For each material to analyse
  for material in orbit_materials:
    # Create a mask for that material which covers both orbits
    material_mask = utils.mask_from_material(material, orbit_materials[material], bbox_both)
    # then intersect the material with each orbit
    for side in sides:
      orbit_name = f"{side}_Orbital Volume"
      orbit = mimics.data.masks.find(orbit_name)
      temp_mask = mimics.segment.boolean_operations(orbit, material_mask, operation="Intersect")
      temp_mask.name = f"{side}_{material}"
      
      volumes[side][material] = temp_mask.volume
    # Clean up the material mask
    mimics.data.masks.delete(material_mask)

  # Check for manual segmentation masks and return their volume
  # The manual masks are called "Manual_Muscle+Nerve", "Manual_Haematoma" and "Haematoma_in_Orbital_Volume" (mosly - some variations)
  manual_muscle = mimics.data.masks.find("(?i)Manual.?Muscle.*Nerve.*", regex=True)
  volumes['manual']['muscle'] = mask_volume(manual_muscle)

  manual_haematoma = mimics.data.masks.find("(?i)Manual.?Haematoma.*", regex=True)
  volumes['manual']['haematoma'] = mask_volume(manual_haematoma)

  manual_haematoma_orbit = mimics.data.masks.find("(?i)Haematoma.*Orbit.*", regex=True) 
  volumes['manual']['haematoma_orbit'] = mask_volume(manual_haematoma_orbit)
  
  return volumes

def write_results(results_file, volumes):
  # Collapse the volume data into a single level dict for logging, 
  # encoding side (from the dict) in the label
  collapsed = {side + "|" + k: v for side, d in volumes.items() for k, v in d.items()}
  # The header will be all the keys in this dict, and the data will be all the values
  headers = list(collapsed.keys())
  results = list(collapsed.values())
  # On first call this will create the file and write the headers, then the results.
  # Subsequent calls wil only write the results.
  print(f"writing\n{results}")
  result_logger.log_to_file(results_file, headers, results)

if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.

  root = r'D:\Projects & Research\Enophthalmos Study\to_measure'

  results_file = result_logger.Path(os.path.join(root, 'volumes.csv'))
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
      # set up blank volumes dict
      volumes = {'project': {'name': project_name[:-4]}} # cut the .mcs extension
      for side in ['left', 'right']:
        volumes[side] = {mat:0 for mat in orbit_materials}
      volumes['manual'] = {part : 0 for part in ['muscle', 'haematoma', 'haematoma_orbit']}

      try:
        volumes = measure_project()
        # Cut the .mcs extension then add a suffix and put the extension back
        measured_name = f"{p[:-4]}_measured.mcs"
        mimics.file.save_project(measured_name)
      except KeyboardInterrupt:
        exit
      except Exception as e:
        print("Error measuring volumes: ", e)
      
      write_results(results_file, volumes)
      mimics.file.close_project()
  
  print("Finished.")