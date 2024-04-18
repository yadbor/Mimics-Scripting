# Analyse all Mimics projects in given folder
# Enophthalmous study 2024

import os # for scandir() etc
import re # for regexp matching

import numpy as np # Only used for np.array() in bounding box construction

from const import * # CONSTANT definitions for this project (* is safe as only consts)
import utils # Utility functions to simplify code

from orbital_analysis import make_orbit_mask
from result_logger import log_to_file         # Write measurements to a .CSV file 

import mimics # API to access mimics

import materials # Contains definitions of all materials for analysis
# Segment orbital contents into these Materials & measure their volumes
# Define a new material in materials.py and add here to include in the analysis.
orbit_materials = {
  'air'    : materials.MATL_AIR, 
  'fat'    : materials.MATL_FAT,
  'muscle' : materials.MATL_MUSCLE
}

# Useful (?) information in the mimics project
project_info_fields = ("height, width, "
                       "slice_increment, slice_thickness, number_of_slices, "
                       "obliqueness, orientation, algorithm, gantry_tilt, "
                       "pixel_size, pixel_units", 
                       "patient_id, patient_name, project_path, study_date"
                       )
info_fields = dict(zip(('image', 'study'), [re.split(r',\s*', name) for name in project_info_fields]))

def extract_info(info, fields):
  return(dict([(att, getattr(info, att) ) for att in fields]))

# Ge tthe names of all projects analyzed by Sam
root = r'D:\Projects & Research\Enophthalmos Study'
projects = [f for f in os.scandir(root) if re.match(r'.*SS\s+\d+\.mcs', f.name)]

for p in projects:
  mimics.file.open_project(filename = p.path, read_only_mode=True)
  
  # Extract project information
  project_info = mimics.file.get_project_information()
  info_dict = {k:extract_info(project_info, v) for k, v in info_fields.items()}

  # 1. make mask_bone and part_bone

  # Ignore any previously defined bone mask and make a new one 
  # to use the fill_holes and keep_largest options
  mask_bone = utils.mask_from_material("bone", materials.MATL_BONE)
  mimics.segment.fill_holes(mask_bone)
  mimics.segment.keep_largest(mask_bone)
  # Create a matching Part
  part_bone = utils.part_from_mask(mask_bone)

  # Also need a background mask for *everything* that is air
  mask_air = utils.mask_from_material("Air Mask", materials.MATL_AIR)

  # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
  #    X+ to the patient's Left
  #    Y+ to the patient's Posterior
  #    Z+ to the patient's Superior
  # Determine which geometry is on which side, based JUST on the location of the eyeballs
  num_eyes = len(mimics.data.spheres)
  if num_eyes == 2:
    # Look for eye on right hand side of the scan. If not it must be on the left.
    if mimics.data.spheres[0].center[X] < mimics.data.spheres[1].center[X]:
      eyes = 'right', 'left'
    else:
      eyes = 'left', 'right'
  elif num_eyes == 1:
    # No eye to compare with, so compare against approximate centreline
    if mimics.data.spheres[0].center[X] < 0:
      eyes = 'right',  # note the trailing comma to make this a singleton tuple
    else:
      eyes = 'left', 
  else:
    print(f'Wrong number of eyes! Expects one ot two, have {num_eyes}')
    break # ERROR - too many or too few eyes, so move to next project

  # This could be more robust, as it assumes that each eye component is measured
  # in the same order - i.e. all the left components before the right components
  # So could do left rim, left globe, left point, right rim, right globe, right point
  # or left rim, right rim, left globe, right globe, left point, right point, etc.
  # Could test each indiviual component and either save in a dict or 
  # name them in mimics.data.{spheres|splines|points}

  volumes = {} # Create a blank dict to hold measured volumes for each eye
  # Analyze each eye in the project
  for eye, side in enumerate(eyes):
    rim = mimics.data.splines[eye]
    globe = mimics.data.spheres[eye]
    point = mimics.data.points[eye]
    side_label = side + ' ' # add a trailing space to make nice labels
    
    # Get the rim extents 
    rim_geometry = utils.spline_geometry(rim)
    min_pt = rim_geometry.min_point
    max_pt = rim_geometry.max_point
  
    # make an orbit mask, given the rim and globe
    mask_union = make_orbit_mask(rim, globe)

    # Union the orbit mask with Bone and Air masks
    mask_union = utils.masks_unite(mask_union, mask_bone)
    mask_union = utils.masks_unite(mask_union, mask_air)



    # Fill the combined masks
    mask_smartfill = mimics.segment.smart_fill_global(mask_union, 7)
  
    # Expand the bbox beyond the rim by +/- 10 X and +/- 15 Z
    # with Y extending -5 and out 80. Could go to max_pt[Y] + 80
    expand = (10, 5, 15) # I think X should be more

    box_orig = np.array(min_pt) - np.array(expand)
    vector_x = (max_pt[X] - min_pt[X] + (2 * expand[X]), 0, 0)
    vector_y = (0, 80, 0)
    vector_z = (0, 0, max_pt[Z] - min_pt[Z] + (2 * expand[Z]))
    
    orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)

    # Crop the filled mask and create a part
    mask_smartfill = mimics.segment.crop_mask(mask_smartfill, orbit_ROI)
    part_smartfill = utils.part_from_mask(mask_smartfill)
    
    # Here could use code to improve the bone surface, as at
    # https://github.com/Pythonsegmenter/Orbital-floor-maxillar-sinus-reconstruction
    # also at https://gist.github.com/Pythonsegmenter/bf9c0df43c9d6260e35ad4b786faf90c
    
    part_wrap = mimics.tools.wrap(part_smartfill, 0.2, 10, False, True, True)
    mask_wrapped = mimics.segment.calculate_mask_from_part(part_wrap, None)

    orbit_vol = utils.masks_subtract(mask_wrapped, mask_smartfill)
    orbit_vol = mimics.segment.morphology_operations(orbit_vol, 'Erode', 1, 8, None, None)
    orbit_vol = mimics.segment.region_grow(orbit_vol, orbit_vol, globe.center, 'Axial', False, True, connectivity='6-connectivity') 
    orbit_vol.name = side_label + "Orbital Volume"

    # Make a mask from the globe
    mask_globe = utils.sphere_to_mask(globe)
    # Subtract the globe from the orbit mask
    mask_intersect_vol = utils.masks_subtract(orbit_vol, mask_globe)
    mask_intersect_vol.name = side_label + "Intersect Mask"
    # and make into a Part
    part_orbital_vol = utils.part_from_mask(mask_intersect_vol)

    # Create a list of masks and corresponding parts for each material 
    # in the orbit_materials dict.
    # Put the orbital volume first in the parts list
    parts = {'orbital': part_orbital_vol}
    masks = {'orbital': mask_intersect_vol}
    for matl in orbit_materials:
      masks[matl] = utils.mask_from_material(matl + ' mask', orbit_materials[matl])
      masks[matl] = utils.masks_intersect(masks[matl], mask_intersect_vol)
      if masks[matl].number_of_pixels > 0:
        parts[matl] = utils.part_from_mask(side_label + matl, masks[matl])

    
    # # Delete all masks (except the Orbital Volume mask...)
    # objects_to_keep = ("Bone", "Spline 1", "Sphere 1", "Orbital Volume")
    # for m in mimics.data.objects:
    #   if m.name in objects_to_keep:
    #     print(f"Keeping {m.name} {m.type}")
    #   else:
    #     mimics.data.objects.delete(m)

    # extract the results for this eye
    volumes[side] = {side_label + name: part.volume for name, part in parts.items()}

  # Write the results from this project
  
  mimics.file.close_project # move to the next project
  
# Finished