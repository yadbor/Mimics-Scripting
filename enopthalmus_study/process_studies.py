# Analyse all Mimics projects in given folder
# Enophthalmous study 2024

import os # for scandir() etc
import re # for regexp matching

import numpy as np # Only used for np.array() in bounding box construction

from const import * # CONSTANT definitions for this project (* is safe as only consts)
import utils # Utility functions to simplify code

#from orbital_analysis import make_orbit_mask, make_orbit_ROI
import orbital_analysis
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

def find_eyes():
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


# Get the names of all projects analyzed by Sam
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
    
    # Mask off everything anterior to the orbital rim, given the rim and globe
    mask_union = orbital_analysis.make_orbit_mask(rim, globe)
    # Make an ROI based on the rim, extending past that to cover the whole orbit
    orbit_ROI = orbital_analysis.make_orbit_ROI(rim)

    # Crop the air and bone masks to this ROI and unite them to give a "not orbit contents" mask.
    mask_bone_ROI = mimics.segment.crop_mask(mask_bone, orbit_ROI)
    mask_air_ROI = mimics.segment.crop_mask(mask_air, orbit_ROI)
    mask_not_orbit = utils.mask.unite(mask_bone_ROI, mask_air_ROI)
    mimics.data.masks.delete(mask_bone_ROI) # Clean up the list of masks
    mimics.data.masks.delete(mask_air_ROI) # Clean up the list of masks
    
    # Union the orbit mask with Bone and Air masks
    mask_union = utils.unite(mask_union, mask_not_orbit)

    # Repair Orbit Walls + Floor
    # Could use code to improve the bone surface, as at
    # https://github.com/Pythonsegmenter/Orbital-floor-maxillar-sinus-reconstruction
    # also at https://gist.github.com/Pythonsegmenter/bf9c0df43c9d6260e35ad4b786faf90c

    # Fill the combined masks to give everything *except* the orbital contents
    mask_smartfill = mimics.segment.smart_fill_global(mask_union, 7)
    # Crop the filled mask with the orbit ROI and convert to a part
    # Already cropped above, but redo after uniting the anterior mask
    mask_smartfill = mimics.segment.crop_mask(mask_smartfill, orbit_ROI)
    part_smartfill = utils.part_from_mask(mask_smartfill)

    # Wrap the part to filter small inclusions > 0.2 or close small holes < 10
    part_wrapped = mimics.tools.wrap(part_smartfill, 0.2, 10, False, True, True)
    #  and convert back into a mask.
    mask_wrapped = mimics.segment.calculate_mask_from_part(part_wrapped, None)
    
    # Now create a mask of everything in the ROI
    # Make a mask for all thresholds, cropped to the bounding box
    mask_temp_orbit = mimics.segment.threshold(
                         mask = mimics.segment.create_mask(), 
                         threshold_min = MIN_GV, 
                         threshold_max = MAX_GV,
                         bounding_box = orbit_ROI
                         ) #change this to landmarks

    # The orbit is the wrapped mask of everything minus the filled version ???
    mask_orbit_vol = utils.difference(mask_wrapped, mask_smartfill)
    # Erode that to separate the orbit volume from the surroundings
    mask_orbit_vol = mimics.segment.morphology_operations(mask_orbit_vol, 'Erode', 1, 8, None, None)
    # Grow the separated region to get a mask of only the orbital contents, starting from the centre of the globe
    mask_orbit_vol = mimics.segment.region_grow(mask_orbit_vol, mask_orbit_vol, globe.center, 'Axial', False, True, connectivity='6-connectivity') 
    mask_orbit_vol.name = side_label + "Orbital Volume"
  
    # Make a mask from the globe
    mask_globe = utils.sphere_to_mask(globe)
    # Subtract the globe from the orbit mask
    mask_intersect_vol = utils.difference(mask_orbit_vol, mask_globe)
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
      masks[matl] = utils.intersect(masks[matl], mask_intersect_vol)
      if masks[matl].number_of_pixels > 0:
        parts[matl] = utils.part_from_mask(side_label + matl, masks[matl])
 
    # extract the results for this eye
    volumes[side] = {side_label + name: part.volume for name, part in parts.items()}

  # Write the results from this project
  
  mimics.file.close_project # move to the next project
  
# Finished

### From version 1 (which worked?)

# Mask air once, then crop that to each eye

# Make air mask from thresholds ** cropped to a bounding box **
mask_temp_air = mimics.segment.threshold(
                       mask = mimics.segment.create_mask(), 
                       threshold_min = mimics.segment.HU2GV(-1024), 
                       threshold_max = mimics.segment.HU2GV(-200),
                       bounding_box= mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector)
                       )
# Retain the largest connected segment
mimics.segment.keep_largest(mask_temp_air)
# Make it 2 pixels bigger all around
mask_temp_morph = utils.mask_dilate(mark_temp_air, number_of_pixels = 2, connectivity = 8)
# Add to the bone mask
mask_bone_boolean = utils.unite(mask_temp_morph, mask_bone)
# Smart fill the combined mask
mask_bone_filled = mimics.segment.smart_fill_global(mask=mask_bone_boolean, hole_closing_distance=7)
part_bone_repaired = mimics.segment.calculate_part(mask_bone_filled, quality='High')
# Should be part_bone_filled, but really is part_not_orbit?

# Smooth the part and wrap it, then turn that into a mask
part_bone_smoothed = mimics.tools.smooth(part_bone_repaired, 0.5, 5, False, False)
part_bone_wrapped = mimics.tools.wrap(part_bone_smoothed, 0.2, 10, False, True, False)
part_bone_wrapped.name = "Smoothed & Wrapped Orbit"
subtraction_mask = mimics.segment.calculate_mask_from_part(part_bone_wrapped, None)
subtraction_mask.name = "Subtraction Mask"

# Make a mask for all thresholds, cropped to the bounding box
mask_temp_orbit = mimics.segment.threshold(
                         mask = mimics.segment.create_mask(), 
                         threshold_min = mimics.segment.HU2GV(-1024), 
                         threshold_max = mimics.segment.HU2GV(3071),
                         bounding_box = mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector)
                         ) #change this to landmarks
#mask_intersect = mimics.segment.boolean_operations(mask_a= mask_temp_orbit, mask_b= subtraction_mask, operation='Minus')
mask_intersect = utils.intersect(mask_temp_orbit, subtraction_mask)
#
mask_orbit_volume = mimics.segment.region_grow(mask_intersect, None, sphere1.center, "Axial", keep_original_mask=False, multiple_layer=True, connectivity='6-connectivity')
mask_orbit_volume.name = "Orbit Volume"

