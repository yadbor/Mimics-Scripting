<<<<<<< HEAD
# Analyse all Mimics projects in given folder
# Enophthalmous study 2024

import logging
import os  # for scandir() etc
import re  # for regexp matching
import time  # for timeing processes

import numpy as np  # Was only used for np.array() in bounding box construction; now for several more

# CONSTANT definitions for this project (* is safe as only consts)
from const import *
import utils  # Utility functions to simplify code

# from orbital_analysis import make_orbit_mask, make_orbit_ROI
import orbital_analysis
# Write measurements to a .CSV file
from result_logger import Path, log_to_file

import mimics  # API to access mimics

import materials  # Contains definitions of all materials for analysis
# Segment orbital contents into these Materials & measure their volumes
# Define a new material in materials.py and add here to include in the analysis.
orbit_materials = {
    'air': materials.MATL_AIR,
    'fat': materials.MATL_FAT,
    'muscle': materials.MATL_MUSCLE
}

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


# Get the names of all projects analyzed by Sam
root = r'D:\Projects & Research\Enophthalmos Study'
projects = [f.path for f in os.scandir(root) if re.match(r'.*SS\s+\d+\.mcs', f.name)]

TESTING = False
# For testing
if TESTING == True:
    projects = [f.path for f in os.scandir(r'D:\Projects & Research\Enophthalmos Study\test data')
                if re.match(r'.*for_test.mcs', f.name)]

# Just the chosen good patients
good = ('Evas', 'Dabbs', 'Hayes', 'Madden', 'Umbulgarri', 'Rogers')
re_good = re.compile('('+'|'.join(good)+').*\\d{2}\\.mcs$', flags=re.IGNORECASE)
projects = [f.path for f in os.scandir(root) if re.match(re_good, f.name)]

## Just do the extra scan Sam digitised 2024-05-14
#projects = [r'D:/Projects & Research/Enophthalmos Study/Umbulgarri 2023.12.24 Merged SS 01.mcs']

# Count the files we are processing
num_projects = len(projects)

# Define a helper function to return mask information as a formatted string for printing etc.
def mask_info(mask):
    info = f'{mask.name} {mask.number_of_pixels} pixels from {mask.minimum_value} to {mask.maximum_value}\n {mimics.measure.get_bounding_box(mask)}'
    return info

# Create a results file in the same folder as the project files
results_file = Path(os.path.join(os.path.dirname(projects[0]), 'results.csv'))

# Set up logging in the same folder as the project files
logger = logging.getLogger(__name__)
log_file = os.path.join(os.path.dirname(projects[0]), 'process_studies.log')
logging.basicConfig(filename=log_file, level=logging.DEBUG)


# ## Test with one scan ##
# projects = [r'D:/Projects & Research/Enophthalmos Study/Umbulgarri 2023.12.24 Merged SS 01.mcs']
#projects = projects[4] # just one file that has failed previously

# Time events
events = utils.Events() # Initialise an event timer. There is a seperate one in segment_orbit
events.add('started loop')

logger.info('#####' * 20)

for i, p in enumerate(projects):
    project_start = time.process_time()
    msg = f'\n\nstarted project {i} of {num_projects}:  {p} at {round(events.elapsed(), 2)}.\n'
    logger.info(msg)
    print(msg)

    mimics.file.open_project(filename=p, read_only_mode=True)
    events.add('opened project')

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

    # # 1. make mask_bone and part_bone
    # # Ignore any previously defined bone mask and make a new one
    # # to use the fill_holes and keep_largest options.
    # # Use duplicate() to copy this mask later
    # mask_bone = mimics.segment.keep_largest(mimics.segment.fill_holes(
    #     utils.mask_from_material('bone', materials.MATL_BONE)))
    # print(mask_info(mask_bone))
    # # Also need a background mask for *everything* that is air (maybe?)
    # mask_air = utils.mask_from_material('Air Mask', materials.MATL_AIR)
    # print(mask_info(mask_air))

    # Need to check why I added this. Looks like they _should_ work.
    # if project_info.orientation != 'RAB':
    #     logger.exception(f"Can't handle image orientation {project_info.orientation}. Abort this project.")
    #     events.add('skip project - wrong orientation')
    #     continue # Go to next project

    num_eyes = len(mimics.data.spheres)
    num_rims = len(mimics.data.splines)
    num_pts = len(mimics.data.points)

    if num_eyes != num_rims:
        print(
            f'Number of globes {num_eyes} does not match number of rims {num_rims}.')
        logger.error(
            f'Aborting project - number of globes {num_eyes} does not match number of rims {num_rims}.')
        mimics.file.close_project()  # move to the next project
        continue  # ERROR - too many eyes or too few eyes, so move to next project

    if num_eyes > 2 or num_eyes < 1:
        print(f'Wrong number of eyes! Expected one or two, found {num_eyes}')
        logger.error(
            f'Aborting project - expected one or two eyes, found {num_eyes}.')
        mimics.file.close_project()  # move to the next project
        continue  # ERROR - too many eyes or too few eyes, so move to next project

    # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
    #    X+ to the patient's Left, so -ve to the Right
    #    Y+ to the patient's Posterior, so -ve Anterior
    #    Z+ to the patient's Superior, so -ve Inferior
    # For each eye component, determine which geometry is on which side, based on the location
    # If there are two of a given component, compare them; otherwise compare to 0
    if num_eyes == 1:
        # Just need to determine the side once, and can use the globe for that
        if mimics.data.spheres[0].center[X] < 0:
            side = ('right', 'left')
        else:
            side = ('left', 'right')
        # Note that some may not have an apex point
        eyes = {side[0]: {'globe': mimics.data.spheres[0],
                          'rim':   mimics.data.splines[0],
                          'point': mimics.data.points[0] if len(mimics.data.points) > 0 else None
                            },
                # Empty entry to fill out results table
                side[1]: {'globe': None, 'rim': None, 'point': None}
                }
    else:
        # num_eyes must be 2 here.
        # Need to find the side of each component indivdually, as could be entered in random order.
        # Find which index is left and which is right for each component by checking the X coordinate.

        # Set up a blank dict
        eyes = {
            'right': {'globe': None, 'rim': None, 'point': None},
            'left':  {'globe': None, 'rim': None, 'point': None}
            }

        # globes - centers of the first two spheres
        temp = [o.center[X] for o in [mimics.data.spheres[i] for i in (0,1)]]
        idx = 0 if temp[0] < temp[1] else 1
        eyes['right']['globe'] = mimics.data.spheres[idx] # idx is either 0 or 1
        eyes['left']['globe']  = mimics.data.spheres[1 - idx] # Opposite of idx 
       
        # rims - centroids of the first two splines
        temp = [utils.spline_center(o)[X] for o in [mimics.data.splines[i] for i in (0,1)]]
        idx = 0 if temp[0] < temp[1] else 1
        eyes['right']['rim'] = mimics.data.splines[idx] # idx is either 0 or 1
        eyes['left']['rim']  = mimics.data.splines[1 - idx] # Opposite of idx 

        # points - more complicated as may be missing one or both
        # Have set up the dict with with both missing, so just add any that are found
        if num_pts == 1:
            side = 'right' if mimics.data.points[0][X] < 0 else 'left'
            eyes[side]['point'] = mimics.data.points[0]
        elif num_pts >= 2:
            temp = [o[X] for o in [mimics.data.points[i] for i in (0,1)]]
            idx = 0 if temp[0] < temp[1] else 1
            eyes['right']['point'] = mimics.data.points[idx] # idx is either 0 or 1
            eyes['left']['point']  = mimics.data.points[1 - idx] # Opposite of idx 

    # Name the inputs in mimics.data.{spheres|splines|points}, so can find them via name
    for side, d in eyes.items():
        for part, obj in d.items():
            obj.name = side + '_' + part
            obj.visible = True

    volume_start = time.process_time()

    # Create a blank dict to hold measured volumes for each eye
    volumes = {}
    # And a blank dict to hold summary input information (spline bounds, globe & point location)
    input_info = {}
    # Analyze each eye in the project
    for side, eye_parts in eyes.items():
        rim = eye_parts['rim']
        globe = eye_parts['globe']
        point = eye_parts['point']

        side_label = side + '|'  # add a seperator to enable splitting out the side later

        ############################################################################################
        # Need to do something here to deal with single eye studies.
        # If parts are None for this eye skip calculations & add 0 volumes to keep order of results correct.
        ############################################################################################

        msg = f'\n########## segment {side} orbit at {round(events.elapsed(), 2)}.\n'
        logger.info(msg)
        print(msg)

        events.add(f'start segment_orbit on {side}')
        # Orbit volume mask extraction moved to orbital_anaylsis.py
        try:
            with mimics.disabled_gui_update():
                orbit_vol, orbit_ROI = orbital_analysis.segment_orbit(rim, globe, point, side)
            orbit_vol.name = f'{side}_Orbital Volume'
            print(mask_info(orbit_vol))
        except Exception as e:
           print("Failed to segment the orbit, with error: ", e)
           logger.exception(f"Failed to segment the orbit, with error: {e}. Aborting this project.")
           events.add('skip project - error in segment_orbit')
           print("skipping this project")
           continue # Go to next project

        # Skip deleting old masks

        # Convert globe (which is a Sphere) to a mask
        m_globe = utils.sphere_to_mask(globe)
        m_globe.name = f'{side}_m_globe'
        print(mask_info(m_globe))

        m_intersect_vol = mimics.segment.boolean_operations(orbit_vol, m_globe, 'Minus')
        m_intersect_vol.name = f'{side}_m_intersect_vol'
        print(mask_info(m_intersect_vol))
        check_volume = m_intersect_vol.volume
        # and make into a Part
        p_orbital_vol = mimics.segment.calculate_part(m_intersect_vol, quality='High')
        p_orbital_vol.name = f'{side}_p_orbital_vol'

        orbit_vol_ROI = mimics.measure.get_bounding_box(m_intersect_vol)  # Material masks to insersect only need to be this bigs

        events.add('made orbital volume & globe')

        # Create a list of masks and corresponding parts for each material
        # in the orbit_materials dict. Put the orbital volume first as it should always exist.
        masks = {'orbital': m_intersect_vol}
        parts = {'orbital': p_orbital_vol}
        for matl in orbit_materials:
            # Mask of where this material overlaps with intersect_vol_mask
            masks[matl] = utils.mask_from_material('m_' + side + '_' + matl, orbit_materials[matl], bounding_box=orbit_vol_ROI)
            masks[matl] = utils.masks_intersect(masks[matl], m_intersect_vol)
            masks[matl].name = f'{side}_m_{matl}'

        # Can get volume directly from mask, but it's different to parts, so use both methods to compare.
        # Initialise the vols dict using the mask & part that should always have a volume.
        # This dict has two items: a dict of volumes from the masks and one from the parts (where they exist)
        vols = {'mask': {'orbital': m_intersect_vol.volume},
                'part': {'orbital': p_orbital_vol.volume}}
        # Can't just make everything a part as some masks may be empty, so catch that.
        for name, mask in masks.items():
            vols['mask'][name] = mask.volume
            if mask.number_of_pixels == 0:
                parts[name] = None
                vols['part'][name] = 0
            else:
                parts[name] = utils.part_from_mask(side_label + name, mask)
                vols['part'][name] = parts[name].volume

        # Store the results for this eye, adding the source (mask or part) to each material label
        volumes[side] = {source + '_' + k: v for source,
                         d in vols.items() for k, v in d.items()}

        # Store the input info for this eye
        bbox_rim = mimics.measure.get_bounding_box([rim])
        p1, p2 = utils.bbox_to_points(bbox_rim)
        input_info[side] = {
            **utils.labelled_point(prefix=side_label, name='geom_rim1', point=p1),
            **utils.labelled_point(prefix=side_label, name='geom_rim2', point=p2),
            **utils.labelled_point(prefix=side_label, name='geom_globe', point=globe.center),
            side_label + 'geom_radius': globe.radius,
            **utils.labelled_point(prefix=side_label, name='geom_apex', point=point),
            side_label + 'geom_rim.w': p2[X] - p1[X],
            side_label + 'geom_rim.d': p2[Y] - p1[Y],
            side_label + 'geom_rim.h': p2[Z] - p1[Z]
        }

        events.add(f'{side} eye done')


    # Done all sides. Write the results from this project

    # Collapse the study info to a single dict for logging, as don't care whether they are image or subject info
    collapsed_study = {k: v for d in info_dict.values() for k, v in d.items()}
    print(collapsed_study)
    # Collapse the input information into a single dict for logging, as side is encoded in the label
    collapsed_inputs = {k: v for d in input_info.values() for k, v in d.items()}
    # Collapse the volume data into a single dict for logging, encoding side (from the dict) in the label
    collapsed_volumes = {side + "|" + k: v for side, d in volumes.items() for k, v in d.items()}
    print(collapsed_volumes)
    # Combine the two sets for logging
    combined_for_log = {**collapsed_study, **collapsed_inputs, **collapsed_volumes}
    # The header will be all the keys in this dict, and the data will be all the values
    headers = list(combined_for_log.keys())
    results = list(combined_for_log.values())

    # On first call this will create the file and write the headers, then the results.
    # Subsequent calls wil only write the results.
    log_to_file(results_file, headers, results)

    events.add('volumes done')

    print(f'elapsed time for this project was {events.elapsed()}')
    
    # mimics.file.open_project(filename=p.path, read_only_mode=True)
    processed_path = re.sub(f'\.mcs', '.processed.mcs', project_info.project_path)
    mimics.file.save_project(filename=processed_path,)
    print('Closing Project')
    mimics.file.close_project()  # move to the next project
    events.add('project closed.\n')

events.add('all projects done')
print(f'total elapsed time for {len(projects)} projects was {events.elapsed()}')
logging.shutdown()
# Finished
=======
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
    mask_union = utils.masks_unite(mask_union, mask_not_orbit)

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
    mask_orbit_vol = utils.masks_subtract(mask_wrapped, mask_smartfill)
    # Erode that to separate the orbit volume from the surroundings
    mask_orbit_vol = mimics.segment.morphology_operations(mask_orbit_vol, 'Erode', 1, 8, None, None)
    # Grow the separated region to get a mask of only the orbital contents, starting from the centre of the globe
    mask_orbit_vol = mimics.segment.region_grow(mask_orbit_vol, mask_orbit_vol, globe.center, 'Axial', False, True, connectivity='6-connectivity') 
    mask_orbit_vol.name = side_label + "Orbital Volume"
  
    # Make a mask from the globe
    mask_globe = utils.sphere_to_mask(globe)
    # Subtract the globe from the orbit mask
    mask_intersect_vol = utils.masks_subtract(mask_orbit_vol, mask_globe)
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
mask_bone_boolean = utils.masks_unite(mask_temp_morph, mask_bone)
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
mask_intersect = utils.masks_intersect(mask_temp_orbit, subtraction_mask)
#
mask_orbit_volume = mimics.segment.region_grow(mask_intersect, None, sphere1.center, "Axial", keep_original_mask=False, multiple_layer=True, connectivity='6-connectivity')
mask_orbit_volume.name = "Orbit Volume"

>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
