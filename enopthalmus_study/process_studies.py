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

## Just do the extra scan Sam digitised 2024-05-14
projects = [r'D:/Projects & Research/Enophthalmos Study/Umbulgarri 2023.12.24 Merged SS 01.mcs']

# Just the chosen good patients
good = ('Evas', 'Dabbs', 'Hayes', 'Madden', 'Umbulgarri', 'Rogers')
re_good = re.compile('('+'|'.join(good)+').*\\d{2}\\.mcs$', flags=re.IGNORECASE)
projects = [f.path for f in os.scandir(root) if re.match(re_good, f.name)]

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
            eyes[side]['point'] = mimics.data.spheres[0]
        elif num_pts >= 2:
            temp = [o[X] for o in [mimics.data.points[i] for i in (0,1)]]
            idx = 0 if temp[0] < temp[1] else 1
            eyes['right']['point'] = mimics.data.splines[idx] # idx is either 0 or 1
            eyes['left']['point']  = mimics.data.splines[1 - idx] # Opposite of idx 

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
            orbit_vol.name = 'Orbital Volume'
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
