# Measure the intersection volumes for various tissues, 
# given folders of projects with already segmented orbits.
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching

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

def write_results(study_info, input_info, volumes, results_file):
  # Collapse the study info to a single dict for logging, as don't care whether they are image or subject info
  collapsed_study = {k: v for d in study_info.values() for k, v in d.items()}
  # Collapse the input information into a single dict for logging, as side is encoded in the label
  collapsed_inputs = {k: v for d in input_info.values() for k, v in d.items()}
  # Collapse the volume data into a single dict for logging, encoding side (from the dict) in the label
  collapsed_volumes = {side + "|" + k: v for side, d in volumes.items() for k, v in d.items()}
  
  # Combine the two sets for logging
  combined_for_log = {**collapsed_study, **collapsed_inputs, **collapsed_volumes}
  # The header will be all the keys in this dict, and the data will be all the values
  headers = list(combined_for_log.keys())
  results = list(combined_for_log.values())

  # On first call this will create the file and write the headers, then the results.
  # Subsequent calls wil only write the results.
  log_to_file(results_file, headers, results)

def snapshot_3D(objects, file_name):
  '''Create a snapshot to the 3D window showing only the objects listed and write to file_name.'''
  # Show only what was requested
  # Collect the current visibility of everything
  state = [o.visible for o in mimics.data.objects]
  # Then hide them all
  for o in mimics.data.objects: o.visible = False
  # Show the requested objects
  try:
    for o in objects: 
        o.visible = True # Make sure all the objects are shown
        o.selected = True # Need to select masks for them to show up in the 3D view
  except TypeError:
        objects.visible = True # if we only have one object then it's not iterable
        objects.selected = True
        
  picture_bb = mimics.measure.get_bounding_box(objects=objects)
  # Zoom each view to cover all the given objects
  for v in mimics.data.views: 
      view_cam = mimics.view.get_camera(v)
      view_settings = view_cam.get_settings()
      view_settings.zoom_to_bounding_box(picture_bb, zoom_factor=0.8)
      view_cam.set_settings(view_settings)
  
  # Use the 3D view
  mimics.view.enable_mask_3d_preview() # MAke sure preview is on so can see the masks
  settings = mimics.view.get_camera(view = mimics.data.views['3D']).get_settings()
  settings.zoom_to_bounding_box(picture_bb, zoom_factor=1)
  mimics.file.export_view_by_type(filename=file_name, view='3D', image_type = 'jpg', camera_settings=settings)

  # restore the saved visiblity
  for i, v in enumerate(state): mimics.data.objects[i].visible = v

def things_to_see():
  '''Create a list of objects we want to show in the snapshot'''
  things  = mimics.data.masks.filter('(?i)(left|right)_Orbital Volume$', regex=True)
  things += mimics.data.splines.filter("rim$", regex=True)
  things += mimics.data.spheres.filter("globe$", regex=True)
  things += mimics.data.points.filter("left|right", regex=True)
  return things

def make_project_neat():
    '''Hide all but relevant objects'''
    # First hide everything
    for o in mimics.data.objects: o.visible = False
    # Now show the thing we want
    for o in things_to_see():
        o.visible = True
    # Add some others
    for o in mimics.data.masks.filter('(?i)(left|right)_(muscle|fat|air)$', regex=True):
        o.visible = True
    
def find_eyes():
  '''Return a dict with the globe, rim and (optional) point for each eye, labelled by side.'''
  num_eyes = len(mimics.data.spheres)
  num_rims = len(mimics.data.splines)
  num_pts  = len(mimics.data.points)

  if num_eyes != num_rims:
    print(f'Number of globes {num_eyes} does not match number of rims {num_rims}.')
    raise IndexError

  if num_eyes > 2 or num_eyes < 1:
    print(f'Wrong number of eyes! Expected one or two, found {num_eyes}')
    raise IndexError

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
  
  return eyes

def measure_project():
  '''Process a open mimics project, analysing as many eyes as exist within it.'''

  study_info, dicom_tags = gather_project_info()
 
  try:
    eyes = find_eyes()
            
    # Create dicts to hold measured volumes & inputs (spline bounds, globe & point location)
    # for each eye
    volumes = {}
    input_info = {}
    
    # Analyze each eye in the project
    for side, eye_parts in eyes.items():

        print(f'********** process {side}')

        rim = eye_parts['rim']
        globe = eye_parts['globe']
        point = eye_parts['point']
        side_label = side + '|'  # add a seperator to enable splitting out the side later

        # Find the Orbital Volume mask for this side. 
        # End the regex with '$' to skip any experimental or trial masks (e.g. with/without sinus)
        m_orbit_vol = mimics.data.masks.find(f'(?i){side}_Orbital Volume$', regex=True) # Use perl style ignore case flag

        if m_orbit_vol is None:
            # couldn't find the mask, so raise an insex error & bail
            raise (IndexError, ValueError)
            # Huston, we have a problem. Bail without returning results
            return

        # NEW - as we have changed the image set for the volume masks need to make sure that the image set they are linked to is active 
        mimics.data.images.set_active(m_orbit_vol.image)

        # Convert globe (which is a Sphere) to a mask
        m_globe = utils.sphere_to_mask(globe)
        m_globe.name = f'{side}_m_globe'

        # NEW - write the subtracted masks back to the project so that everything is in the same state
        # and subtract it from the orbital volume (already done for many).
        m_intersect_vol = mimics.segment.boolean_operations(m_orbit_vol, m_globe, 'Minus')
        # Put the new subtracted mask in place of the existing orbit mask
        m_intersect_vol.name = m_orbit_vol.name
        m_orbit_vol.name = f'{m_orbit_vol.name}_orig'
        m_orbit_vol = m_intersect_vol

        m_orbit_vol.name = f'{side}_m_intersect_vol'
        check_volume = m_orbit_vol.volume

        # and make into a Part
        p_orbit_vol = mimics.segment.calculate_part(m_orbit_vol, quality='High')
        p_orbit_vol.name = f'{side}_p_intersect_vol'

        orbit_vol_ROI = mimics.measure.get_bounding_box(m_orbit_vol)  # Material masks to insersect only need to be this bigs

        # Create a list of masks and corresponding parts for each material
        # in the orbit_materials dict. Put the orbital volume first as it should always exist.
        masks = {'orbital': m_orbit_vol}
        parts = {'orbital': p_orbit_vol}
        for matl in orbit_materials:
            # Mask of where this material overlaps with intersect_vol_mask
            masks[matl] = utils.mask_from_material('m_' + side + '_' + matl, orbit_materials[matl], bounding_box=orbit_vol_ROI)
            masks[matl] = utils.masks_intersect(masks[matl], m_orbit_vol)
            masks[matl].name = f'{side}_m_{matl}'

        # Can get volume directly from mask, different to parts, so use both methods to compare.
        # Initialise the vols dict using the mask & part that should always have a volume.
        # This dict has two items: a dict of volumes from the masks and one from the parts (where they exist)
        vols = {'mask': {'orbital': m_orbit_vol.volume},
                'part': {'orbital': p_orbit_vol.volume}}
        # Can't just make everything a part as some masks may be empty, so catch that.
        for name, mask in masks.items():
            vols['mask'][name] = mask.volume
            if mask.number_of_pixels == 0:
                parts[name] = None
                vols['part'][name] = 0
            else:
                parts[name] = utils.part_from_mask(side_label + name, mask)
                vols['part'][name] = parts[name].volume

        # Record the results for this eye, adding the source (mask or part) to each material label
        volumes[side] = {source + '_' + k: v for source,
                         d in vols.items() for k, v in d.items()}

        # Record the input info for this eye
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
    # Having processed as many eyes as exist, return the results for logging
    return study_info, input_info, volumes

  except (IndexError, ValueError):
    # Huston, we have a problem. Bail without returning results
    return

##############################################################################
if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  # This version has one folder with the segmenting person int he file name
  
  root = r'D:\Projects & Research\Enophthalmos Study\re-do_DICOM'
  # Put the combined results in the root, rather than one file per user
  results_file = Path(os.path.join(root, 'results.csv'))

  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]

  num_projects = len(projects)

  for i, p in enumerate(projects): 
    try:
        user = re.match(pattern=r'.*\.(\w+)\.mcs$', string=p)[1]

        print(f'********** project {i+1} of {num_projects} \t user {user} filename {p}')
        mimics.file.open_project(filename=p, read_only_mode=True)
        # Process the current project file and return the results for any eyes it contains.
        # This resupposes that a f'{side}_Orbital Volume' mask exists for each eye to be measured.
        study_info, input_info, volumes = measure_project() 

        # Save a snapshot
        snapshot_3D(things_to_see(), p.replace('.mcs', '.snapshot.jpg'))

        make_project_neat() # Show only relevant objects

        # This version saves changes. Some projects already have the globes subtracted, but
        # want them all top be the same for easier comparison, so save the final version.
        mimics.file.close_project()

        # Add the user name to the study_info. Needs to be a dict to match rest of study_info structure.
        study_info['analysis'] = {'user ': user}
        # Having processed as many eyes as exist, write the results
        write_results(study_info, input_info, volumes, results_file)
        
    except (IndexError, ValueError):
        mimics.file.close_project()  # close the currently open project
        # Move to the next project
        continue

  print('********** Finished all files **********')