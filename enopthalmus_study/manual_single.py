from const import * # Safe to use * as only CONSTANT variables
import utils
import materials

import numpy as np
# standard packages
import re
import os
import csv

import mimics

def get_centre(part):
  '''Get the geometric centre of a mimics sphere, point or spline.'''
  try: # See it if has a center already
    return part.center
  except AttributeError:
    # It didn't have a .center attribute, so not a sphere
    pass
  try: # Maybe it's a point?
    return (part.X, part.Y, part.Z)
  except AttributeError:
    # No X,Y,Z so it's not a point
    pass
  # Fall back on using the bounding box. This *should* work on anything.
  try:
    bbox = mimics.measure.get_bounding_box(part)
  except TypeError: 
    bbox = mimics.measure.get_bounding_box([part]) # Wrap a single object in a list to work on CAD objects
    
  p1 = np.array(bbox.origin)
  span = np.array(bbox.first_vector) + np.array(bbox.second_vector) + np.array(bbox.third_vector)
  return p1 + (span / 2)


def manual_fix():
    # Clean up the display, starting by turning everything off
    for o in mimics.data.objects:
      o.visible = False
    
    for m in mimics.data.masks:
      m.visible = False
    # and then show the bits we want to see
    for p in [mimics.data.spheres, mimics.data.points, mimics.data.splines]:
      for o in p:
        o.visible = True

    soft_orbits = mimics.data.masks["soft_orbits"]
    
    ## Call Edit Mask - to manually fix segmentation leaks.
    soft_clean = mimics.data.masks.duplicate(object=soft_orbits)
    soft_clean.name = "soft_clean"

    # Show the mask we want and turn on 3D preview mode
    for m in mimics.data.masks:
      m.selected = False
    soft_clean.selected = True
    mimics.view.enable_mask_3d_preview()

    soft_clean = mimics.segment.activate_edit_mask(mask = soft_clean, edit_mode = "Erase", edit_type = "Ellipse")

    eyes_left = mimics.data.spheres["left_globe"]
    eyes_right = mimics.data.spheres["right_globe"]

    ## Create the globe masks and add to the soft_clean mask
    globe_left = utils.sphere_to_mask(eyes_left)
    globe_left.name ="left_globe"
    globe_left.visible = False
    globe_right = utils.sphere_to_mask(eyes_right)
    globe_right.name ="right_globe"
    globe_right.visible = False

    ## Region Grow Volume (L) from L globe centre. 
    # Add the globe to the mask then grow from the globe centre. Need this in case the centre is in front of the anterior surface.
    temp_union = mimics.segment.boolean_operations(globe_left, soft_clean, 'Unite')
    temp_union.name = "temp_union_1"
    globe_point = get_centre(eyes_left)
    temp_left = mimics.segment.region_grow(input_mask=temp_union, target_mask=None, point=globe_point, slice_type="Axial", keep_original_mask=False, multiple_layer=True, connectivity='26-connectivity')
    temp_left.name = "temp_left"
    mimics.data.masks.delete(temp_union) # need to delete here as re-using the name below creates a new copy of the mask
    ## Subtract Globe again to leave the orbital volume
    orbit_left = mimics.segment.boolean_operations(temp_left, globe_left, 'Minus')
    orbit_left.name = "left_Orbital Volume"

    ## Repeat for Right side. 
    temp_union = mimics.segment.boolean_operations(globe_right, soft_clean, 'Unite')
    temp_union.name = "temp_union_2"
    globe_point = get_centre(eyes_right)
    temp_right = mimics.segment.region_grow(input_mask=temp_union, target_mask=None, point=globe_point, slice_type="Axial", keep_original_mask=False, multiple_layer=True, connectivity='26-connectivity')
    temp_right.name = "temp_right"
    mimics.data.masks.delete(temp_union) # delete here to be consistent with above
    orbit_right = mimics.segment.boolean_operations(temp_right, globe_right, 'Minus')
    orbit_right.name = "right_Orbital Volume"

    mimics.data.masks.delete([globe_left, globe_right, temp_left, temp_right]) # Delete a temporary masks - comment out for debugging

    ## Manually check that segmentation is good & do any needed clean-up, then continue. ui_blockng=False allows use of normal tools.
    for m in mimics.data.masks:
      m.selected = False
    orbit_left.selected = True
    orbit_right.selected = True
    mimics.view.enable_mask_3d_preview()
    
    answer = mimics.dialogs.question_box(message="Check the masks and clean up as needed, then continue.", buttons="Continue", ui_blocking=False)

    ## Threshold New Mask - Air, Fat, Muscle
    ## Boolean intersect each L Volume with each new mask 
    ## Same for Right
    ## Generate parts
    ## Calculate volume from parts and masks
    ## Write to .CSV

if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
  # Operate on the currently open project

  manual_fix()