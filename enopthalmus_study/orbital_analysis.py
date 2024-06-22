
<<<<<<< HEAD
from const import * # CONSTANT definitions for this project (* is safe as only consts)
import utils # Utility functions to simplify code
import materials # Material definitions (used to make masks)

import logging # To add messages to the Mimics log
import os

import mimics # Mimics API

DEBUG = False

def detect_eyes():
    '''Find the sphere, spline and point objects for each eye in the project'''
    num_eyes = len(mimics.data.spheres)
    num_rims = len(mimics.data.splines)
    num_pts  = len(mimics.data.points)

    if num_eyes != num_rims:
        raise ValueError(f'Aborting project - number of globes {num_eyes} does not match number of rims {num_rims}.')
        return   # ERROR - number of eye parts doesn't match

    if num_eyes > 2 or num_eyes < 1:
        raise ValueError(f'Wrong number of eyes! Expected one or two, found {num_eyes}')
        return  # ERROR - too many eyes or too few eyes, so move to next project

    # Set up a blank dict to hold the eye components
    eyes = {
        'right': {'globe': None, 'rim': None, 'point': None},
        'left':  {'globe': None, 'rim': None, 'point': None}
        } 
        
    # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
    #    X+ to the patient's Left, so -ve to the Right
    #    Y+ to the patient's Posterior, so -ve Anterior
    #    Z+ to the patient's Superior, so -ve Inferior
    # For each eye component, determine which geometry is on which side, based on the location
    # If there are two of a given component, compare them; otherwise compare to 0
    if num_eyes == 1:
        # Just need to determine the side once, and can use the globe for that
        if mimics.data.spheres[0].center[X] < 0:
            side = 'right'
        else:
            side = 'left'
        # Note that some may not have an apex point
        try:
            eyes[side]['globe'] = mimics.data.spheres[0]
            eyes[side]['rim']   = mimics.data.splines[0]
            eyes[side]['point'] = mimics.data.points[0] if len(mimics.data.points) > 0 else None
        except Exception as e:
             raise
    else:
        # num_eyes must be 2 here.
        # Need to find the side of each component indivdually, as could be entered in random order.
        # Find which index is left and which is right for each component by checking the X coordinate.

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


def make_anterior_mask(rim, globe):
  # This function is slow, as it creates a lot of geometry. 
  # Consider usinge with mimics.disabled_gui(): to turn the GUI off and speed things up.

  # Get the BoundingBox3D of the rim to calculate a mask to cover it
  bbox_rim = mimics.measure.get_bounding_box([rim]) # ([rim]) works, but (rim) doesn't
  
  delta_z = bbox_rim.third_vector[Z] # The height of the orbital rim
  spacing = delta_z / NUM_PLANES     # Distance between planes to cover the whole rim
  min_z = bbox_rim.origin[Z]         # The bottom of the orbital rim
  min_z = min_z + spacing / 2.0      # Start half a spacing above the bottom to ensure intersection   
  
  # Make a list of plane origin points, matching the globe in the X,Y plane 
  # and spaced in Z from min_z to (min_z + delta_z) (approximately max_z)
  plane_origins = [(globe.center[X], globe.center[Y], min_z + (z * spacing))
                    for z in range(NUM_PLANES)]

  # Create planes using these origins, parallel to the X,Y plane (normal is Z+)
  #norm_z = [0, 0, 1] # This is not always true, so find the correct normal from the data
  (u, v, w) = utils.mimics_image_vectors()
  norm_z = w
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) 
              for orig in plane_origins]
  
  boxes = [] # Start with a blank list for the boxes
  # For each plane, find the two intersections with the spline and make a bounding box 
  # based on those intersection points.
  for plane in z_planes:
    pt_up, pt_down = None, None
    plane_z = plane.origin[Z]

    # The spline is a continous loop. Check if each pair of points intersects this plane
    for p1, p2 in utils.looped_pairwise(rim.geometry_points):
      # If this segment of spline crosses this plane then one endpoint will be above the plane and one below. 
      # Most lines will not cross this plane, in which case both these will be False, so this is fast.
      if (p1[Z] >= plane_z and p2[Z] < plane_z): # crosses from above
        line_int = mimics.analyze.create_line(p1, p2) # temp line to get intersection point
        pt_up = mimics.analyze.create_point_as_line_and_plane_intersection(line_int, plane)
        mimics.data.lines.delete(line_int) # remove the temp line
      if (p1[Z] <= plane_z and p2[Z] > plane_z): # crosses from below
        line_int = mimics.analyze.create_line(p1, p2) # temp line to get intersection point
        pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line_int, plane)
        mimics.data.lines.delete(line_int) # remove the temp line

      # Each segment will only cross a given plane once. Stop when have found one in each direction.
      if (pt_up is not None) and (pt_down is not None):
        # Create a bounding box and add it to the list
        bb = utils.bbox_from_intersections(pt_up, pt_down, MULT_XY, spacing, SIZE_Y)
        boxes.append(bb)
        # Clean up as have used the points. Could delete here, but keep for debugging
        pt_up.visible = False
        pt_down.visible = False
=======
from const import * # CONSTANT definitions for this project (* is safe as only consts in file)
import utils # Utility functions to simplify code
import materials # Material definitions (used to make masks)

import mimics # API to access mimics (not needed inside mimics)
import mimics.analyze # API to access mimics (not needed inside mimics)

# Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
def make_crop_box(pt_up, pt_down, spacing_z):
  """Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
  Input:  two points for the posterior edge of the box and a thickness.
          The posterior edge will be extended by MULT_XY in both directions and
          the box will extend SIZE_Y anteriorly.
  Output: A mimcs.BoundingBox3d to sue for cropping a mask."""

  # Align the crop box along the XY line between the two intesection points, 
  # expanded slightly to ensure it covers the whole orbit
  vector_x = [ MULT_XY * (pt_up.x - pt_down.x), MULT_XY * (pt_up.y - pt_down.y), 0]
  vector_y = [0, SIZE_Y, 0]    # -Y is anterior, so point towards front of face
  vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
  # Put the origin above the down intesection point, midway between the z planes
  origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
  return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)

def make_orbit_mask(rim, globe):
  """Calculate Orbital contents volumes.
  Input:  a mimics.Spline delmiting the rim, and a mimics.Sphere for the globe.
  Output: a Mask from the 'surface' of the rim, extending forward, for cropping."""

  # Decompose the spline into straight lines between each pair of points.
  # Use these to intersect with planes, as mimics doesn't have a spline intersect plane function.
  # Use utils.looped_pairwise to join the last point back to the first point.
  spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) 
                  for a, b 
                  in utils.looped_pairwise(rim.points)]
  rim_geom = utils.spline_geometry(rim)
  # How tall is the orbital rim?
  delta_z = round(rim_geom.span[Z], 0) # Round the nearest whole number
  spacing = delta_z/NUM_PLANES
  
  # Make a list of plane origin points, matching the globe in the X,Y plane 
  # and spaced in Z from min_z to (min_z + delta_z) (approximately max_z)
  globe_x = globe.center[X]
  globe_y = globe.center[Y]
  min_z = rim_geom.min[Z]
 
  plane_origins = [(globe_x, globe_y, min_z + (z + 1) * spacing)
                    for z
                    in range(NUM_PLANES - 1)]

  # Create planes using these origins paralle to the X,Y plane (normal is Z+)
  norm_z = [0, 0, 1]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) 
              for orig 
              in plane_origins]
  
  # For each plane, find intersections with the spline lines 
  # and create a bounding box based on the intersection points.
  boxes = [] # Start with an empty list for the boxes
  for plane in z_planes:
    pt_up, pt_down = None, None
    for line in spline_lines:
      # If a line crosses this plane it will have one endpoint above the plane and one below. 
      # Most lines will not cross this plane, in which case both these will be False
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z]) # TODO: should this be >=
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z]) # TODO: should this be =<
      if (from_above):
        pt_up = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)
      if (from_below):
        pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)
    
      # Each line will only cross any given plane once, so can stop
      # when there is an intersection in each direction.
      if (pt_up is not None) and (pt_down is not None):
        boxes.extend(make_crop_box(pt_up, pt_down)) # add the current crop box to the list
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
        break # We have found that line for this plane   
    else:
      # Fell through the loop without breaking
      print(f"did not find intersection for plane {plane}")
<<<<<<< HEAD
  
  # Thicken the first and last bounding boxes to ensure full overlap.
  boxes[0]  = utils.bbox_thicken(boxes[0], -EXTRA_Z) # first goes down
  boxes[-1] = utils.bbox_thicken(boxes[-1], EXTRA_Z) # last goes up
    
  # Create a series of masks, covering all materials, clipped by each bounding box
  # As each mask is created, unite it into a combined mask to cut away tissue anterior to the rim.
  mask_ant = mimics.segment.create_mask(select_new_mask=False)
  for b in boxes:
      mask_ab = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                        threshold_min=materials.MIN_GV,
                                        threshold_max=materials.MAX_GV, 
                                        bounding_box= b)
      old_ant = mask_ant # unite creates a new mask, so save the current one to delete after.
      mask_ant = utils.masks_unite(mask_ant, mask_ab)
      mimics.data.masks.delete(mask_ab)
      mimics.data.masks.delete(old_ant)
  
  return mask_ant

# start_off = time.perf_counter()
# with mimics.disabled_gui_update(): 
#   mask_anterior = make_anterior_mask(rim, globe)

# time_off = time.perf_counter() - start_off
# print(f"GUI on took {time_off}")

def make_orbit_ROI(rim, point = None):
  # Make a BoundingBox3D based on the orbital rim; a bit bigger in X & Z and much deeper in Y
  # If there is an orbital apex point use that, otherwise use a standard amount to extend Y.
  
  # Get the BoundingBox3D of the rim to work out where to stat the ROI
  if point is None:
    bbox_rim = mimics.measure.get_bounding_box((rim, )) # need to give ist a list (or tuple)
    extra_y = POST_EXPAND
  else:
    bbox_rim = mimics.measure.get_bounding_box((rim, point)) # Bounding box of all objects in list
    extra_y = BOX_EXPAND[Y]
=======

  # Adjust the first and last bounding box to ensure full overlap.
  boxes[0].third_vector[2] = -MULT_Z * boxes[0].third_vector[2]  # first goes down
  boxes[-1].third_vector[2] = MULT_Z * boxes[-1].third_vector[2] # and last goes up

  # Now use the bounding boxes to create a Region Of Interest
  mask_union = mimics.segment.create_mask() # Start with an empty mask
  # Loop through each box, create a mask, crop it with the box, and Union it together.
  for bbox in boxes:
    # Create a temporary mask, thresholding it to cover everything and
    # crop it with BoundingBox for this plane.
    m = mimics.segment.threshold(mask = mimics.segment.create_mask(), 
                                 threshold_min = materials.MIN_GV, 
                                 threshold_max = materials.MAX_GV,
                                 bounding_box = bbox)
    # Union with the existing mask
    mask_union = utils.masks_unite(mask_union, m)
    # then delete it
    mimics.data.masks.delete(m)

  # The union of all the masks is now everything in front of the orbital rim
  return(mask_union)

def make_orbit_ROI(rim):
  # Get the rim extents 
  rim_geometry = utils.spline_geometry(rim)
  min_pt = rim_geometry.min_point
  max_pt = rim_geometry.max_point
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734

  # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
  #    X+ to the patient's Left
  #    Y+ to the patient's Posterior
  #    Z+ to the patient's Superior
<<<<<<< HEAD
  # Add BOX_EXPAND to both sides in all directions, then a bit more for POST_EXPAND
  # Move the origin back by BOX_EXPAND 
  new_origin = [(a - b) for a, b in zip(bbox_rim.origin, BOX_EXPAND)]
  # Make the extents bigger by 2 * expand as we have moved the orgin back by 1 * expand
  extents = (bbox_rim.first_vector[X], bbox_rim.second_vector[Y], bbox_rim.third_vector[Z])
  vector_x = (extents[X] + 2 * BOX_EXPAND[X], 0, 0)
  vector_z = (0, 0, extents[Z] + 2 * BOX_EXPAND[Z])
  vector_y = (0, extents[Y] + BOX_EXPAND[Y] + extra_y, 0) # extra for POST
  # A perhaps more pythonic but harder to understand way?
  # new_extents = [(a + 2 * b) for a, b in zip(extents, BOX_EXPAND)]
  # vector_x = (new_extents[X], 0, 0)
  # vector_y = (0, new_extents[Y] + POST_EXPAND, 0)
  # vector_z = (0, 0, new_extents[Z])
  orbit_ROI = mimics.BoundingBox3d(new_origin, vector_x, vector_y, vector_z)

  return orbit_ROI

from utils import Events

def segment_orbit_orig(rim, globe, point, side):

    t = Events()

    # Make the anterior blocking mask and the Orbital ROI
    with mimics.disabled_gui_update():
        m_anterior = make_anterior_mask(rim, globe) # functoon in orbital_analysis
        if DEBUG: m_anterior.name = 'm_anterior'
        # Make a bigger one as well to plug gaps
        m_ant_big = mimics.segment.morphology_operations(m_anterior, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_ant_big', limited_to_mask=None)
        if DEBUG: m_ant_big.name = 'm_ant_big'
        # Hide all the planes as they get in the way of seeing the rim etc.
        for p in mimics.data.planes:
            p.visible = False

    t.add('make anterior mask')

    orbit_ROI = make_orbit_ROI(rim, point) # functoon in orbital_analysis
    t.add('make orbit ROI')

    # Create the bone mask, limited to the orbit_ROI
    #m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(148), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    # Make the bone lower threshold split the difference between Bone (CT) [226 : 3071] and Spongial BOne (Adult CT) [148 : 661]?
    #m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(187), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)

    # A different approach is to use the bone thresholds and then use a local threshold to add lower density material.
    m_bone_1 = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    if DEBUG: m_bone_1.name = 'm_bone_1'
    m_bone = mimics.segment.local_threshold(mask=m_bone_1, threshold_min=mimics.segment.HU2GV(148), threshold_max=mimics.segment.HU2GV(3071), search_distance=2, isolate=True)
    if DEBUG: m_bone.name = 'm_bone'
    # Delete the initial bone mask (mimics doesn't reuse objects passed as arguments)
    if not(DEBUG): mimics.data.masks.delete(m_bone_1)

    # Could use keep largest (on uncropped mask), but that runs the risk of losing disconnected bone fragments.
    # Instead, clean up by removing small areas with 'Open' 1 x 8
    #m_bone = mimics.segment.morphology_operations(m_bone, operation='Open', number_of_pixels=1, connectivity=8, target_mask_name='m_bone', limited_to_mask=None)

    # Better appears to be to use smooth_mask
    m_bone = mimics.segment.smooth_mask(m_bone)
   
    # Don't appear to use this object?
    # p_bone = mimics.segment.calculate_part(mask=m_bone, quality='High')
    # p_bone.name = 'p_bone'

    # Create the air mask, limited to the orbit_ROI
    # Use new 'Air' threshold, based on mimics skin = [-718 : -177] (air was [-1024 : -200])
    m_air = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200), bounding_box=orbit_ROI)
    if DEBUG: m_air.name = 'm_air'

    t.add('thresholding')

    # Use erode to seperate the parts
    m_air_erode = mimics.segment.morphology_operations(m_air, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name='m_air_erode', limited_to_mask=None)
    if DEBUG: m_air_erode.name = 'm_air_erode'

    # Want the air outside of the face. Should be the largest, and could make bbox longer anterior to make that more likely.
    # m_air_external = mimics.segment.keep_largest(mimics.data.masks[-1])

    # Or, could region_grow using the front of the bbox to set a seed point.
    # The most antero-lateral point on left eye is away from origin, on the right it is the origin. 
    # The point is offset a small amount (defaults delta = 2 mm) into the box to make sure it is inside the mask.
    lateral_pt = utils.antero_lateral(mimics.measure.get_bounding_box(m_air_erode), side)
    m_air_external_1 = mimics.segment.region_grow(input_mask=m_air_erode, target_mask=None, point=lateral_pt,
                                                slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    if DEBUG: m_air_external_1.name = 'm_air_external_1'
    # Dilate back to previous size, but should only have the external part of the air now
    m_air_external = mimics.segment.morphology_operations(m_air_external_1, operation='Dilate', number_of_pixels=10, connectivity=26, target_mask_name='m_air_external', limited_to_mask=m_air)
    if DEBUG: m_air_external.name = 'm_air_external'

    # Subtract the external from all the air in the ROI, to leave the internal air
    m_air_internal = mimics.segment.boolean_operations(m_air, m_air_external, 'Minus')
    if DEBUG: m_air_internal.name = 'm_air_internal'

    t.add('split air mask')

    # Open internal air to remove small regions caused by low density tissue overlapping with our definition of air.
    m_air_int_open = mimics.segment.morphology_operations(m_air_internal, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_open', limited_to_mask=None)
    if DEBUG: m_air_int_open.name = 'm_air_int_open'

    # Dilate internal air 2 x 8 and unite with bone (maybe 1 x 8 would be better?)
    m_air_int_dilate = mimics.segment.morphology_operations(m_air_int_open, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_dilate', limited_to_mask=None)
    if DEBUG: m_air_int_dilate.name = 'm_air_int_dilate'

    m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
    m_unite_bone_air.name = 'm_bone+air'

    # Clean up the list of maskes used to create the 'Opened' air mask
    if not(DEBUG): mimics.data.masks.delete((m_bone, m_air, m_air_erode, m_air_external_1, m_air_external, m_air_internal, m_air_int_open, m_air_int_dilate))

    t.add('dilate and unite')

    # Sadly, local_threshold doesn't seem to work when called from the API, although it does using the GUI.
    # Use Local Threshold to pick up thin/low density bone close to mask
    # m_local_thresh = mimics.segment.local_threshold(mask=m_unite_bone_air, 
    #                                       threshold_min=mimics.segment.HU2GV(148), 
    #                                       threshold_max=mimics.segment.HU2GV(3071), search_distance=2, isolate=True)
    # Work around that with the following kludge:
    m_fake_local_1 = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                            threshold_min=mimics.segment.HU2GV(148), threshold_max=mimics.segment.HU2GV(3071), 
                                            bounding_box=mimics.measure.get_bounding_box((m_unite_bone_air))
                                            )
    if DEBUG: m_fake_local_1.name = 'm_fake_local_1'

    m_fake_local = mimics.segment.boolean_operations(mask_a=m_unite_bone_air, mask_b=m_fake_local_1, operation='Unite')
    if DEBUG: m_fake_local.name = 'm_fake_local'
    #m_fake_local = mimics.segment.smooth_mask(m_fake_local)
    #p_fake_local = imics.segment.calculate_part(m_fake_local, quality = 'High')
    #p_fake_smoothed = mimics.tools.smooth(p_fake_local, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=False)
    #m_wrap_last = mimics.tools.wrap(object_to_wrap= p_fake_smoothed, smallest_detail=0.2, gap_closing_distance=5, dilate_result=False, protect_thin_walls=True, keep_originals=True)
    ## Clean up temp objects
    #mimics.data.objects.delete(m_fake_local, p_fake_local, p_fake_smoothed)

    # Smartfill the 'orbital boundary' mask to close small holes
    #m_fill = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=7) # Is this distance too big?
    m_fill = mimics.segment.smart_fill_global(mask=m_fake_local, hole_closing_distance=3) # Is this distance too big?
    m_fill = mimics.segment.smooth_mask(m_fill)
    if DEBUG: m_fill.name = 'm_fill'

    t.add('smart fill')

    # Create a part from the filled mask, smooth it and wrap it to close more gaps in orbital boundary, then convert back to a mask
    p_fill = mimics.segment.calculate_part(mask=m_fill, quality='High')
    if DEBUG: p_fill.name = 'p_fill'

    t.add('calculate part')

    # Note keep_originals=True for debugging - not used later so could be deleted here
    p_smooth = mimics.tools.smooth(p_fill, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
    if DEBUG: p_smooth.name = 'p_smooth'

    t.add('smooth')

    p_wrap = mimics.tools.wrap(p_smooth, smallest_detail=0.2, gap_closing_distance=10,
                                dilate_result=False, protect_thin_walls=True, keep_originals=True)
    p_wrap.name = 'p_wrap'

    t.add('wrap')

    ## Clean up temp objects
    if not(DEBUG): mimics.data.objects.delete((m_fake_local_1, m_fake_local, m_fill, p_fill, p_smooth))

    # Make the smoothed, wrapped part back into a mask. It should now define the orbit boundary
    m_wrap = mimics.segment.calculate_mask_from_part(part=p_wrap, target_mask=None)
    m_wrap.name = 'm_wrap'

    t.add('mask from part')

    # Remove tissue in front of the orbit to isolate the contents. This should be just subtracting the 'anterior' mask,
    # but due to the resolution of the blocks there can still be gaps between the bone and the anterior mask.

    # Make two masks - one expanded to segment the contents and one to size to expand back into
    # Add the anterior block to the subtract mask, to remove tissue in front of the orbital rim
    m_subtract = mimics.segment.boolean_operations(m_wrap, m_anterior, 'Unite')
    if DEBUG: m_subtract.name = 'm_subtract'
    # Now make the bigger one
    m_oversize = mimics.segment.morphology_operations(m_subtract, operation='Close', number_of_pixels=2, connectivity=26, target_mask_name='m_oversize', limited_to_mask=None)
    if DEBUG: m_oversize.name = 'm_oversize'
    
    # Create a block of everything (should be soft tissue + air?)
    m_temp_orbit = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    if DEBUG: m_temp_orbit.name = 'm_temp_orbit'

    t.add('make temp orbit')

    # Subtract the to-size orbital boundary from above, leaving the orbital contents + extra tissue areas
    m_vol_to_fill = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract, operation='Minus')
    if DEBUG: m_vol_to_fill.name = 'm_vol_to_fill'

    # Subtract the oversize orbital boundary from above, hopefully seperating the orbit from anything else
    m_segmented_1 = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_oversize, operation='Minus')
    if DEBUG: m_segmented_1.name = 'm_segmented_1'
    # then shrink what's left to get seperate regions
    m_segmented_2 = mimics.segment.morphology_operations(m_segmented_1, operation='Open', number_of_pixels=6, connectivity=8, target_mask_name='m_segmented', limited_to_mask=None)
    if DEBUG: m_segmented_2.name = 'm_segmented_2'
    # Delete temp maskS
    if not(DEBUG): mimics.data.masks.delete((m_subtract, m_oversize, m_temp_orbit, m_segmented_1))
    
    t.add('m_subtract minus temp_orbit')

    #mimics.logging.log_system_message(level=logging.DEBUG, message="making seed point")
    mimics.logging.log_system_message(level=logging.INFO, message="making seed point")
    #mimics.logging.log_user_message(level=logging.DEBUG, message="making seed point")
    mimics.logging.log_user_message(level=logging.INFO, message="making seed point")

    # THIS INITIALLY FAILED. ####################################################
    # seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])
    
    # # Region grow the soft tissue - boundary from the back of the globe, to seperate intra-orbital component
    # m_core = mimics.segment.region_grow(input_mask=m_segmented, target_mask=None, point=seed,
    #                                             slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    # mimics.data.masks[-1].name = 'm_core'

    # Convert globe (which is a Sphere) to a mask
    m_globe = utils.sphere_to_mask(globe)
    # Add to the m_segmented mask
    m_segmented = mimics.segment.boolean_operations(mask_a=m_segmented_2, mask_b=m_globe, operation='Unite')
    if DEBUG: m_segmented.name = 'm_segmented'
        # Smart fill this combined mask to close any interior holes
    m_segmented_fill = mimics.segment.smart_fill_global(mask=m_segmented, hole_closing_distance=7) # Is this distance too big?
    if DEBUG: m_segmented_fill.name = 'm_segmented_fill'

    # Delete temp masks
    if not(DEBUG): mimics.data.masks.delete((m_segmented_2, m_segmented))

    # Trying to select just the interior contents by region growing from the globe fails, even using
    # a point at the back of the globe. Can try to fix that by adding the globe to the volume mask
    # and then region growing, so will pick up via any point of contact. Then subtract globe again.

    # THIS ALSO FAILED
    #seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])
    
    
    # Try instead to add the globe to the orbital contents, then set a seed at the centre of the globe, which should always succeed.
    seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])
    # Region grow the soft tissue - boundary from the globe, to seperate intra-orbital component
    m_core_1 = mimics.segment.region_grow(input_mask=m_segmented_fill, target_mask=None, point=seed,
                                                slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    if DEBUG: m_core_1.name = 'm_core_1'
    t.add('region grow core')

    # Subtract the globe again
    m_core = mimics.segment.boolean_operations(mask_a=m_core_1, mask_b=m_globe, operation='Minus')
    if DEBUG: m_core.name = 'm_core'

    # Delete temp masks
    if not(DEBUG): mimics.data.masks.delete((m_core_1, m_segmented_fill, m_globe))
    
    # Dilate the core region back out to fill the volume, limited to the original volume to fill mask
    m_orbit_1 = mimics.segment.morphology_operations(m_core, operation='Dilate', number_of_pixels=10, connectivity=8, target_mask_name='m_orbit', limited_to_mask=m_vol_to_fill)
    if DEBUG: m_orbit_1.name = 'm_orbit_1'
    # This will extrude anteriorly if there are gaps between the bone & anterior mask, so subtract the elarged anterior mask
    m_orbit_2 = mimics.segment.boolean_operations(mask_a=m_orbit_1, mask_b=m_ant_big, operation='Minus')
    if DEBUG: m_orbit_2.name = 'm_orbit_2'
    # Subtract the wrapped bone again for good measure
    m_orbit_3 = mimics.segment.boolean_operations(mask_a=m_orbit_2, mask_b=m_wrap, operation='Minus')
    if DEBUG: m_orbit_3.name = 'm_orbit_3'
    
    # There may be bits outside fo the orbit (not sure how), so region grow to get onlu connected parts, then smooth the mask
    m_orbit = mimics.segment.region_grow(input_mask=m_orbit_3, target_mask=None, point=seed,
                                                slice_type='Axial', keep_original_mask=False, multiple_layer=True, connectivity='6-connectivity')
    m_orbit = mimics.segment.smooth_mask(m_orbit)
    if DEBUG: m_orbit,name = 'm_orbit'

    # Clean up any temp or unneeded masks
    # It looks like mimics adds new masks to the list every time one is returned,
    # which is to say if a mask is made by m_a = fn(m_a) we get two copies of m_a.
    # To clean up the masks list need to use temoprary names, or maybe just delete everything
    # that is not on the list to keep?
    if not(DEBUG): mimics.data.masks.delete((m_core, m_orbit_1, m_orbit_2, m_orbit_3))

    print(f'Total time taken {t.elapsed()} seconds')
    print(t.as_dict())

    return m_orbit, orbit_ROI

def segment_orbit(rim, globe, point, side):

    t = Events()

    s_id = side[0] # Just the first letter to use for identifying masks etc.
    s_id = side # Or in fact use the whole word

    # Make the anterior blocking mask and the Orbital ROI
    with mimics.disabled_gui_update():
        m_anterior = make_anterior_mask(rim, globe) # functoon in orbital_analysis
        m_anterior.name = f'{s_id}_m_anterior'
        # Make a bigger one as well to plug gaps
        m_ant_big = mimics.segment.morphology_operations(m_anterior, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_ant_big', limited_to_mask=None)
        m_ant_big.name = f'{s_id}_m_ant_big'
        # Hide all the planes as they get in the way of seeing the rim etc. (Should probably delete them since this is working)
        for p in mimics.data.planes:
            p.visible = False

    t.add('make anterior mask')

    orbit_ROI = make_orbit_ROI(rim, point) # functoon in orbital_analysis
    t.add('make orbit ROI')

    # Create the bone mask, limited to the orbit_ROI
    # Make the bone lower threshold split the difference between Bone (CT) [226 : 3071] and Spongial BOne (Adult CT) [148 : 661]?
    m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(187), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    m_bone.name = f'{s_id}_m_bone'
    m_bone = mimics.segment.smooth_mask(m_bone)

    # Get air by thresholding, erode to seperate into internal and external regions
    m_air = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200), bounding_box=orbit_ROI)
    m_air.name = f'{s_id}_m_air'
    m_air_erode = mimics.segment.morphology_operations(m_air, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name=f'{s_id}_m_air_erode', limited_to_mask=None)

    # Want the air outside of the face, so region_grow using the front of the bbox to set a seed point.
    # The most antero-lateral point on left eye is away from origin, on the right it is towards the origin. 
    # The point is offset a small amount (defaults delta = 2 mm) into the box to make sure it is inside the mask.
    lateral_pt = utils.antero_lateral(mimics.measure.get_bounding_box(m_air_erode), side)
    # Files 11 & 12 failed at line 424
    # Get a point on the antero-lateral corner of the air mask to start region growing in the external part.
    m_air_ext_1 = mimics.segment.region_grow(input_mask=m_air_erode, target_mask=None, point=lateral_pt,
                                                    slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')

    # Dilate only the external part back to original size, then subtract that from the total to get the internal part
    m_air_external = mimics.segment.morphology_operations(m_air_ext_1, operation='Dilate', number_of_pixels=10, connectivity=26, target_mask_name=f'{s_id}_m_air_external', limited_to_mask=m_air)
    # Subtract external from total air to get internal
    m_air_internal = mimics.segment.boolean_operations(m_air, m_air_external, 'Minus')
    m_air_internal.name = f'{s_id}_m_air_internal'

    t.add('split air mask')

    # Open internal air to remove small regions caused by low density tissue overlapping with our definition of air.
    m_air_int_open = mimics.segment.morphology_operations(m_air_internal, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name=f'{s_id}_m_air_int_open', limited_to_mask=None)
    # Dilate what remains to extend to cover thin bone regions, then unite with bone
    m_air_int_dilate = mimics.segment.morphology_operations(m_air_int_open, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=f'{s_id}_m_air_int_dilate', limited_to_mask=None)
    m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
    m_unite_bone_air.name = f'{s_id}_m_bone+air'

    # Clean up temporary masks. Need them as a tuple or list to give to delete(), hence the two step process. 
    #temp_masks = (m_air, m_air_erode, m_air_ext_1, m_air_int_open, m_air_int_dilate)
    #mimics.data.masks.delete(temp_masks) # Could just call delete with a second set of brackets & list the names there

    t.add('dilate and unite')

    # Plan was to use local threshold to add more bone around the edges of the selected mask, but sadly, 
    # local_threshold doesn't seem to work when called from the API, although it does using the GUI.
    # Work around that by just using the lower threshold in the first place. This may get extra non-bone soft tissue, but that's hard to fix.

    # Fill holes in the bone+air mask, then smooth it. Use keep_largest as some scans have air pockets in front of globe that breaks wrap later.                                         
    m_fill_ba = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=3) # Is this distance too big?
    m_fill_ba = mimics.segment.smooth_mask(m_fill_ba)
    m_fill_ba = mimics.segment.keep_largest(m_fill_ba)
    m_fill_ba.name = f'{s_id}_m_fill_ba'

    t.add('smart fill')

    # Next turn the mask into a part, smooth and wrap it to get a hopefully watertight orbital border all around.
    p_fill_ba = mimics.segment.calculate_part(mask=m_fill_ba, quality='High')
    p_fill_ba.name = f'{s_id}_p_fill_ba'
    # Note keep_originals=True for debugging - not used later so could be deleted here
    p_smooth_ba = mimics.tools.smooth(p_fill_ba, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
    p_smooth_ba.name = f'{s_id}_p_smooth_ba'
    p_wrap_ba = mimics.tools.wrap(p_smooth_ba, smallest_detail=0.2, gap_closing_distance=6, dilate_result=False, protect_thin_walls=True, keep_originals=True)
    p_wrap_ba.name = f'{s_id}_p_wrap_ba'

    # Make the smoothed, wrapped part back into a mask. It should now define the orbit boundary
    m_wrap_ba = mimics.segment.calculate_mask_from_part(part=p_wrap_ba, target_mask=None)
    m_wrap_ba.name = f'{s_id}_m_wrap_ba'

    t.add('mask from part')

    temp_start = len(mimics.data.masks) - 1 # Note start of temporary masks to delete
    
    # Convert globe (which is a Sphere) to a mask
    m_globe = utils.sphere_to_mask(globe)
    
    # Make a temporary 'all tissues' mask and subtract each boundary mask in turn
    t1 = mimics.segment.threshold(mimics.segment.create_mask(), threshold_min=mimics.segment.HU2GV(-1024), threshold_max=mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    t2 = mimics.segment.boolean_operations(t1, m_wrap_ba, operation='Minus')
    t2 = mimics.segment.boolean_operations(t2, m_fill_ba, operation='Minus')
    t2 = mimics.segment.boolean_operations(t2, m_air_external, operation='Minus')
    t2 = mimics.segment.boolean_operations(t2, m_ant_big, operation='Minus')
    
    # Temporarily add in the globe, so that the globe centre point will alwaye be in the mask
    t3 = mimics.segment.boolean_operations(t2, m_globe, operation='Unite')
    # Region grow from the centre of the globe, then subtract the globe again
    t4 = mimics.segment.region_grow(input_mask=t3, target_mask=None, point=globe.center, slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    t5 = mimics.segment.boolean_operations(t4, m_globe, operation='Minus')
    # The resulting mask may be a bit smaller than the bony orbit, so expand it a bit, then remove the filled bony orbit
    t6 = mimics.segment.morphology_operations(input_mask=t5, operation='Dilate', number_of_pixels=2, target_mask_name='t6', connectivity=8,limited_to_mask=None)
    t7 = mimics.segment.boolean_operations(t6, m_fill_ba, operation='Minus')
    
    # Finally subtract the globe, as the dilate will have expanded into that region as well
    m_orbit = mimics.segment.boolean_operations(t7, m_globe, operation='Minus')
    # And smooth it for good measure
    m_orbit = mimics.segment.smooth_mask(m_orbit) # file 17 failed here
    m_orbit.name = f'{s_id}_m_orbit'

    # Clean up temporary masks.
    temp_masks = [mimics.data.masks[i] for i in range(temp_start, len(mimics.data.masks)-1)] # All masks from temp_start to just before last mask
    mimics.data.masks.delete(temp_masks) # Could just call delete on th elist comprehension, or a list of named masks.

    print(f'Total time taken {t.elapsed()} seconds')
    print(t.as_dict())

    return m_orbit, orbit_ROI
=======

  # Make a bounding box that extends beyond the rim by +/- 10 X and +/- 15 Z
  # with Y extending ant by -5 and posterior 80. Could go to max_pt[Y] + 80 ?
  expand = (10, 5, 15) # I think Y could be more?

  box_orig = (min_pt[X] - expand[X], min_pt[Y] - expand[Y], min_pt[Z] - expand[Z])
  vector_x = (max_pt[X] - min_pt[X] + (2 * expand[X]), 0, 0)
  vector_y = (0, 80, 0)
  vector_z = (0, 0, max_pt[Z] - min_pt[Z] + (2 * expand[Z]))
  
  orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
  return orbit_ROI
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
