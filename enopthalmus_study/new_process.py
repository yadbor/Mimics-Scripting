# New insights - version based on the bounding box of all objects

from const import * # Safe to use * as only CONSTANT variables
import utils
import materials

import numpy as np
# standard packages
import re
import os
import csv

import mimics

# Write measurements etc. to a .CSV file.
from result_logger import Path, log_to_file

import time

DEFAULT_BASIS = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

def v_hat(v):
    '''Return a normalised unit vector.''' 
    mag = np.linalg.norm(v)
    return([i/mag for i in v])

def mimics_basis_vectors(img):
    '''return the three unit vectors that descrivbe this image volume.'''
    p0 = img.get_voxel_center([0, 0, 0])
    d = img.get_voxel_buffer().shape
    x = img.get_voxel_center([d[0]-1, 0, 0])
    y = img.get_voxel_center([0, d[1]-1, 0])
    z = img.get_voxel_center([0, 0, d[2]-1])
    span = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
    (i, j, k) = [v_hat(v) for v in span]

    return (i,j,k)

def active_image():
    '''Return the current active image, or Nothing if there are no active images.'''
    for i in mimics.data.images:
        if i.active:
            return i

def basis_vectors():
  return mimics_basis_vectors(active_image())

def bbox_from_points(p1, p2, basis=DEFAULT_BASIS):
    '''Given two points p1, p2 create a mimics.BoundingBox3D between them.'''
    span = np.array(p2) - np.array(p1)
    vectors = span * basis
    return mimics.BoundingBox3d(p1, vectors[0], vectors[1], vectors[2])

def bbox_to_points(bbox):
    """Find the extreme points p1 and p2 of a mimics.BoundingBox3D, with p1 at the origin."""
    p1 = bbox.origin
    span = np.array(bbox.first_vector) + np.array(bbox.second_vector) + np.array(bbox.third_vector)
    p2 = np.array(p1) + span
    return p1, p2

def expand_bbox_vector(bbox, expand, basis=DEFAULT_BASIS):
  """Expand a mimics.BoundingBox3D by adding a vector = (X_left, X_right), (Y_ant, Y_post), (Z_inf, Z_sup)."""
  # Rearrange the expansion values for easier calculation
  # so that they are ordered (min(X, Y, X), max(X, Y, Z))
  exp_min, exp_max = [idx for idx in zip(* expand)]
  
  p1, p2 = bbox_to_points(bbox)

  # Subtract min(expand) from the origin, by X,Y,Z component
  # Multiply the expansion vector by the basis to allow for skew scans
  new_p1 = np.array(p1) - np.diag(np.array(exp_min) * np.array(basis))
  new_p2 = np.array(p2) + np.diag(np.array(exp_max) * np.array(basis))

  return bbox_from_points(new_p1, new_p2)

def expand_bbox_points(bbox, exp_min, exp_max, basis=DEFAULT_BASIS):
  """Expand a mimics.BoundingBox3D by adding two offset vectors (X_left, Y_ant, Z_inf) & (X_right, Y_post, Z_sup)."""

  p1, p2 = bbox_to_points(bbox)

  # Subtract min(expand) from the origin, by X,Y,Z component (or add to the max point)
  # Multiply the expansion vector by the basis to allow for skew scans
  new_p1 = np.array(p1) + np.diag(np.array(exp_min) * np.array(basis))
  new_p2 = np.array(p2) + np.diag(np.array(exp_max) * np.array(basis))

  return bbox_from_points(new_p1, new_p2)

def make_crop_box(pt_a, pt_b, xy_scale, depth, thickness, basis=DEFAULT_BASIS):
  """Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
  Input:  two points for the posterior edge of the box, an xy scale, a thickness and a y depth.
          The posterior edge will be extended by xy_scale in both "X" directions and by depth anteriorly (Y).
          The box will be thickness in the Z direction.
          If the basis vectors are given, use those to adjust the directions of X, Y & Z.
  Output: A mimcs.BoundingBox3d to use for cropping a mask."""

  # numpy version
  k = (xy_scale - 1)/2.0
  delta = np.array(pt_a) - np.array(pt_b)
  p1 = np.array(pt_a) + (k * delta)
  p2 = np.array(pt_b) - (k * delta)
  origin = np.array(p1) - np.array(basis[Z]) * (thickness/2) # offset down to centre thickness on the plane
  first_vector = np.array(p2) - np.array(p1)
  second_vector = depth * np.array(basis[Y])
  third_vector = thickness * np.array(basis[Z])
  
  return mimics.BoundingBox3d(origin, first_vector, second_vector, third_vector)

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

def get_sides(parts):
  '''Given a list of 0 or more mimics objects return 0 to 2 of them allocated to the correct side of the head.'''
  sides = {'left': None, 'right': None}
  try:
    p0 = get_centre(parts[0])[0]
  except IndexError:
    # there are none of this part, so return both sides as None
    return sides
  
  try:
    p1 = get_centre(parts[1])[0]
  except IndexError:
    # There was one of this part, so compare centre.X to 0
    if p0 < 0:
      sides['right'] = parts[0]
    else:
      sides['left'] = parts[0]
  else:
    # Given two centres, so compare them to work out sides (0 = False, 1 = True)
    left_idx = int(p0 < p1) # If p0 < p1 then p0 is right and p1 is left, so return 1, and vice-versa.
    right_idx = 1 - left_idx
    sides['left'] = parts[left_idx]
    sides['right'] = parts[right_idx]
    
  return sides

def label_eyes():
    '''Check if there are two each of globes, rims and apexes'''
    globes = get_sides(mimics.data.spheres)
    for k,v in globes.items():
        v.name = f"{k}_globe"
    
    rims = get_sides(mimics.data.splines)
    for k,v in rims.items():
        v.name = f"{k}_rim"

    apexes = get_sides(mimics.data.points)
    for k,v in apexes.items():
        v.name = f"{k}_apex"
    
    return {"globe": globes, "rim": rims, "apex": apexes}

# To get all parts of one side use 
# eyes = label_eyes()
# {g['right'] for g in eyes.values()}

def draw_bbox(bbox, id = "bb"):
    '''Draw the outline of the given bounding box.'''
    mimics.analyze.create_point(bbox.origin, name=f"{id}_origin", color=(0.3,0.7,0.8))
    mimics.analyze.create_line(point1=bbox.origin, point2=np.array(bbox.origin) + np.array(bbox.first_vector),  name=f"{id}_first")
    mimics.analyze.create_line(point1=bbox.origin, point2=np.array(bbox.origin) + np.array(bbox.second_vector), name=f"{id}_second")
    mimics.analyze.create_line(point1=bbox.origin, point2=np.array(bbox.origin) + np.array(bbox.third_vector),  name=f"{id}_third")

def expand_bbox_old(bb: mimics.BoundingBox3d, p1_offset, p2_offset):
    '''Expand a bounding box bb by adding p1_offset to the origin and p2_offset to the opposite corner.'''
    side_extra = np.array(p2_offset) - np.array(p1_offset) # Change in length of each side
    sides = (bb.first_vector, bb.second_vector, bb.third_vector) # existing sides
    sizes = (np.linalg.norm(v) for v in sides) # the length of each existing side

    m = 1 + np.array(side_extra)/np.array(sizes) # multiplier for new vector lengths

    new_bb = mimics.BoundingBox3d()
    new_bb.origin = bb.origin + np.array(p1_offset)

    new_bb.first_vector  = m[0] * np.array(bb.first_vector)
    new_bb.second_vector = m[1] * np.array(bb.second_vector)
    new_bb.third_vector  = m[2] * np.array(bb.third_vector)

    return new_bb

def boolean_list(mask_list, op="Unite"):
    del_list = list()
    mask_a = mask_list[0]
    for mask_b in mask_list[1:]:
        mask_a = mimics.segment.boolean_operations(mask_a=mask_a, mask_b=mask_b, operation=op)
        del_list.append(mask_a)

    del del_list[-1]
    for m in del_list:
        mimics.data.masks.delete(m)

    return mask_a

def make_masks(bb):
    air_lo = mimics.segment.HU2GV(-1024)
    air_hi = mimics.segment.HU2GV(-200)
    air_all = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                       threshold_min=air_lo, threshold_max=air_hi, 
                                       bounding_box=bb)
    air_all.name = "air_all"
    air_all.visible = False

    # The exterior air should be the biggest connected region of air
    air_big = mimics.segment.keep_largest(mimics.data.masks.duplicate(object=air_all))

    air_big.name = "air_ext"
    # The interior air is everything else, but can include disconnected pockets in the injured eye
    air_int = mimics.segment.boolean_operations(mask_a=air_all, mask_b=air_big, operation='Minus')

    # Don't include the diconnected pockets in the interior air.
    air_int = mimics.segment.keep_largest(air_int)
    air_int.name = "air_int"

    air_shell_pixels = 1
    # Make a shell around the air, to represent very thin bone that doesn't show up on CT very well
    air_shell = mimics.segment.morphology_operations(input_mask=air_int, 
                                                 operation='Dilate', number_of_pixels=air_shell_pixels, connectivity=26, 
                                                 target_mask_name="air_shell")
    bone_lo = mimics.segment.HU2GV(226)
    bone_hi = mimics.segment.HU2GV(2976)
    bone_all = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                        threshold_min=bone_lo, threshold_max=bone_hi, 
                                        bounding_box=bb)
    bone_all.name = "bone_all"

    # Add the shell around the air to the thresholded bone
    bone_smooth = mimics.segment.boolean_operations(mask_a=bone_all, mask_b=air_shell, operation='Unite')
    # Remove floating bits from the bone + air shell mask
    mimics.segment.smooth_mask(bone_smooth) # NB - doesn't create a new mask, so no temp/delete step
    mimics.segment.keep_largest(bone_smooth)
    bone_smooth.name = "bone_smooth"

    mimics.data.masks.delete([air_all, air_shell, bone_all])  # Delete temporary masks
    
    # Add low density bone using local_threshold
    ### NB - local_threshold takes threshold values in HU, not GV like everything else ###
    bone_search_dist = 2
    bone_local = mimics.segment.local_threshold(mask=bone_smooth, threshold_min=148, threshold_max=661, search_distance=bone_search_dist, isolate=True, bounding_box=bb)
    bone_temp = mimics.segment.boolean_operations(mask_a=bone_smooth, mask_b=bone_local, operation='Unite')
    mimics.data.masks.delete(bone_local) # Delete a temporary mask - comment out for debugging
    
    # and smartfill to close some holes
    bone_orbit = mimics.segment.smart_fill_global(mask=bone_temp, hole_closing_distance=3)
    mimics.data.masks.delete(bone_temp) # Delete a temporary mask
    bone_orbit.name = "bone_orbit"

    # Convert the bone to a part, then Smooth and Wrap it to close more gaps
    bone_temp = mimics.segment.calculate_part(mask=bone_orbit, quality='High')
    bone_temp = mimics.tools.smooth(bone_temp, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=False)
    bone_temp.name = "bone_smooth"
    bone_part = mimics.tools.wrap(bone_temp, smallest_detail=0.2, gap_closing_distance=6, dilate_result=False, protect_thin_walls=True, keep_originals=False)
    bone_part.name = "bone_part"
    bone_part.visible = False
    #mimics.data.parts.delete(bone_temp)
    bone_wrap = mimics.segment.calculate_mask_from_part(part=bone_part, target_mask=None)
    bone_wrap.name = "bone_wrap"

    return air_big, air_int, bone_smooth, bone_orbit, bone_wrap # External and Internal air and the augmented bone mask (bone + internal air + extra shell)

def make_anterior_mask(rim, globe):
  """Calculate Orbital contents volumes.
  Input:  a mimics.Spline delmiting the rim, and a mimics.Sphere for the globe.
  Output: a Mask from the 'surface' of the rim, extending forward, for cropping."""
  
  # This function is slow, as it creates a lot of geometry. 
  # Consider usinge with mimics.disabled_gui(): to turn the GUI off and speed things up.

  # Get the BoundingBox3D of the rim to calculate a mask to cover it
  bbox_rim = mimics.measure.get_bounding_box([rim]) # ([rim]) works, but (rim) doesn't
  
  delta_z = bbox_rim.third_vector[Z]       # The height of the orbital rim
  min_z = np.floor(bbox_rim.origin[Z]) - 1 # The bottom of the orbital rim, minus one to ensure coverage
  n_planes = int(np.ceil(delta_z)) + 2     # place a plane every mm, plus one below and one more above the rim

  # spacing = delta_z / n_planes     # Distance between planes to cover the whole rim
  thickness = 2 # not spacing, but the height of each mask, which is 1 unit above and below the plane
 
  # Make a list of plane origin points, matching the globe in the X,Y plane 
  # and spaced in Z from min_z to max_z
  plane_origins = [(globe.center[X], globe.center[Y], min_z + z) for z in range(n_planes)]
  # Create planes using these origins, parallel to the X,Y plane (normal is Z+, [0, 0, 1])
  # To allow for skewed data get the correct normal from the active image set
  norm_z = basis_vectors()[Z]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(o, norm_z) for o in plane_origins]
  # Hide all the planes for speed of rendering
  for p in mimics.data.planes:
    p.visible = False
  
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
        #mimics.data.lines.delete(line_int) # remove the temp line
      if (p1[Z] <= plane_z and p2[Z] > plane_z): # crosses from below
        line_int = mimics.analyze.create_line(p1, p2) # temp line to get intersection point
        pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line_int, plane)
        #mimics.data.lines.delete(line_int) # remove the temp line

      # Each segment will only cross a given plane once. Stop when have found one in each direction.
      if (pt_up is not None) and (pt_down is not None):
        # Create a bounding box and add it to the list
        #-#bb = utils_2.bbox_from_intersections(pt_up, pt_down, MULT_XY, thickness, SIZE_Y)
        bb = make_crop_box(pt_up, pt_down, MULT_XY, SIZE_Y, thickness, basis)
        boxes.append(bb)
        # Clean up the used points. Do this here to avoid deleting the apex points.
        mimics.data.objects.delete([pt_up, pt_down])
        break # We have found that line for this plane   
    else:
      # Fell through the loop without breaking
      print(f"WARNING: did not find intersection for plane {plane}")
  
  # Clean up CAD objects
  #mimics.data.objects.delete(z_planes) # delete only the planes in a list
  mimics.data.objects.delete(mimics.data.lines)  # delete all the lines that exist
  mimics.data.objects.delete(mimics.data.planes) # delete all the planes that exist

  # Thicken the first and last bounding boxes to ensure full overlap.
  add_z = np.array((0, 0, EXTRA_Z))
  boxes[0]  = expand_bbox_points(boxes[0], -add_z, (0,0,0), basis) # first goes down
  boxes[-1] = expand_bbox_points(boxes[-1], (0,0,0), add_z, basis) # last goes up

  # Create a series of masks, covering all materials, clipped by each bounding box
  # As each mask is created, unite it into a combined mask to cut away tissue anterior to the rim.
  mask_ant = mimics.segment.create_mask(select_new_mask=False)
  for b in boxes:
      mask_ab = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                        threshold_min=materials.MIN_GV,
                                        threshold_max=materials.MAX_GV, 
                                        bounding_box= b)
      old_ant = mask_ant # unite creates a new mask, so save the current one to delete after.
      mask_ant = mimics.segment.boolean_operations(mask_ant, mask_ab, operation="Unite")
      mimics.data.masks.delete(mask_ab)
      mimics.data.masks.delete(old_ant)
  
  return mask_ant


def process_project():
    '''Draw the expanded bonding box for the current Project and save a snapshot.'''

    t0 = time.perf_counter()
    print(f"start_process {t0}")

    # Clean up the display, starting by turning everything off
    ### For some reason they are visible when this runs? Not sure why.
    for o in mimics.data.objects:
      o.visible = False
    
    for m in mimics.data.masks:
      m.visible = False

    # Find the two eyes, label them and make them visible
    eyes = label_eyes()
    # to get all parts of one eye use {g['right'] for g in eyes.values()}
    for part in eyes.values():
      for side in part.values():
          side.visible = True

    # The eyes are defined by the first two splines, the first two globes, and the first two points.
    # There were some scans with no points defined, but they are very useful so I went back and added them where missing.

    # Get the bounding box of the eye objects
    bb = mimics.measure.get_bounding_box([
        mimics.data.splines[0], mimics.data.splines[1],
        mimics.data.spheres[0], mimics.data.spheres[1],
        mimics.data.points[0],  mimics.data.points[1]
    ])

    # Expand this bounding box to make sure all of the orbits are covered.
    # Use a big anterior (-ve Y) expansion to ensure that external_air is the biggest part of the air mask.
    p1_offset = (-10, -50, -15) # added to the origin of the bounding box. Made the Z deeper to get more sinus for a blowout fracture.
    p2_offset = (10, 10, 15)    # added to the opposite point of the bounding box. Made Z taller as first case was just clipping top of orbit.

    bb_all = expand_bbox_points(bb, p1_offset, p2_offset, basis=basis)

    # Turn off the GUI for the autoamated parts to speed it up
    with mimics.disabled_gui_update():
      # Draw a mask for this bounding box
      #all = mimics.segment.threshold(mask=mimics.segment.create_mask(),threshold_min=0,threshold_max=4096,bounding_box=bb_all)
      #all.name = "all"
      #all.visible = False
      
      # Create masks for the externa and internal air and for the bone, smoothed, filled and wrapped
      air_ext, air_int, bone_smooth, bone_orbit, bone_wrap = make_masks(bb_all)

      # Hide most of the masks
      for m in (air_ext, air_int, bone_smooth, bone_orbit, bone_wrap):
        m.visible = False

      ## Run the make_anterior_mask for the left and right side 
      anterior_right = make_anterior_mask(eyes['rim']['right'], eyes['globe']['right'])
      anterior_right.name = "anterior_right"
      anterior_right.visible = False
      anterior_left = make_anterior_mask(eyes['rim']['left'], eyes['globe']['left'])
      anterior_left.name = "anterior_left"
      anterior_left.visible = False

      ## Boolean together anterior masks
      anterior_both = mimics.segment.boolean_operations(anterior_left, anterior_right, 'Unite')
      anterior_both.name = "anterior_both"
      anterior_both.visible = False

      ## Create Soft Tissue masks - Threshold (lower air upper bone)
      air_hi = mimics.segment.HU2GV(-200)
      bone_lo = mimics.segment.HU2GV(226)
    
      soft_all = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                        threshold_min=air_hi, threshold_max=bone_lo, 
                                        bounding_box=bb)
      soft_all.name = "soft_all"
      soft_all.visible = False
      
      ## Boolean subtract Orbital Volume = Soft Tissue - Anterior - (Bone + Air).
      soft_temp = mimics.segment.boolean_operations(soft_all, anterior_both, 'Minus')
      soft_orbits = mimics.segment.boolean_operations(soft_temp, bone_wrap, 'Minus')
      soft_orbits.name = "soft_orbits"
      soft_orbits.visible = False
      mimics.data.masks.delete(soft_temp) # Delete a temporary mask - comment out for debugging
      
      mimics.data.masks.delete([anterior_both, anterior_left, anterior_right, soft_all]) # Delete more temporary masks - comment out for debugging

    # GUI active again

    t1 = time.perf_counter()
    print(f"end_process {t1}\telapsed {t1-t0}")
    
    ## Call Edit Mask - to manually fix segmentation leaks.
    soft_clean = mimics.data.masks.duplicate(object=soft_orbits)
    soft_clean.name = "soft_clean"

    # Show the mask we want and turn on 3D preview mode
    for m in mimics.data.masks:
      m.selected = False
    soft_clean.selected = True
    mimics.view.enable_mask_3d_preview()

    soft_clean = mimics.segment.activate_edit_mask(mask = soft_clean, edit_mode = "Erase", edit_type = "Ellipse")

    ## Create the globe masks and add to the soft_clean mask
    globe_left = utils.sphere_to_mask(eyes['globe']['left'])
    globe_left.name ="left_globe"
    globe_left.visible = False
    globe_right = utils.sphere_to_mask(eyes['globe']['right'])
    globe_right.name ="right_globe"
    globe_right.visible = False

    ## Region Grow Volume (L) from L globe centre. 
    # Add the globe to the mask then grow from the globe centre. Need this in case the centre is in front of the anterior surface.
    temp_union = mimics.segment.boolean_operations(globe_left, soft_clean, 'Unite')
    temp_union.name = "temp_union_1"
    globe_point = get_centre(eyes['globe']['left'])
    temp_left = mimics.segment.region_grow(input_mask=temp_union, target_mask=None, point=globe_point, slice_type="Axial", keep_original_mask=False, multiple_layer=True, connectivity='26-connectivity')
    temp_left.name = "temp_left"
    mimics.data.masks.delete(temp_union) # need to delete here as re-using the name below creates a new copy of the mask
    ## Subtract Globe again to leave the orbital volume
    orbit_left = mimics.segment.boolean_operations(temp_left, globe_left, 'Minus')
    orbit_left.name = "left_Orbital Volume"

    ## Repeat for Right side. 
    temp_union = mimics.segment.boolean_operations(globe_right, soft_clean, 'Unite')
    temp_union.name = "temp_union_2"
    globe_point = get_centre(eyes['globe']['right'])
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
    
    ##Create empty masks called "Manual_Muscle+Nerve" and "Manual_Haematoma"


    # Create

    # Create an "all soft tissue" mask, then subtract the orbit_surrounds mask



def clean_project_name(project_name):
    '''Create a consistent name from the project name.'''
    (stem, date, series) = re.match(pattern="(.*) ([0-9\\.]+) ([A-Za-z]+) SS 01.mcs", string=project_name).groups()
    series = series.casefold()
    if series[0:5] == "merge":
      series = "merged"
    orbit_name = f"{stem}_{date}_{series}_orbit.mcs"
    return orbit_name

def point_info(pt):
    return (pt.name, pt[0], pt[1], pt[2])

def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def show_3D_obj(objs):
      # Get the combined bounding box to zoom to
    orbital_bb = mimics.measure.get_bounding_box(objects=objs)
    # Zoom each view to cover both masks
    for v in mimics.data.views: 
        view_cam = mimics.view.get_camera(v)
        view_settings = view_cam.get_settings()
        view_settings.zoom_to_bounding_box(orbital_bb, zoom_factor=0.8)
        view_cam.set_settings(view_settings)
    
    # Save a screenshot for this project
    # Use the 3D view
    settings = mimics.view.get_camera(view = mimics.data.views['3D']).get_settings()
    settings.zoom_to_bounding_box(orbital_bb, zoom_factor=1)
    pic_filename = filename.replace('processed.mcs', 'snapshot.jpg')
    mimics.file.export_view_by_type(filename=os.path.join(folder, pic_filename), view='3D', image_type = 'jpg', camera_settings=settings)


if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  # This version has one folder with the segmenting person in the file name
  
  #root = r'D:\Projects & Research\Enophthalmos Study\re-do_DICOM'
  root = r'D:\Projects & Research\Enophthalmos Study'
  
  # Put the combined results in the root, rather than one file per user.
  results_file = Path(os.path.join(root, 'project_points.csv'))
  # Where to put the output files?
  output_folder = root

  # Get a list of all the .mcs files in root
  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]
  
  # Get a list of all the Mimics files that have not been done yet.
  # column 1 is the file name, column 2 is whether the file has been procesessed.
  with open(os.path.join(root, 'files_to_analyse.csv'), newline = '') as csvfile:
    files_to_analyse = [row[0] for row in csv.reader(csvfile) if row[1] == "FALSE"]
  
  # Convert to full paths
  projects = [os.path.join(root, f) for f in files_to_analyse]

  # Only look at a subset to test
  projects = projects[0:3]

  num_projects = len(projects)
  print(f'processing {num_projects} project files from {root}')

  # Make logging quieter to hopefully speed things up?


  #with mimics.disabled_gui_update():
  for i, p in enumerate(projects):
      mimics.file.open_project(filename=p, read_only_mode=True)
      
      print(f'project {i+1}/{num_projects}\t"{os.path.basename(p)}"')
      
      #log_to_file(results_file, headers = ["file", "p0", "x", "y", "z", "p1", "x", "y", "z"],
      #                          results = flatten((os.path.basename(p), two_points))
      #)

      # Set a global with the basis vectors of the current image set, in case it is skewed.
      basis = basis_vectors()
      

      process_project()

      # Make a picture file name by changing the extewnsion, then save the snapshot
      #picture_path = re.sub(pattern="mcs$", repl="bbox.png", string=p)
      #mimics.file.save_screenshot(picture_path)

      path_name = os.path.join(output_folder, clean_project_name(p)) # full path for output file
      mimics.file.save_project(path_name)

      mimics.file.close_project()
