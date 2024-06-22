
# Simple 3D geometry functions
from math import sqrt # needed to calculate vector magnitude

def dot(u, v):
  return [u[0] * v[0], 
          u[1] * v[1], 
          u[2] * v[2]
         ]

def mag(u):
  return sqrt(u[0]**2 + u[1]**2 + u[2]**2)

def cross(u, v):
  return [u[1]*v[2] - u[2]*v[1],
          u[2]*v[0] - u[0]*v[2],
          u[0]*v[1] - u[1]*v[0]
         ]

def norm_from_points(A, B, C):
  u = [C[0] - A[0], C[1] - A[1], C[2] - A[2]]
  v = [B[0] - A[0], B[1] - A[1], B[2] - A[2]]
  return cross(u, v)

def clockwise(curve):
  ofs = len(curve) // 3
  (A, B, C) = [curve[p] for p in (0, ofs, 2 * ofs)] # get spaced points
  norm = norm_from_points(A, B, C)
  return norm[1] < 0

# def side(point):
#   if point[1] < 0:
#     side = 'left'
#   elif point[1] > 0:
#     side = 'right'
#   else:
#     side = 'ambiguous'
#   return side

## Only needed for dummy routines
import mimics
from mimics import segment
from mimics import analyze
from mimics import data

import const # Contains definitions of all materials

# Define parameters for planes and crop boxes to trim front of orbit
NUM_PLANES = 10	# Number of planes to divide up orbit
MULT_XY = 1.2   # Factor to expand the XY extents of the crop box
MULT_Z = 1.5    # Factor to expand the Z extent of the first & last crop boxes
SIZE_Y = -20    # Y extent of crop boxes

CLOSING_DIST = 7 # Smartfill 

# Segment orbital contents into these Materials & measure volume
materials = {
  'air'    : const.MATL_AIR, 
	'fat'    : const.MATL_FAT,
	'muscle' : const.MATL_MUSCLE
}

from utils import *


# Main function -  called if file is run and not imported
if __name__ == '__main__':
  
  # If there is no mask called "Bone Mask" then make one & a corrresponding part.
  if mimics.data.objects.find("Bone Mask") is None:
  	mask_bone = material_mask("Bone Mask", const.MATL_BONE)
  	part_bone = part_from_mask("Bone", mask_bone)
  
  #User Inputs
  if mimics.data.objects.find("Spline 1", False) == None:
    spline1 = mimics.analyze.indicate_spline(
      message='Indicate spline on the orbital rim, passing through landmarks', 
      show_message_box=True, confirm=False, title=None
     )
  if mimics.data.objects.find("Sphere 1", False) == None:
     orbit = mimics.analyze.indicate_sphere(
      message="Indicate the globe using 3pts on the axial view", 
      show_message_box=True, confirm=True, title=None
    )

  
  spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) for a, b in looped_pairwise(spline1.points)]

  (max_x, max_y, max_z) = [max(idx) for idx in list(zip(* spline1.points))]
  (min_x, min_y, min_z) = [min(idx) for idx in list(zip(* spline1.points))]
 
  print("Spline limits are:")
  print(f"X {min_x} to {max_x}\nY {min_y} to {max_y}\nZ {min_z} to {max_z}") 
 
  delta_z = round(max_z - min_z, 0)
  spacing_z = delta_z/NUM_PLANES
  print(f"Z span {delta_z} divided by {const.NUMPLANES} = spacing_z {spacing_z}")
  
  # Create a list of plane origin points, each aligned to the orbit in the X,Y plane 
  # and spaced in Z from min_z to min_z + delta_z (approximately max_z)
  orig_x =  orbit.center[0]
  orig_y =  orbit.center[1]
  plane_origins = [[orig_x, orig_y, min_z + (z + 1) * spacing_z] for z in range(const.NUM_PLANES - 1)]
  # Create planes using these origins paralle to the X,Y plane (normal is Z+)
  norm_z = [0, 0, 1]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) for orig in plane_origins]

  # Function to create a mimics.BoundingBox3D given two plane intersection points.
  def make_crop_box(pt_up, pt_down):
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ MULT_XY * (pt_up.x - pt_down.x), MULT_XY * (pt_up.y - pt_down.y), 0]
    vector_y = [0, SIZE_Y, 0]    # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # Put the origin above the down intesection point, midway between the z planes
    origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)
  
  # Create a blank list for the boxes
  boxes = []
  # For each plane, find intersections with the spline lines and create a 
  # bounding box based on the intersection points.
  for plane in z_planes:
    pt_up, pt_down = None, None
    for line in spline_lines:
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z])
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z])
      if (from_above):
        pt_up = mimics_analyze_create_point_as_line_and_plane_intersection(line, plane)
      if (from_below):
        pt_down = mimics_analyze_create_point_as_line_and_plane_intersection(line, plane)

    if from_above and from_below:
      boxes.append(make_crop_box(pt_up, pt_down))
    else:
      print(f"did not find intersection for plane {plane}")

  # Adjust the first and last bounding box to ensure full overlap.
  boxes[0].third_vector[2] = -MULT_Z * boxes[0].third_vector[2]  # first goes down
  boxes[-1].third_vector[2] = MULT_Z * boxes[-1].third_vector[2] # last goes up
   
  # Create the combined mask
  united_masks = mimics.segment.create_mask()
  # Loop through each plane again, create a mask, crop it, and Union it together.
  for bbox in boxes:
    # Create a mask, threshold to cover everything
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    # Crop it with BoundingBox for this plane
    mimics.segment.crop_mask(m, bbox)
    # Union with the existing mask
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    #united_masks = add_masks(united_masks, m)
    
    pause 
            
    mask_bone = material_mask("Bone Mask", const.MATL_BONE) # Already made this above???
    mask_air = material_mask("Air Mask", const.MATL_AIR)
    
    # Union complete crop mask with air and bone
    united_masks = mimics.segment.boolean_operations(united_masks, mask_bone, 'Unite')
    united_masks = mimics.segment.boolean_operations(united_masks, mask_air,  'Unite')
    #united_masks = add_masks(united_masks, mask_bone)
    #united_masks = add_masks(united_masks, mask_air)
    
    smartfill_mask = mimics.segment.smart_fill_global(united_masks, CLOSING_DIST)
    
    # Define the extents of the region to analyse, based on which side is being analyzed.
    # This assumes a LPS coordinate system and that the midline is x == 0.
    if   orbit.center[1] < 0: # Right side
      x_lim = min_x - 10
      x_delta = (max_x - min_x) + 20
    elif orbit.center[1] > 0: # Left side
      x_lim = max_x + 10
      x_delta = (min_x - max_x) - 20
    else:                     # On midline
      print(f"ERROR: Can't tell which side to analyse (orbit at { orbit.x,  orbit.y,  orbit.z}")
    box_orig = (x_lim, min_y - 5, min_z -15)
    vector_x = [x_delta, 0, 0]
    vector_y = [0, 80, 0]
    vector_z = [0, 0, (max_z - min_z) + 30]

    orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    
    smartfill_mask = mimics.segment.crop_mask(smartfill_mask, orbit_ROI)
    
    part_temp = mimics.segment.calculate_part(smartfill_mask, 'High')
    part_wrap = mimics.tools.wrap(part_temp, 0.2, 10, False, True, True)
    wrapped_mask = mimics.segment.calculate_mask_from_part(part_wrap, None)
    
    orbit_vol = mimics.segment.boolean_operations(wrapped_mask, smartfill_mask, 'Minus')
    orbit_vol = mimics.segment.morphology_operations(orbit_vol, 'Erode', 1, 8, None, None)
    orbit_vol = mimics.segment.region_grow(orbit_vol, orbit_vol, sphere1.center, 'Axial', False, True, connectivity='6-connectivity') 
    orbit_vol.name = "Orbital Volume"

    # Delete all masks (except the Orbital Volume mask...)
    objects_to_keep = ("Bone", "Spline 1", "Sphere 1", "Orbital Volume")
    for m in mimics.data.objects:
      if m.name in objects_to_keep:
        print(f"Keeping {m.name} {m.type}")
      else:
        mimics.data.objects.delete(m)
    
    # Remove the Globe of the eye from the volume to analyse
    # Amazingly, there does not appears to be a way to convert an object to a mask.
    # convert sphere to mask, manually. 
    mimics.dialogs.message_box("Convert the sphere object into a mask and rename the mask 'Globe'. Click okay to continue", title=None, ui_blocking=False)
    globe_mask = mimics.data.masks.find("Globe", False)
    if globe_mask is not None:
        print("Globe Exists!")
        intersect_vol_mask = mimics.segment.boolean_operations(orbit_vol, globe_mask, 'Minus')
        intersect_vol_mask.name = "Intersect Mask"
    else:
        print("ERROR: Cant Find Globe Mask")
 
