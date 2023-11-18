
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
  	part_bone = mask_to_part("Bone", mask_bone)
  
  if TESTING: ##########################################################
    # Fake mimics geometry objects
    class Sphere:
      def __init__(self, x, y, z, radius):
        self.radius = radius
        self.center = [x, y, z]
    
    class Spline:
      def __init__(self, points):
        self.points = points

    # Define the geometry that the user would normally enter
    orbit = Sphere(11.776, -29.9621, -36.1569, 76.0900)
    spline_pts = [[-22.31, -46.94, 94.00],
                  [-32.22, -46.54, 92.10],
                  [-41.65, -42.97, 87.38],
                  [-46.08, -38.57, 83.18],
                  [-45.93, -35.38, 75.76],
                  [-46.55, -35.92, 68.80],
                  [-44.12, -40.42, 61.95],
                  [-36.31, -44.23, 58.72],
                  [-26.31, -44.84, 60.14],
                  [-16.61, -46.71, 64.92],
                  [-11.09, -48.23, 70.07],
                  [-8.86, -46.85, 76.10],
                  [-10.93, -44.75, 82.40],
                  [-13.85, -45.99, 88.07],
                  [-17.94, -46.54, 92.42,]
                 ]
    spline1 = Spline(points = spline_pts)
  else: ################################################################
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
  ######################################################################
  
  
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

    #delete all masks (except the Orbital Volume mask...)
    objects_to_keep = ("Bone", "Spline 1", "Sphere 1", "Orbital Volume")
    for m in mimics.data.objects:
      if m.name in objects_top keep:
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
 



###############################

# Ryan version
  intersect_points = {} # here len(intersect_points) == 0, so using len()+1 below starts indices from 1
  for n in range(len(z_planes)):
    for m in range(len(spline_lines)):
      if spline_lines[m].point1[2] > z_planes[n].origin[2] and spline_lines[m].point2[2] < z_planes[n].origin[2]:
        intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n])
      if spline_lines[m].point1[2] < z_planes[n].origin[2] and spline_lines[m].point2[2] > z_planes[n].origin[2]:
        intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n])

# Rob version
  intersect_points_RED = [] # List to store intersection points.
  for n, plane in enumerate(z_planes):
    # Go around the spline lines and for each line segment check if the y values span the plane.
    # If yes, then calculate the intercept. Going around the whole line gets both the upwards
    # and downwards intersections, one following the other so each pair of points is for the same plane
    for line in spline_lines:
      pt1_z = line.point1[Z]
      pt2_z = line.point2[Z]
      plane_z = plane.origin[Z]
      from_above = (pt1_z > plane_z and pt2_z < plane_z)
      from_below = (pt1_z < plane_z and pt2_z > plane_z)
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z] )
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z] )
      if (from_above or from_below):
        intersect_points_RED.append(mimics_analyze_create_point_as_line_and_plane_intersection(line, plane))

# Compare versions. 
# Ryan's is a dict so need to unpack the values to a list first
  all([i.__dict__ == j.__dict__ for i, j in zip([v for v in intersect_points.values()], intersect_points_RED)])
## Identical (all True)

# Rob version - by planes
  intersect_points_RED2 = len(z_planes) * [None] # List to store intersection points.
  for n, plane in enumerate(z_planes):
    # Go around the spline lines and for each line segment check if the y values span the plane.
    # If yes, then calculate the intercept. Going around the whole line gets both the upwards
    # and downwards intersections, one following the other so each pair of points is for the same plane
    for line in spline_lines:
      # Does ths line cross the curent plane? If so, in which direction?
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z])
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z])
      # Add the corresponding point to the list
      if from_above:
        pt_down = intersect_points_RED.append(mimics.analyze.create_point_as_line_and_plane_intersection(line, plane))
      if from_below:
        pt_up = intersect_points_RED.append(mimics.analyze.create_point_as_line_and_plane_intersection(line, plane))
    intersect_points_RED2[n] =   


# Ryan version 
  crop_masks = {} 
  for n in range(len(intersect_points)):   
    if n % 2 == 1:
   # If I read this correctly, this makes a line from points (1,2), (3, 4) etc.
   # What happened to point zero? ANS: we lost it in creating the list, as used [n+1] and n starts from 0
      vector_x = [(intersect_points[n].x-intersect_points[n+1].x)*1.2,(intersect_points[n].y-intersect_points[n+1].y)*1.2,0]
      print(f"vector_x from ({intersect_points[n].x}, { intersect_points[n+1].x}) and ({intersect_points[n].y}, {intersect_points[n+1].y})")
      vector_y = [0,-20,0]    
      if n == 1:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z+spacing_z/2)         
        vector_z = [0,0,-1.5*spacing_z]        
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)) 
      elif n == len(intersect_points)-1:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)         
        vector_z = [0,0,1.5*spacing_z]        
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
      else:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)          
        vector_z = [0,0,spacing_z] 
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
  
  print("Ryan result")      
  print('\n'.join(map(str, [i.__dict__ for i in crop_masks.values()])))

# Rob version
  united_masks = mimics.segment.create_mask()
  crop_masks_RED = []
  point_pairs = list(batched(intersect_points_RED, 2))
  last_pair = len(point_pairs) - 1
  for i, (pt_up, pt_down) in enumerate(point_pairs):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0] # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # extend first and last masks in z beyond their planes to cover whole orbit
    if i == 0:
      vector_z = [0, 0, -1.5 * spacing_z] # extends down past start
    elif i == last_pair:
      vector_z = [0, 0, +1.5 * spacing_z] # extends up past end
    
    # Set the bounding box origin above the down intesection point
    # and midway between the z planes
    box_orig = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    crop_box = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    
    # create the mask, threshold and crop it, then append to the list
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, crop_box)
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    crop_masks_RED.append(m)
    
    
  print("Rob Result")
  print('\n'.join(map(str, [m.__dict__ for m in crop_masks_RED])))


  point_pairs = list(batched(intersect_points_RED, 2))
  def make_crop_box(pt_up, pt_down):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0]       # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # Put the origin above the down intesection point, midway between the z planes
    origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)
  
  boxes = [make_crop_box(pt_up, pt_down) for pt_up, pt_down in point_pairs]
  # Make the First and Last box a bit bigger to ensure coverage
  boxes[0].third_vector[2] = -1.5 * boxes[0].third_vector[2]
  boxes[-1].third_vector[2] = 1.5 * boxes[-1].third_vector[2]

  # Create the united mask
  united_masks = mimics.segment.create_mask()
  for bbox in boxes:
    # Create a mask, threshold to cover everything
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, bbox) # Crop it with BoundingBox for this plane
    # Then Union with the existing mask
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    
  
  last_pair = len(point_pairs) - 1
  for i, (pt_up, pt_down) in enumerate(point_pairs):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0] # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # extend first and last masks in z beyond their planes to cover whole orbit
    if i == 0:
      vector_z = [0, 0, -1.5 * spacing_z] # extends down past start
    elif i == last_pair:
      vector_z = [0, 0, +1.5 * spacing_z] # extends up past end
    
    # Set the bounding box origin above the down intesection point
    # and midway between the z planes
    box_orig = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    crop_box = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    
    # create the mask, threshold and crop it, then append to the list
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, crop_box)
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    crop_masks_RED.append(m)
    
