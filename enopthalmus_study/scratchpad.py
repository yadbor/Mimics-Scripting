X, Y, Z = 0, 1, 2

# The DICOM co-ordinate system is defined (for a BIPED) as patient LPS. That is:
#    X+ to the Left hand side of the patient
#    Y+ to the Posterior (towards the back)
#    Z+ to the Superior (towards the head)
#
# Mimics appears to use the DICOM LPS co-ordinate system.
#
# Which eye is being analysed? The head is usually scanned in the centre 
# of the CT co-ordinate system, so X=0 is roughly the midline, meaning 
# the Right eye is X < 0 and the Left eye is X > 0
# TODO - a possibly more robust midline is the mean image X co-ordinate
if sphere1.x < 0:
  side = 'right'
elif sphere1.x > 0:
  side = 'left'
else:
  side = 'ambiguous'

import mimics # for sytax checker and dummy routines

###############################
# Old test code to create dummy data
# Supesceded by creating dummy mimics.analyze.indicate_{spline|sphere} code

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
    


# Part from triangles
# in Mimics API help under mimcs.segment
import numpy as np
p = mimics.data.parts[0]
v,t = p.get_triangles()
v = np.array(v)
t = np.array(t)
for i in range(len(v)):
    v[i] = v[i]+100
mimics.segment.create_part(v,t)

import numpy as np
# Create a Scale matrix
def scale_matrix(s):
    return np.array([[s, 0, 0],[0, s, 0], [0, 0, s]])
# Scale an object by s (in our case, r/10 for a 10 unit sphere)
def scale_object(p, s):
    tx = scale_matrix(s)
    v,t = p.get_triangles()
    v = np.array(v)
    t = np.array(t) # Don't actually need this line 
    vt = v.dot(tx.T)
    return mimics.segment.create_part(vt,t)


# Text crop box calculation
# Make a tuple of spans (min_x, max_x)

spans = ((10, 20), (-20, -10))
for mn, mx in spans:
  box_orig_pos = mn - 10
  vector_x_pos = mx - mn + 20
  box_orig_neg = mx + 10
  vector_x_neg = mn - mx - 20
  print(f"box ({mn}, {mx}) neg from {box_orig_neg} by {vector_x_neg} to {box_orig_neg + vector_x_neg}, pos from {box_orig_pos} by {vector_x_pos} to {box_orig_pos + vector_x_pos}")



def find_eyes():
  """Find the correct side for each eye digitised component (globe, rim & apex point)."""
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



import numpy as np

def value_safe(v, default = 0):
  if v:
    return v
  else:
    return default

def expand_bbox(bbox, X_exp=None, Y_exp=None, Z_exp=None):
  """Given a mimics.BoungingBox, expand it by the given amounts (lo, hi) in each direction."""
  lo = [0,0,0]
  if X_exp[0]:
    lo[0] = X_exp[0]
  if Y_exp[0]:
    lo[0] = Y_exp[0]
  if Z_exp[0]:
    lo[0] = Z_exp[0]

  return bbox

DEFAULT_BASIS= ((1, 0, 0), (0, 1, 0), (0, 0, 1))

i = [1, 0, 0]; j = [0, 1, 0]; k = [0, 0, 1]
basis = (i, j, k)

def make_bbox(base_point, offset, extents, basis=DEFAULT_BASIS):
  vectors = (np.array(extents) * np.array(basis))
  bbox = mimics.BoundingBox3d(origin = base_point - np.array(offset),
                              first_vector = vectors[0],
                              second_vector = vectors[1],
                              third_vector = vectors[2]
                              )
  return bbox

def bbox_to_points_old(bbox):
    """Find the extreme points p1 and p2 of a mimics.BoundingBox3D, with p1 at the origin."""
    p1 = bbox.origin
    # Add each vector to the origin in turn to get from P1 to P2
    p2 = [(a + b) for a, b in zip(bbox.origin, bbox.first_vector)]
    p2 = [(a + b) for a, b in zip(p2, bbox.second_vector)]
    p2 = [(a + b) for a, b in zip(p2, bbox.third_vector)]
    return p1, p2

def bbox_to_points(bbox):
    """Find the extreme points p1 and p2 of a mimics.BoundingBox3D, with p1 at the origin."""
    p1 = bbox.origin
    span = np.array(bbox.first_vector) + np.array(bbox.second_vector) + np.array(bbox.third_vector)
    p2 = np.array(p1) + span
    return p1, p2

def bbox_from_points(p1, p2):
    """Given two points p1 and p2 crerate a mimics.BoundingBox3D between them."""
    span = np.array(p1) - np.array(p1)
    bbox = mimics.BoundingBox3d(p1, span[0]], span[1], span[2])
    return bbox

def expand_points(p1, p2, expand, basis=DEFAULT_BASIS):
  return None

def expand_bbox(bbox, expand, basis=DEFAULT_BASIS):
  """Expand a mimics.BoundingBox3D by adding a vector = (X_left, X_right), (Y_ant, Y_post), (Z_inf, Z_sup)."""
  # Rearrange the expansion values for easier calculation
  # so that they are ordered (min(X, Y, X), max(X, Y, Z))
  exp_min, exp_max = [idx for idx in zip(* expand)]
  
  p1, p2 = bbox_to_points(bbox)

  # Subtract min(expand) from the origin, by X,Y,Z component
  # Multiply the expansion vector by the basis to allow for skew scans
  new_p1 = np.array(p1) - (np.array(exp_min) * np.array(basis))
  new_p2 = np.array(p2) + (np.array(exp_max) * np.array(basis))

  return bbox_from_points(new_p1, new_p2)

############################################

img = mimics.data.images[-1]

def get_basis_vectors(img, get_origin = False):
  origin = img.get_voxel_center([0, 0, 0])
  dims = img.get_voxel_buffer().shape
  i = img.get_voxel_center([dims[0]-1, 0, 0])
  j = img.get_voxel_center([0, dims[1]-1, 0])
  k = img.get_voxel_center([0, 0, dims[2]-1])
  span = [np.asarray(v) - np.asarray(origin) for v in (i, j, k)]
  basis = [component / np.linalg.norm(component, ord=1) for component in span]
  if get_origin:
    return basis, origin
  else:
    return basis


def get_centre(part):
  """Get the geometric centre of a mimics sphere, point or spline."""
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
  bbox = mimics.measure.get_bounding_box([part])
  p1 = np.array(bbox.origin)
  span = np.array(bbox.first_vector) + np.array(bbox.second_vector) + np.array(bbox.third_vector)
  return p1 + (span / 2)

def get_sides(parts):
  """Given a list of 0 to 2 mimics objects return them allocated to side of the head."""
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
    # Managed to retrieve both centres, so compare them to work out sides
    left_idx = int(p0 < p1)
    right_idx = 1 - left_idx
    sides['left'] = parts[left_idx]
    sides['right'] = parts[right_idx]
    
  return sides

def find_eyes(spheres, splines, points):
  """Given some spheres (globes), splines (rims) and points (apical points) find all eyes"""

  # make an empty dict of eye components
  eyes = {'num_eyes': 0,
          'left': {'globe': None, 'rim': None, 'apex': None},
          'right':{'globe': None, 'rim': None, 'apex': None}}
  
  n_eyes = len(spheres)
  if n_eyes != len(splines):
    # mis-match, so exit with an error
    raise ValueError(f'Number of globes {n_eyes} and rims {len(splines)} do not match.')
    return Nine

  if n_eyes < 1 or n_eyes > 2:
    # wrong number of eyes, so exit with an error
    raise ValueError(f'{n_eyes} is not 1 or 2.')
    return None

  sides = get_sides(spheres)
  eyes['left']['globe'] = sides['left']
  eyes['right']['globe'] = sides['right']

  sides = get_sides(splines)
  eyes['left']['rim'] = sides['left']
  eyes['right']['rim'] = sides['right']

  sides = get_sides(points)
  eyes['left']['apex'] = sides['left']
  eyes['right']['apex'] = sides['right']

  eyes['num_eyes'] = n_eyes

  return eyes

# loop version that may be better, but maybe less obvious?
def find_eyes_loop(spheres, splines, points):
  """Given a list of spheres (globes), splines (rims) and points (apical points) find all eyes"""

  # make an empty dict to hold the eye components
  eyes = {'num_eyes': 0,
          'left': {'globe': None, 'rim': None, 'apex': None},
          'right':{'globe': None, 'rim': None, 'apex': None}}  
  
  n_eyes = len(spheres)

  if n_eyes != len(splines):
    # mis-match, so exit with an error
    raise ValueError(f'Number of globes {n_eyes} and rims {len(splines)} do not match.')
    return Nine

  if n_eyes < 1 or n_eyes > 2:
    # wrong number of eyes, so exit with an error
    raise ValueError(f'{n_eyes} is not 1 or 2.')
    return None
  # Gather the parts to check into a dict
  parts = {'globe': spheres, 
           'rim': splines,
           'apex': points}
  
  for part_name, part in parts.items():
    sides = get_sides(part) # Return the correct part or None for each side
    for side in ('left', 'right'):
      eyes[side][part_name] = sides[side]
      
  eyes['num_eyes'] = n_eyes

  return eyes

def flatten_eyes(eyes):
  """turn the eyes dict into a simple list of mimics objects to use with get_bounding_box()."""
  return [v for k, d in eyes.items() if k != 'num_eyes' for v in d.values()]
