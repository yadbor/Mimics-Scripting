# Utilities to help scripting Mimics from python
# 2.0  Rob Day  2024-04-15

# Specialised iterators. Most are in recent version of itertools
from itertools import tee, zip_longest, islice
import numpy as np 

# Need to define pairwise() as itertools for python 3.7 doesn't have it.
# Slight variation to return a list() and not a zip object.
def pairwise(iterable):
  """Return consecutive overlapping pairs from a list. pairwise('ABCDEFG') --> AB BC CD DE EF FG."""
  # pairwise('ABCDEFG') --> AB BC CD DE EF FG
  a, b = tee(iterable)
  next(b, None)
  return list(zip(a, b))
  
# Variant of pairwise that loops to the first value (like a ring buffer).
def looped_pairwise(iterable):
  """Return consecutive overlapping pairs from a list. The final pair is (item, first)."""
  # looped_pairwise('ABCDEFG') --> AB BC CD DE EF FG GA
  a, b = tee(iterable)
  return zip_longest(a, b, fillvalue=next(b, None))
  
# Return a list in n sized groups (last group may be incomplete).
def batched(iterable, n):
  """Batch data into tuples of length n. The last batch may be shorter."""
  # batched('ABCDEFG', 3) --> ABC DEF G
  it = iter(iterable)
  while True:
    batch = tuple(islice(it, n))
    if not batch:
      return
    yield batch
  
# Helper functions for Mimics objects
import mimics # done automatically by Mimics scripting environment

def mask_from_thesholds(name, thresh_lo, thresh_hi):
  """Convenience function to create & name a mimics.segment mask with given thresholds."""
  mask = mimics.segment.create_mask()
  mask.name = name
  mimics.segment.threshold(mask, thresh_lo, thresh_hi)
  return mask

def mask_from_material(name, material):
  """Convenience function to create & name a mimics.segment mask from a material definition."""
  if material.units == "HU":
    lo_gv = mimics.segment.HU2GV(material.lo)
    hi_gv = mimics.segment.HU2GV(material.hi)  
  else:
    lo_gv = material.lo
    hi_gv = material.hi  
  return mask_from_thesholds(name, lo_gv, hi_gv)

def part_from_mask(name, mask, quality = 'High'):
  """Convenience function to create & name a mimics.segment part from an existing mask."""
  part = mimics.segment.calculate_part(mask = mask, quality = quality)
  part.name = name
  return part
  
def unite(mask1, mask2):
  """Convenience function to do a boolean Unite on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Unite')

def minus(mask1, mask2):
  """Convenience function to do a boolean Difference on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Difference')

def intersect(mask1, mask2):
  """Convenience function to do a boolean Intersect on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Intersect')

def mask_dilate(mask, number_of_pixels = 2, connectivity = 8):
  return mimics.segment.morphology_operations(mask, 
                                              operation='Dilate', 
                                              number_of_pixels=number_of_pixels, 
                                              connectivity=connectivity, 
                                              target_mask_name=None, 
                                              limited_to_mask=None)

# Functions for assiging geometry (sphere, spline & point) to the L or R eye

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

def flatten_eyes(eyes):
  """turn the eyes dict into a simple list of mimics objects to use with get_bounding_box()."""
  return [v for k, d in eyes.items() if k != 'num_eyes' for v in d.values()]


# Functions for creating and manipulating mimcs.BoundingBox3D objects
# Define a deault set of basis vectors to use if get_basis_vector() not called

DEFAULT_BASIS = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

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
    bbox = mimics.BoundingBox3d(p1, span[0], span[1], span[2])
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


# Geometry    
from collections import namedtuple

def spline_geometry(spline):
  """Calculate the min, max, centroid and (optionally) span) of a given mimics.spline."""
  Geometry = namedtuple("Desc", ["max", "min", "mean", "span"])
  max_point = [max(idx) for idx in list(zip(* spline.points))]
  min_point = [min(idx) for idx in list(zip(* spline.points))]
  delta = [(a - b)       for a, b in zip(max_point, min_point)]
  mean  = [(a + b) / 2.0 for a, b in zip(max_point, min_point)]

  return Geometry(max_point, min_point, mean, delta)

# There is no API to make a mask from the globe Spheres, so need a workaround
# Before starting the study:
# 1. Create a unit sphere at the origin
# 2. Create a Part from that using the GUI
# 3. Save the part as an stl file
#
# To make a mask from the globe:
# 4. Load the stl as a part (it will be at the origin)
# 5. Scale it to the radius of the globe
# 6. Move it to the center of the globe 
# 
# Could scale and translate in one transform, but this is cleaner
# There is no API call to scale an object, so roll our own from scratch
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

# Definitions of the Unit Sphere file, created manually.
UNIT_RADIUS = 1
UNIT_SPHERE_FILE = 'unit_sphere.stl'

def sphere_to_mask(s):
  # Load STL as a part
  p = mimics.file.import_stl(filename = UNIT_SPHERE_FILE)
  ps = scale_object(p, s.radius / UNIT_RADIUS)
  mimics.move_object(ps, s.center) # start from (0,0,0) so offset is just the target point
  m = mimics.segment.mask_from_part(p)
  m.name = "globe_mask"
  return(m)
