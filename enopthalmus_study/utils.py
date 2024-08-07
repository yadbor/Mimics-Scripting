# Utilities to help scripting Mimics from python
# 2.0  Rob Day  2024-04-15
<<<<<<< HEAD

from const import * # CONSTANT definitions (* is safe as only consts)
=======
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734

# Specialised iterators. Most are in recent version of itertools
from itertools import tee, zip_longest, islice

<<<<<<< HEAD
import sys # for system functions

=======
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
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

<<<<<<< HEAD
def mask_from_thesholds(name, thresh_lo, thresh_hi, bounding_box=None):
=======
def mask_from_thesholds(name, thresh_lo, thresh_hi):
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
  """Convenience function to create & name a mimics.segment mask with given thresholds."""
  mask = mimics.segment.create_mask()
  mask.name = name
  mimics.segment.threshold(mask, thresh_lo, thresh_hi, bounding_box=bounding_box)
  return mask

<<<<<<< HEAD
def mask_from_material(name, material, bounding_box=None):
=======
def mask_from_material(name, material):
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
  """Convenience function to create & name a mimics.segment mask from a material definition."""
  if material.units == "HU":
    lo_gv = mimics.segment.HU2GV(material.lo)
    hi_gv = mimics.segment.HU2GV(material.hi)  
  else:
    lo_gv = material.lo
    hi_gv = material.hi  
<<<<<<< HEAD
  return mask_from_thesholds(name, lo_gv, hi_gv, bounding_box)
=======
  return mask_from_thesholds(name, lo_gv, hi_gv)
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734

def part_from_mask(name, mask, quality = 'High'):
  """Convenience function to create & name a mimics.segment part from an existing mask."""
  part = mimics.segment.calculate_part(mask = mask, quality = quality)
  part.name = name
  return part
  
def masks_unite(mask1, mask2):
  """Convenience function to do a boolean Unite on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Unite')

def masks_subtract(mask1, mask2):
  """Convenience function to do a boolean Difference on two masks."""
<<<<<<< HEAD
  return mimics.segment.boolean_operations(mask1, mask2, 'Minus')
=======
  return mimics.segment.boolean_operations(mask1, mask2, 'Difference')
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734

def masks_intersect(mask1, mask2):
  """Convenience function to do a boolean Intersect on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Intersect')

def mask_dilate(mask, number_of_pixels = 2, connectivity = 8):
  return mimics.segment.morphology_operations(mask, 
                                              operation='Dilate', 
                                              number_of_pixels=number_of_pixels, 
                                              connectivity=connectivity, 
                                              target_mask_name=None, 
                                              limited_to_mask=None)

# Geometry    
from collections import namedtuple

def spline_geometry(spline):
  """Calculate the min, max, centroid and (optionally) span) of a given mimics.spline."""
  Geometry = namedtuple("Desc", ["max", "min", "mean", "span"])
  max_point = [max(idx) for idx in list(zip(* spline.points))]
  min_point = [min(idx) for idx in list(zip(* spline.points))]
  delta = [(a - b) for a, b in zip(max_point, min_point)]
  mean = [(a + b) / 2.0 for a, b in zip(max_point, min_point)]

  return Geometry(max_point, min_point, mean, delta)

<<<<<<< HEAD
def spline_center(spline):
  # Use the hack that measure.get_bounding_box works for a list of analytical objects, even a list of one.
  bbox = mimics.measure.get_bounding_box([spline])
  mid_pt = bbox_center(bbox)
  return mid_pt

=======
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
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

# Definitions of the Unit Sphere file, created manually.
UNIT_RADIUS = 1
<<<<<<< HEAD
UNIT_SPHERE_FILE = r'D:\Projects & Research\Enophthalmos Study\Mimics-Scripting\enopthalmus_study\unit_sphere.stl'

def sphere_to_mask(s):
  # Load STL as a part
  part = mimics.file.import_stl(filename = UNIT_SPHERE_FILE)
  part_scaled = scale_object(part, s.radius / UNIT_RADIUS)
  mimics.move_object(part_scaled, s.center) # start from (0,0,0) so offset is just the target point
  m = mimics.segment.calculate_mask_from_part(part_scaled)
  m.name = "globe_mask"
  return(m)

  
def bbox_from_intersections(pt_a, pt_b, mult, thickness, depth):
  """Make a bounding box from two points in the XY plane.
  The baseline is mult times wider, and the box has specified thickness (Z) and depth (Y).""" 
  k = (mult - 1)/2.0
  delta = [(a-b) for a, b in zip(pt_a, pt_b)]
  ofs = [ k * i for i in delta]

  p1 = [a + o for a, o in zip(pt_a, ofs)]
  p2 = [b - o for b, o in zip(pt_b, ofs)]

  origin = p1[X], p1[Y], p1[Z] - (thickness/2)
  first_vector=[p2[X]-p1[X], p2[Y]-p1[Y], 0]
  second_vector = [0, depth, 0]
  third_vector = [0, 0, thickness]
  bbox_ab = mimics.BoundingBox3d(origin, first_vector, second_vector, third_vector)
  #mask_ab = mimics.segment.threshold(mask=mimics.segment.create_mask(select_new_mask=False), threshold_min=materials.MIN_GV, threshold_max=materials.MAX_GV, bounding_box=bbox_ab)
  
  return bbox_ab

def bbox_thicken(bbox, extra):
  """Given a bounding box, add extra Z either above (+ve) or below (-ve) the original."""
  if extra < 0:
    bbox.origin = [bbox.origin[X], bbox.origin[Y], bbox.origin[Z] + extra]
    bbox.third_vector = [bbox.third_vector[X], bbox.third_vector[Y], bbox.third_vector[Z] - extra]
  else:
    bbox.third_vector = [bbox.third_vector[X], bbox.third_vector[Y], bbox.third_vector[Z] + extra]
  return bbox

def bbox_to_points(bbox):
  p1 = bbox.origin
  # Add the three vectors
  span = [(a + b + c) 
          for a, b, c 
          in zip(bbox.first_vector, bbox.second_vector, bbox.third_vector)]
  p2 = tuple((p + s) for p, s in zip(p1, span)) # Return a tuple so both points the same
  return (p1, p2)

def bbox_center(bbox):
  # Add the three vectors & divide by 2 to get vectors to centre of bbox
  span = [(a + b + c) / 2
          for a, b, c 
          in zip(bbox.first_vector, bbox.second_vector, bbox.third_vector)]
  # then add the spans to to origin
  mid_pt = tuple((p + s) for p, s in zip(bbox.origin, span))
  return mid_pt

def antero_lateral(bbox, side, delta = 2):
  # The antero-lateral point on left eye is away from origin, on the right it is the origin.
  # This will be on the very edge of the bounding box, so move back toward the centre.
  # delta is the amount in mm to nudge to point into the ROI
  pt = bbox.origin
  if side == 'left':
    pt  = (pt[X] + bbox.first_vector[X] - delta, pt[Y] + delta, pt[Z] + delta)
  elif side == 'right':
    pt  = (pt[X] + delta, pt[Y] + delta, pt[Z] + delta)
  else:
    print(f"Called with side == '{side}', but it must be 'left' or 'right'.")

  return pt

# Find the unit basis vectors for the project.
# These may not be (1,0,0), (0,1,0), (0,0,1) if the project is tilted
def v_hat(v):
  '''Return a normalised unit vector.''' 
  mag = sum([i**2 for i in v])**0.5
  return([i/mag for i in v])

def mimics_image_vectors():
  '''return the three unit vectors that descrivbe this image volume.'''
  active_img = [i for i in mimics.data.images if i.active][0]
  p0 = active_img.get_voxel_center([0, 0, 0])
  d = active_img.get_voxel_buffer().shape
  x = active_img.get_voxel_center([d[0]-1, 0, 0])
  y = active_img.get_voxel_center([0, d[1]-1, 0])
  z = active_img.get_voxel_center([0, 0, d[2]-1])
  if 'numpy' in sys.modules:
    # use the faster neater version
    basis = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
    (i, j, k) = [b_i / np.linalg.norm(b_i, ord=1) for b_i in basis]
  else:
    basis = [tuple(b-a for a,b in zip(p0,v)) for v in (x, y, z)]
    (i, j, k) = [v_hat(v) for v in basis]

  return (i,j,k)


# R equivalent functions for lists
# if use numpy use argmin(vals) and argmax(vals)
def which_max(vals):
    return(vals.index(max(vals)))

def which_min(vals):
    return(vals.index(min(vals)))


def labelled_point(prefix = '', name = '', point=None):
   if name != '':
    name = name + '.'  # add a spacer if needed
    
   if point is None:  # deal with an empty imput gracefully
    point = (None, None, None)

   return {prefix + name + axis: val for axis, val in list(zip(('X','Y','Z'), point))}

import time

class Events:
  """A simple event logger"""
  def __init__(self):
    # init a pair of lists to store the event labels and times
    self.labels = ['init']
    self.times = [time.perf_counter()]
  
  def add(self, label = ''):
    self.labels.append(label)
    self.times.append(time.perf_counter())

  def elapsed(self):
    return self.times[-1] - self.times[0]

  def as_lists(self):
    return self.labels, self.times

  def as_dict(self):
    return dict(zip(self.labels, self.times))

=======
UNIT_SPHERE_FILE = 'unit_sphere.stl'

def sphere_to_mask(s):
  # Load STL as a part
  p = mimics.file.import_stl(filename = UNIT_SPHERE_FILE)
  ps = scale_object(p, s.radius / UNIT_RADIUS)
  mimics.move_object(ps, s.center) # start from (0,0,0) so offset is just the target point
  m = mimics.segment.mask_from_part(p)
  m.name = "globe_mask"
  return(m)
>>>>>>> 94527e54fc0373010a82208260b31cf0f3ac4734
