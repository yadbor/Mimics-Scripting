# Utilities to help scripting Mimics from python
# 2.0  Rob Day  2024-04-15

# Specialised iterators. Most are in recent version of itertools
from itertools import tee, zip_longest, islice

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
UNIT_SPHERE_FILE = 'unit_sphere.stl'

def sphere_to_mask(s):
  # Load STL as a part
  p = mimics.file.import_stl(filename = UNIT_SPHERE_FILE)
  ps = scale_object(p, s.radius / UNIT_RADIUS)
  mimics.move_object(ps, s.center) # start from (0,0,0) so offset is just the target point
  m = mimics.segment.mask_from_part(p)
  m.name = "globe_mask"
  return(m)
