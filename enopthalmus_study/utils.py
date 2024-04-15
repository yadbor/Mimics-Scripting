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
  return create_mask(name, lo_gv, hi_gv)

def mask_to_part(name, mask, quality = 'High'):
  """Convenience function to create & name a mimics.segment part from an existing mask."""
  part = mimics.segment.calculate_part(mask = mask, quality = quality)
  part.name = name
  return part
  
def masks_unite(mask1, mask2):
  """Convenience function to do a boolean Unite on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Unite')

def masks_subtract(mask1, mask2):
  """Convenience function to do a boolean Difference on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Difference')

def masks_intersect(mask1, mask2):
  """Convenience function to do a boolean Intersect on two masks."""
  return mimics.segment.boolean_operations(mask1, mask2, 'Intersect')

# Geometry    
from collections import namedtuple

def spline_geometry(spline):
  """Calculate the min, max, centroid and (otionally) span) of a given mimics.spline."""
  Desc = namedtuple("Desc", ["max", "min", "mean", "span"])
  max_point = [max(idx) for idx in list(zip(* spline.points))]
  min_point = [min(idx) for idx in list(zip(* spline.points))]
  span = [(a - b) for a, b in zip(max_point, min_point)]
  mean = [(a + b) / 2.0 for a, b in zip(max_point, min_point)]

  return Desc(max_point, min_point, mean, span)
