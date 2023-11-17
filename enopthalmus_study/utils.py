# Utilities to help scripting Mimics from python
# 1.0  Rob Day  2023-11-13

# Specialised iterators. Most are in recent version of itertools
from itertools import tee, zip_longest, islice

# Need to define pairwise() itertools for python 3.7 doesn't have it.
# Slight variation to return a list() and not a zip object.
def pairwise(iterable):
  """Return consecutive overlapping pairs from a list. pairwise('ABCDEFG') --> AB BC CD DE EF FG."""
  # pairwise('ABCDEFG') --> AB BC CD DE EF FG
  a, b = tee(iterable)
  next(b, None)
  return list(zip(a, b))
  
# Variant of pairwise that loops to the first value (like a ring buffer).
def circular_pairwise(iterable):
  """Return consecutive overlapping pairs from a list. The final pair is (item, first)."""
  # circular_pairwise('ABCDEFG') --> AB BC CD DE EF FG GA
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

def create_mask(name, thesh_lo, thresh_hi):
  """Convenience function to create & name a mimics.segment mask with given thresholds."""
  mask = mimics.segment.create_mask()
  mask.name = name
  mimics.segment.threshold(mask, thresh_lo, thresh_hi)
  return mask

def material_mask(name, material):
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
  
# A dataclass is like a record or c struct, with named attributes
from dataclasses import dataclass
# Define a class to hold the thresholds for each material
@dataclass
class Material:
  """Class to hold a material definition, given by hi & lo thresholds. units can be 'HU' (default) or 'GV'"""
  lo: int
  hi: int
  units: str = 'HU' # can be HU or GR for Grey Values
