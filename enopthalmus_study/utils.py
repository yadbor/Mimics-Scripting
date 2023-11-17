# Utilities to help scripting Mimics from python
# 1.0	Rob Day	2023-11-13

from itertools import tee, zip_longest, islice
# import numpy as np

# Need to define pairwise() itertools for python 3.7 doesn't have it
# Slight variation to return a list() and not a zip object
def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return list(zip(a, b))
	
# Variant that loops to the first value
def circular_pairwise(iterable):
    # circular_pairwise('ABCDEFG') --> AB BC CD DE EF FG GA
    a, b = tee(iterable)
    return zip_longest(a, b, fillvalue=next(b, None))
	
# return a list in n sized groups ( the las tone may be incomplete)
def batched(iterable, n):
	"Batch data into tuples of length n. The last batch may be shorter."
	# batched('ABCDEFG', 3) --> ABC DEF G
	it = iter(iterable)
	while True:
		batch = tuple(islice(it, n))
		if not batch:
			return
		yield batch
	
import mimics # done automatically by Mimics scripting environment

def create_mask(name, thesh_lo, thresh_hi):
	mask = mimics.segment.create_mask()
	mask.name = name
	mimics.segment.threshold(mask, thresh_lo, thresh_hi)
	return mask

def material_mask(name, material):
	if material.units == "HU":
		lo_gv = mimics.segment.HU2GV(material.lo)
		hi_gv = mimics.segment.HU2GV(material.hi)	
	else:
		lo_gv = material.lo
		hi_gv = material.hi	

	return create_mask(name, lo_gv, hi_gv)

def mask_to_part(name, mask, quality = 'High'):
	part = mimics.segment.calculate_part(mask = mask, quality = quality)
	part.name = name
	return part
	
from dataclasses import dataclass
# Define a class for materials to hold the thresholds
@dataclass
class Material:
	lo: int
	hi: int
	units: str = 'HU'
	
	
