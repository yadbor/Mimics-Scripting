# Dummy mimics.segment module

def HU2GV(HU):
  return HU

class Mask:
  def __init__(self):
    self.crop = None
    self.lo = None
    self.hi = None
    self.n = 0

def create_mask():
  return Mask()

def crop_mask(mask, bbox):
    mask.crop = bbox
    return mask

def threshold(mask, lo, hi):
  mask.lo = lo
  mask.hi = hi
  return mask

def boolean_operations(mask1, mask2, op):
  if op == 'Unite':
    if mask1.lo == None:
      mask1.lo = mask2.lo
    elif mask2.lo != None:
      mask1.lo = min(mask1.lo, mask2.lo)
    if mask1.hi == None:
      mask1.hi = mask2.hi
    elif mask2.hi != None:
      mask1.hi = min(mask1.hi, mask2.hi)
    if mask1.crop == None:
      mask1.crop = mask2.crop
    elif mask2.crop != None:
      mask1.crop = max(mask1.crop, mask2.crop)
    mask1.n += 1 # increment the boolean count
  return mask1
