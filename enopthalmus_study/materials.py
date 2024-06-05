# mimics material constants
# materials.py

# Define some useful "constants"

# Hounsfield Units (HU) for various substances and their Grey Value (GV)
# conversions to use with mimics.segment.threshold(). Mimics routines all use GV.

# These have all been chosen by Ryan for the enopthalmus study based on input
# from surgeons about accuracy of orbital contents segmentation.

import mimics.segment # For Dummy HU to GV converion routine.

MIN_HU = -1024 # The minimum possible HU
MAX_HU =  3071 # The maximum possible HU

MIN_BONE_HU = 148 # Changed from 226 to match Mimics range for 'Spongial Bone' + 'Cortical Bone'
MAX_BONE_HU = MAX_HU

MIN_AIR_HU = MIN_HU
MAX_AIR_HU = -200

MIN_FAT_HU = MAX_AIR_HU
MAX_FAT_HU = -24

MIN_MUSCLE_HU = MAX_FAT_HU
MAX_MUSCLE_HU = 100

# Pre-calculate GV equivalents of all material thresholds
MIN_GV = mimics.segment.HU2GV(MIN_HU)
MAX_GV = mimics.segment.HU2GV(MAX_HU)

MIN_BONE_GV = mimics.segment.HU2GV(MIN_BONE_HU)
MAX_BONE_GV = mimics.segment.HU2GV(MAX_BONE_HU)

MIN_AIR_GV = mimics.segment.HU2GV(MIN_AIR_HU)
MAX_AIR_GV = mimics.segment.HU2GV(MAX_AIR_HU)

MIN_FAT_GV = mimics.segment.HU2GV(MIN_FAT_HU)
MAX_FAT_GV = mimics.segment.HU2GV(MAX_FAT_HU)

MIN_MUSCLE_GV = mimics.segment.HU2GV(MIN_MUSCLE_HU)
MAX_MUSCLE_GV = mimics.segment.HU2GV(MAX_MUSCLE_HU)

# A dataclass is like a record or c struct, with named attributes
from dataclasses import dataclass
# Define a class to hold the thresholds for each material
@dataclass
class Material:
  """Class to hold a material definition, given by hi & lo thresholds. units can be 'HU' (default) or 'GV'"""
  lo: int
  hi: int
  units: str = 'HU' # can be HU or GR for Grey Values

MATL_AIR = Material(lo = MIN_AIR_GV, hi = MAX_AIR_GV, units = 'GV')
MATL_FAT = Material(lo = MIN_FAT_GV, hi = MAX_FAT_GV, units = 'GV')
MATL_MUSCLE = Material(lo = MIN_MUSCLE_GV, hi = MAX_MUSCLE_GV, units = 'GV')
MATL_BONE = Material(lo = MIN_BONE_GV, hi = MAX_BONE_GV, units = 'GV')

