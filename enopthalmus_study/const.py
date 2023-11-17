# mimics constants
# const.py

# Define some useful "constants"

NUM_PLANES = 10	# Number of planes to divide up orbit

# Hounsfield Units for various substances and 
# Grey Value conversions to use with mimics.segment.threshold()

#import mimics # to use HU to GV converion routine

import mimics.segment

MIN_HU = -1024
MAX_HU = 3072

MIN_BONE_HU = 226
MAX_BONE_HU = MAX_HU - 1

MIN_AIR_HU = MIN_HU
MAX_AIR_HU = -200

MIN_FAT_HU = MAX_AIR_HU
MAX_FAT_HU = -24

MIN_MUSCLE_HU = MAX_FAT_HU
MAX_MUSCLE_HU = 100

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

import utils

MATL_AIR = Material(lo = MIN_AIR_GV, hi = MAX_AIR_GV, units = 'GV')
MATL_FAT = Material(lo = MIN_FAT_GV, hi = MAX_FAT_GV, units = 'GV')
MATL_MUSCLE = Material(lo = MIN_MUSCLE_GV, hi = MAX_MUSCLE_GV, units = 'GV')
MATL_BONE = Material(lo = MIN_BONE_GV, hi = MAX_BONE_GV, units = 'GV')

materials = {
  'air' : MATL_AIR, 
	'fat' : MATL_FAT,
	'muscle' : MATL_MUSCLE
}
			 
