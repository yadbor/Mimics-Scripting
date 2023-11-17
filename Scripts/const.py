# mimics constants
# const.py

# Define some useful "constants"

NUM_PLANES = 10	# Number of planes to divide up orbit

# Hounsfield Units for various substances and Grey Value
# conversions to use with mimics.segment.threshold()

#import mimics # to use HU to GV converion routine
MIN_HU = -1024
MAX_HU = 3072

MIN_BONE_HU = 226
#MIN_BONE_GV = mimics.segment.HU2GV(MIN_BONE_HU)
MAX_BONE_HU = MAX_HU - 1
#MAX_BONE_GV = mimics.segment.HU2GV(MAX_BONE_HU)

MIN_AIR_HU = MIN_HU
#MIN_AIR_GV = mimics.segment.HU2GV(MIN_AIR_HU)
MAX_AIR_HU = -200
#MAX_AIR_GV = mimics.segment.HU2GV(MAX_AIR_HU)

MIN_FAT_HU = MAX_AIR_HU
#MIN_FAT_GV = mimics.segment.HU2GV(MIN_FAT_HU)
MAX_FAT_HU = -24
#MAX_FAT_GV = mimics.segment.HU2GV(MAX_FAT_HU)

MIN_MUSCLE_HU = MAX_FAT_HU
#MIN_MUSCLE_GV = mimics.segment.HU2GV(MIN_MUSCLE_HU)
MAX_MUSCLE_HU = 100
#MAX_MUSCLE_GV = mimics.segment.HU2GV(MAX_MUSCLE_HU)
