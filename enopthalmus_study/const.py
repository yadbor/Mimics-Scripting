# mimics constants
# const.py

# Define parameters for planes and crop boxes to trim front of orbit
NUM_PLANES = 10	# Number of planes to divide up orbit
MULT_XY = 1.2   # Factor to expand the XY extents of the crop box
MULT_Z = 1.5    # Factor to expand the Z extent of the first & last crop boxes
SIZE_Y = -20    # Y extent of crop boxes

CLOSING_DIST = 7 # For smartfill 

# Axis names for points stored as [X, Y, Z]
X, Y, Z = 0, 1, 2

# Amount to expand the bounding box of the orbit to make sure everything is covered
BOX_EXPAND = [10, 5, 15]

