# mimics constants
# const.py

# Axis names for points stored as [X, Y, Z]
X, Y, Z = 0, 1, 2

# Define parameters for planes and crop boxes to trim front of orbit
NUM_PLANES = 12	# Number of planes to divide up orbit
MULT_XY = 2.0   # Factor to expand the XY extents of the crop box
EXTRA_Z = 5     # Amount to extend the Z extent of the first & last crop boxes
EXTRA_Y = 0.5   # Amount to extend the crop box into the orbit to ensure overlap with the rim
SIZE_Y = -20    # Y extent of crop boxes. Negative is Anterior

CLOSING_DIST = 7 # For smartfill 

# Amount to expand the bounding box of the orbit to make sure everything is covered
BOX_EXPAND = [10, 50, 10]
# If there is no orbit apex point make the box arbitrarily deeper posteriorly
POST_EXPAND = 80

