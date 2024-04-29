# mimics constants
# const.py

# Define parameters for planes and crop boxes to trim front of orbit
NUM_PLANES = 10	# Number of planes to divide up orbit
MULT_XY = 2.0   # Factor to expand the XY extents of the crop box
EXTRA_Z = 5     # Amount to extend the Z extent of the first & last crop boxes
SIZE_Y = -40    # Y extent of crop boxes. Negative is Anterior

CLOSING_DIST = 7 # For smartfill 

# Axis names for points stored as [X, Y, Z]
X, Y, Z = 0, 1, 2

# Amount to expand the bounding box of the orbit to make sure everything is covered
BOX_EXPAND = [10, 15, 15]
# Make it deeper posteriorly to cover the whole orbit
POST_EXPAND = 80

