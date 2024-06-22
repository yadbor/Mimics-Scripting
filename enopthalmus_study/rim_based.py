# New insights - version based on the rim bounding box

from const import * # Safe to use * as only CONSTANT variables
import utils
import materials

import mimics

# Get the rim extents 
rim_geometry = utils.spline_geometry(rim)
min_pt = rim_geometry.min
max_pt = rim_geometry.max

# The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
#    X+ to the patient's Left
#    Y+ to the patient's Posterior
#    Z+ to the patient's Superior
orbit_expansion = ((10,10), # X_left, X_right
                   (30,80), # Y_ant, Y_post
                   (15,15)  # Z_inf, Z_sup
                   )

# This doesn't work on a Spline, but does on [Spline]
bbox_rim =  mimics.measure.get_bounding_box([rim, point])
# expand the bbox
bbox_orbit = expand_bbox(bbox_rim, orbit_expansion)
# Make & unite cropped masks for bone and air
mask_bone_roi = mimics.segment.crop_mask(mask_bone, bounding_box = bbox_orbit)

# Close up Orbit Walls + Floor
# Take the air, dilate 2 pixels, combine with bone and smart_fill
mask_air_roi = mimics.segment.crop_mask(mask_air, bounding_box = bbox_orbit)
mask_air_roi = mimics.segment.keep_largest(mask_air_roi)
mask_air_roi = utils.mask_dilate(mask_air_roi, number_of_pixels = 2, connectivity = 8)
mask_air_bone = utils.masks_unite(mask_bone_roi, mask_air_roi)
mask_air_bone = mimics.segment.smart_fill_global(mask_air_bone, hole_closing_distance=7)

# Convert the mask to a part, smooth it and wrap to close up small holes etc.
part_air_bone = utils.part_from_mask(mask_air_bone, quality = 'High')
part_air_bone = mimics.tools.smooth(part_air_bone, 0.5, 5, False, False)
part_air_bone = mimics.tools.wrap(part_air_bone, 0.2, 10, False, True, False)

# Convert the repaired part back to a mask
mask_not_orbit = mimics.segment.calculate_mask_from_part(part_air_bone)
# Make a mask to remove all tissue anterior to the orbital rim
mask_anterior = orbital_analysis.make_orbit_mask(rim, globe)
# Add that to the not_orbit mask
mask_not_orbit = utils.masks_unite(mask_not_orbit, mask_anterior)

# Make a mask of potential orbit contencts (all thresholds in the ROI)
mask_orbit = mimics.segment.threshold(
                         mask = mimics.segment.create_mask(), 
                         threshold_min = materials.MIN_GV, 
                         threshold_max = materials.MAX_GV,
                         bounding_box = bbox_orbit
                         )
# Subtract the not_orbit 
mask_orbit = utils.masks_subtract(mask_orbit, mask_not_orbit)
# And select only tht part adjacent to the globe
back_of_globe = globe.center
back_of_globe[Y] = back_of_globe[Y] + globe.radius
mask_orbit_volume = mimics.segment.region_grow(mask_orbit, 
                                               target_mask = None, 
                                               point = back_of_globe, 
                                               slice_type = "Axial",
                                               keep_original_mask = False, 
                                               multiple_layer = True, 
                                               connectivity = '6-connectivity')

def expand_pts(min_pt, max_pt, expand):
    MIN = 0
    MAX = 1
    # reorder the expansion box for easier calculation
    # so they are ordered (min(X, Y, X), max(X, Y, Z))
    exp_dist = [idx for idx in list(zip(* orbit_expansion))]
    # For each axis, either add (max) or subtract (min) the expansion distance
    max_pt = [(a + b) for a, b in zip(max_pt, exp_dist[MAX])]
    min_pt = [(a - b) for a, b in zip(min_pt, exp_dist[MIN])]
    return min_pt, max_pt

def expand_bbox(bbox, expand):
    """Expand a mimics.BoundingBox3D by adding a vector = (X_left, X_right), (Y_ant, Y_post), (Z_inf, Z_sup)."""
    # reorder the expansion box for easier calculation
    # so they are ordered (min(X, Y, X), max(X, Y, Z))
    exp_min, exp_max = [idx for idx in list(zip(* orbit_expansion))]
    bbox.origin = [(a - b) for a, b in zip(bbox.origin, exp_min)]
    bbox.first_vector = [bbox.first_vector[X] + exp_max[X], 
                         bbox.first_vector[Y],
                         bbox.first_vector[Z]]
    bbox.second_vector = [bbox.second_vector[X], 
                         bbox.second_vector[Y] + exp_max[Y],
                         bbox.second_vector[Z]]
    bbox.third_vector = [bbox.third_vector[X], 
                         bbox.third_vector[Y],
                         bbox.third_vector[Z] + exp_max[Z]]
    return bbox
    
def bounding_box_to_points(bbox):
    """Find the extreme points p1 and p2 of a mimics.BoundingBox3D, with p1 at the origin."""
    p1 = bbox.origin
    # Add each vector to the origin in turn to get from P1 to P2
    p2 = [(a + b) for a, b in zip(bbox.origin, bbox.first_vector)]
    p2 = [(a + b) for a, b in zip(p2, bbox.second_vector)]
    p2 = [(a + b) for a, b in zip(p2, bbox.third_vector)]
    return p1, p2

def points_to_bounding_box(p1, p2):
    """Create a mimics.BoundingBox3D with origin at p1 and vectors to p2."""
    delta = [(a - b) for a, b in zip(p2, p1)]
    bbox = mimics.BoundingBox3d(p1, [delta[X], 0, 0], [0, delta[Y], 0], [0, 0, delta[Z]])

    return bbox

def move_seed_point(sphere):
    p1 = mimics.analyze.create_point(point = sphere.center)
    p1.y = p1.y + sphere.radius
    return p1
    