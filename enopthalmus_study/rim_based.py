# New insights - version based on the rim bounding box

from const import * # Safe to use * as only CONSTANS variables
import utils


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
    bbox.first_vector = [bbox.third_vector[X], 
                         bbox.third_vector[Y],
                         bbox.third_vector[Z] + exp_max[Z]]
    return bbox
    
                   





