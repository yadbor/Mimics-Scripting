#Script: Perform automated segmentation of the orbit and analysis for estimation of enophthalmos in orbital reconstructions - Study with Dieter 2022. 

# Original  Ryan Collier 2022
# This version  Robert Day  2023-11-18

from result_logger import Path, log_to_file

## Only needed for dummy routines
import mimics
from mimics import segment
from mimics import analyze
from mimics import data

# Define parameters for planes and crop boxes to trim front of orbit
NUM_PLANES = 10	# Number of planes to divide up orbit
MULT_XY = 1.2   # Factor to expand the XY extents of the crop box
MULT_Z = 1.5    # Factor to expand the Z extent of the first & last crop boxes
SIZE_Y = -20    # Y extent of crop boxes

CLOSING_DIST = 7 # For smartfill 

# Axis names for points stored as [X, Y, Z]
X, Y, Z = 0, 1, 2

from utils import looped_pairwise

import materials # Contains definitions of all materials
# Define a Material data structure and some helper functions
from utils import Material, material_mask, part_from_mask
# Segment orbital contents into these Materials & measure volume
materials = {
  'air'    : const.MATL_AIR, 
  'fat'    : const.MATL_FAT,
  'muscle' : const.MATL_MUSCLE
}

# If there is no "Bone Mask" mask the make a bone mask and part
if mimics.data.objects.find("Bone Mask") is None:
  mask_bone = material_mask("Bone Mask", const.MATL_BONE)
  part_bone = part_from_mask("Bone", mask_bone)
# else use the exisitng mask

# TODO: What about creating the Air mask here as well?

## This version assumes that the mimcs project has already had the 
## bone mask and at least one eye defined.
## Each eye consists of a spline for the orbital border, a sphere for the globe and
## optionally a point marking the apex of the orbit.

## THe bone mask will be the last mask defined in the project.
## This allows for the basic thresholded bone mask to be refined by, for example, extractng a connected region.
  
## The spline, sphere and point for each eye are identified, based on their X coordinat.
## CT data go from -ve X on the Right to +ve X on the Left, but because the patient might be scanned 

#User Inputs
if mimics.data.objects.find("Spline 1", False) == None:
    orbital_rim = mimics.analyze.indicate_spline(
      message='Mark a spline on the orbital rim, passing through landmarks', 
      show_message_box=True, confirm=False, title=None
    )
if mimics.data.objects.find("Sphere 1", False) == None:
    globe_sphere = mimics.analyze.indicate_sphere(
      message="Indicate the globe using 3pts on the axial view", 
      show_message_box=True, confirm=True, title=None
    )
    
#Check if already calculated orbital volume
if mimics.data.objects.find("Orbital Volume", False) is None:
  # Create Orbital Volume from orbital_rim and globe_sphere
  
  # Decompose the spline into straight lines between each pair of points in the spline.
  # Use these to intersect with planes, as mimics doesn't have a spline intersect plane function.
  # Use utils.looped_pairwise to join the last point back to the first point
  spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) for a, b in looped_pairwise(orbital_rim.points)]
 
  # Find the x,y,z limits of the spline points (i.e. the max & min for each axis)
  # In the following code * unpacks the iterator (here orbital_rim.points), then
  # zip() returns the unpacked lists in parallel (here the [X, Y, Z] components).
  # The list comprehension then maps max() or min() across each of idx = x, y, z
  # and returns a tuple with the max of min for each axis.
  (max_x, max_y, max_z) = [max(idx) for idx in list(zip(* orbital_rim.points))]
  (min_x, min_y, min_z) = [min(idx) for idx in list(zip(* orbital_rim.points))]
  
  print(f"X is {min_x} to {max_x}")
  print(f"Y is {min_y} to {max_y}")
  print(f"Z is {min_z} to {max_z}")
  
  # Save the orbital rim infor for logging
  orbital_rim_info = {}

  delta_z = max_z - min_z
  approx_delta_z = round(delta_z, 0)
  spacing_z = approx_delta_z/const.NUM_PLANES
  
  print(f"Make {const.NUM_PLANES} planes {spacing_z} apart to span {approx_delta_z}")
    
  # Create a list of plane origin points, matching the orbit in the X,Y plane 
  # and spaced in Z from min_z to min_z + approx_delta_z (approximately max_z)
  orig_x =  orbit.center[X]
  orig_y =  orbit.center[Y]
  plane_origins = [[orig_x, orig_y, min_z + (z + 1) * spacing_z] for z in range(const.NUM_PLANES - 1)]
  # Create planes using these origins paralle to the X,Y plane (normal is Z+)
  norm_z = [0, 0, 1]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) for orig in plane_origins]

  # Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
  def make_crop_box(pt_up, pt_down):
    """Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points."""
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ MULT_XY * (pt_up.x - pt_down.x), MULT_XY * (pt_up.y - pt_down.y), 0]
    vector_y = [0, SIZE_Y, 0]    # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # Put the origin above the down intesection point, midway between the z planes
    origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)
  
  # For each plane, find intersections with the spline lines and create a 
  # bounding box based on the intersection points.
  boxes = [] # Start with a blank list for the boxes
  for plane in z_planes:
    pt_up, pt_down = None, None
    for line in spline_lines:
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z]) # TODO: should this be >=
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z]) # TODO: should this be =<
      if (from_above):
        pt_up = mimics_analyze_create_point_as_line_and_plane_intersection(line, plane)
      if (from_below):
        pt_down = mimics_analyze_create_point_as_line_and_plane_intersection(line, plane)

    if pt_up is not None and pt_down is not None:
      boxes.extend(make_crop_box(pt_up, pt_down))
    else:
      print(f"did not find intersection for plane {plane}")

  # Adjust the first and last bounding box to ensure full overlap.
  boxes[0].third_vector[2] = -MULT_Z * boxes[0].third_vector[2]  # first goes down
  boxes[-1].third_vector[2] = MULT_Z * boxes[-1].third_vector[2] # last goes up
   
  # Create an empty combined mask to start
  united_masks = mimics.segment.create_mask()
  # Loop through each box, create a mask, crop it with the box, and Union it together.
  for bbox in boxes:
    # Create a mask, threshold to cover everything
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    # Crop it with BoundingBox for this plane
    mimics.segment.crop_mask(m, bbox)
    # Union with the existing mask
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    
  # united_masks is now everything in front of the orbit rim defined by orbital_rim
  # out to 

  # Unite with Bone and Air masks
  united_masks = mimics.segment.boolean_operations(united_masks, mask_bone, 'Unite')
  
  mask_air = material_mask("Air Mask", const.MATL_AIR)
  united_masks = mimics.segment.boolean_operations(united_masks, mask_air, 'Unite')
  # Fill the combined masks
  smartfill_mask = mimics.segment.smart_fill_global(united_masks, 7)
  
  # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS. That is:
  #    X+ to the patient's Left
  #    Y+ to the patient's Posterior
  #    Z+ to the patient's Superior
  #
  # The head is usually scanned in the centre of the CT co-ordinate system, 
  # so X == 0 is roughly the midline, Right eye is X < 0 and Left eye is X > 0
  # TODO - a possibly more robust midline is the mean image X co-ordinate
  if orbit.center[X] < 0:
    side = 'right'
  elif orbit.center[X] > 0:
    side = 'left'
  else:
    side = 'ambiguous'
  
  orbit.name = side + ' eye'

  # Calculate an expanded bounding box to cover all of the objects of interest.
  # This probably doesn't depend on side, but the origin of the box will.
  # start from beyond the origin by [-10, -5, -15]
  # TODO: replace these numbers with CONSTANTS or pull out to an options file
  box_orig = (min_x - 10, min_y - 5, min_z - 15)
  vector_x = [(max_x - min_x) + 20, 0, 0]
  vector_y = [0, 80, 0]
  vector_z = [0, 0, (max_z - min_z) + 30]
  
  # I think this is always looking for the lateral margin and making a box
  # with either +ve or -ve width to point towards the midline

  orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
  # Crop the filled mask and create a part
  smartfill_mask = mimics.segment.crop_mask(smartfill_mask, orbit_ROI)
  smartfill_part = mimics.segment.calculate_part(smartfill_mask, 'High')
  
  # Here could use code to improve the bone surface, as at
  # https://github.com/Pythonsegmenter/Orbital-floor-maxillar-sinus-reconstruction
  # also at https://gist.github.com/Pythonsegmenter/bf9c0df43c9d6260e35ad4b786faf90c
  
  smartfill_wrap = mimics.tools.wrap(smartfill_part, 0.2, 10, False, True, True)
  wrapped_mask = mimics.segment.calculate_mask_from_part(part_wrap, None)
  
  orbit_vol = mimics.segment.boolean_operations(wrapped_mask, smartfill_mask, 'Minus')
  orbit_vol = mimics.segment.morphology_operations(orbit_vol, 'Erode', 1, 8, None, None)
  orbit_vol = mimics.segment.region_grow(orbit_vol, orbit_vol, globe_sphere.center, 'Axial', False, True, connectivity='6-connectivity') 
  orbit_vol.name = "Orbital Volume"
  
    # Delete all masks (except the Orbital Volume mask...)
  objects_to_keep = ("Bone", "Spline 1", "Sphere 1", "Orbital Volume")
  for m in mimics.data.objects:
    if m.name in objects_to_keep:
      print(f"Keeping {m.name} {m.type}")
    else:
      mimics.data.objects.delete(m)

# End of: if mimics.data.objects.find("Orbital Volume", False) is None
## Have a calcuated Orbital Volume, either this time or on a previous run.

#convert sphere to mask, manually. !??? Surely there's a better way?
mimics.dialogs.message_box("Convert the sphere object into a mask and rename the mask 'Globe'. Click okay to continue", 
                           title=None, 
                           ui_blocking=False)

globe_mask = mimics.data.objects.find("Globe", False)

# DUMMY: create a Globe mask
globe_mask = material_mask("Bone Mask", const.MATL_BONE)

if globe_mask is not None:
    print("Globe Exists!")
    intersect_vol_mask = mimics.segment.boolean_operations(orbit_vol, globe_mask, 'Minus')
    intersect_vol_mask.name = "Intersect Mask"
else:
    print("ERROR: Cant Find Globe Mask")
 
# create parts
part_orbital_vol = mimics.segment.calculate_part(mask=intersect_vol_mask, quality='High')

# Create a list of masks and corresponding parts for each material in the materials dict
# Put the orbital volume first in the parts list
parts = {'orbital': part_orbital_vol}
masks = {'orbital': intersect_vol_mask}
for matl in materials.keys():
  masks[matl] = material_mask(matl + ' mask', materials[matl])
  masks[matl] = mimics.segment.boolean_operations(masks[matl], intersect_vol_mask, 'Intersect')
  if masks[matl].number_of_pixels > 0:
    parts[matl] = part_from_mask(matl, masks[matl])

# # Discrete object version
# mask_air = mimics.segment.boolean_operations(material_mask("Air Mask", const.MATL_AIR), intersect_vol_mask, 'Intersect')
# if mask_air.number_of_pixels > 0:
#   part_air = mask_to_part("Air", mask_air)
# 
# mask_fat = mimics.segment.boolean_operations(material_mask("Fat Mask", const.MATL_FAT), intersect_vol_mask, 'Intersect')
# if mask_fat.number_of_pixels > 0:
#  part_fat = mask_to_part("Fat", mask_fat)
# 
# mask_muscle = mimics.segment.boolean_operations(material_mask("Muscle Mask", const.MATL_MUSCLE), intersect_vol_mask, 'Intersect')
# if mask_muscle.number_of_pixels > 0:
#  part_muscle = mask_to_part("Muscle", mask_muscle)
#  
# print(f"Orbital Volume (excluding globe) = {part_orbital_vol.volume}mm^3")
# if mask_air.number_of_pixels != 0:
#     print(f"Air Volume = {part_air.volume}mm^3")
# if mask_fat.number_of_pixels != 0:
#     print(f"Fat Volume = {part_fat.volume}mm^3")
# if mask_muscle.number_of_pixels != 0:
#     print(f"Muscle Volume = {part_muscle.volume}mm^3")
#
 
# CSV logging below
result_log = Path('name_of_project_log.csv') # TODO: set this at top of program, or make it a configurable CONST

study_name = 'find the name of the project?' # TODO: get from the project name somehow
user = 'me' # TODO: pick from a list

from datetime import date # to get the current date
today = date.today().isoformat()

#results = [user, today, part_orbital_vol.volume, part_air.volume, part_fat.volume, part_muscle.volume]
# This method will automatically cope with changes to the measured volumes
headers = ["user", "date", "study"] + [n + '_volume' for n in parts.keys()] + ["rim_x_min"]
results = [user,    today, study_name] + [p.volume for p in parts.values()]

log_to_file(result_log, headers, results)

["side", "center", "radius"]
[side] + [(s.center, s.radius) in mimics.data.spheres]
# or [side, globe_sphere.center, globe_sphere.radius]