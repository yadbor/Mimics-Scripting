#Script: Perform automated segmentation of the orbit and analysis for estimation of enophthalmos in orbital reconstructions - Study with Dieter 2022. 

# Original  Ryan Collier 2022
# This version  Robert Day  2023-11-18

from result_logger import Path, log_to_file

import const # useful 'constants'
# Axis names for points stored as [X, Y, Z]
X, Y, Z = 0, 1, 2

# Define a Material data structure and some helper functions
from utils import Material, material_mask, mask_to_part
# Materials to analyse
materials = {
  'air'    : const.MATL_AIR, 
  'fat'    : const.MATL_FAT,
  'muscle' : const.MATL_MUSCLE
}

# If there is no "Bone Mask" mask the make a bone mask and part
if mimics.data.objects.find("Bone Mask") is None:
  mask_bone = material_mask("Bone Mask", const.MATL_BONE)
  part_bone = mask_to_part("Bone", mask_bone)
  
# TODO: What about creating the Air mask here as well?

#User Inputs
if mimics.data.objects.find("Spline 1", False) == None:
    spline1 = mimics.analyze.indicate_spline(
      message='Indicate spline on the orbital rim, passing through landmarks', 
      show_message_box=True, confirm=False, title=None
    )
if mimics.data.objects.find("Sphere 1", False) == None:
    sphere1 = mimics.analyze.indicate_sphere(
      message="Indicate the globe using 3pts on the axial view", 
      show_message_box=True, confirm=True, title=None
    )
    
# The DICOM co-ordinate system is defined (for a BIPED) as patient LPS. That is:
#    X+ to the Left hand side of the patient
#    Y+ to the Posterior (towards the back)
#    Z+ to the Superior (towards the head)
#
# Mimics appears to use the DICOM LPS co-ordinate system.
#
# Which eye is being analysed? The head is usually scanned in the centre 
# of the CT co-ordinate system, so X=0 is roughly the midline, meaning 
# the Right eye is X < 0 and the Left eye is X > 0
# TODO - a possibly more robust midline is the mean image X co-ordinate
if sphere1.x < 0:
  side = 'right'
elif sphere1.x > 0:
  side = 'left'
else:
  side = 'ambiguous'

#Check if already calculated orbital volume
if mimics.data.objects.find("Orbital Volume", False) is None:
  # Create Orbital Volume from spline1 and sphere1
  
  # Decompose the spline into straight lines between each pair of points in the spline.
  # Use these to intersect with planes, as mimics doesn't have a spline intersect plane function.
  # Use utils.looped_pairwise to join the last point back to the first point
  spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) for a, b in looped_pairwise(spline1.points)]
 
  # Find the x,y,z limits of the spline points (i.e. the max & min for each axis)
  # In the following code * unpacks the iterator (i.e. spline1_points), 
  # then zip() returns the unpacked lists in parallel (in this case the list of X, Y, Z co-ordinates).
  # The list comprehension then maps max() or min() across each of idx = x, y, z
  # and returns a tuple with the max of min for each axis, which are then assigned.
  (max_x, max_y, max_z) = [max(idx) for idx in list(zip(* spline1.points))]
  (min_x, min_y, min_z) = [min(idx) for idx in list(zip(* spline1.points))]
  
  print(f"X is {min_x} to {max_x}")
  print(f"Y is {min_y} to {max_y}")
  print(f"Z is {min_z} to {max_z}")
  
  delta_z = max_z - min_z
  approx_delta_z = round(delta_z, 0)
  spacing_z = approx_delta_z/const.NUM_PLANES
  print(f"Make {const.NUM_PLANES} planes {spacing_z} apart to span {approx_delta_z}")
    
  # Create a list of plane origin points, matching the orbit in the X,Y plane 
  # and spaced in Z from min_z to min_z + approx_delta_z (approximately max_z)
  orig_x =  orbit.center[0]
  orig_y =  orbit.center[1]
  plane_origins = [[orig_x, orig_y, min_z + (z + 1) * spacing_z] for z in range(const.NUM_PLANES - 1)]
  # Create planes using these origins paralle to the X,Y plane (normal is Z+)
  norm_z = [0, 0, 1]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) for orig in plane_origins]

  # Function to create a mimics.BoundingBox3D given two plane intersection points.
  def make_crop_box(pt_up, pt_down):
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ MULT_XY * (pt_up.x - pt_down.x), MULT_XY * (pt_up.y - pt_down.y), 0]
    vector_y = [0, SIZE_Y, 0]    # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # Put the origin above the down intesection point, midway between the z planes
    origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)
  
  # Create a blank list for the boxes
  boxes = []
  # For each plane, find intersections with the spline lines and create a 
  # bounding box based on the intersection points.
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
      boxes.append(make_crop_box(pt_up, pt_down))
    else:
      print(f"did not find intersection for plane {plane}")

  # Adjust the first and last bounding box to ensure full overlap.
  boxes[0].third_vector[2] = -MULT_Z * boxes[0].third_vector[2]  # first goes down
  boxes[-1].third_vector[2] = MULT_Z * boxes[-1].third_vector[2] # last goes up
   
  # Create the combined mask
  united_masks = mimics.segment.create_mask()
  # Loop through each plane again, create a mask, crop it, and Union it together.
  for bbox in boxes:
    # Create a mask, threshold to cover everything
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    # Crop it with BoundingBox for this plane
    mimics.segment.crop_mask(m, bbox)
    # Union with the existing mask
      united_masks = mimics.segment.boolean_operations(  united_masks, m, 'Unite')

    pause 
            
    # Unite with Bone and Air masks
    united_masks = mimics.segment.boolean_operations(united_masks, mask_bone, 'Unite')
    
    mask_air = mask_bone = material_mask("Air Mask", const.MATL_AIR)
    united_masks = mimics.segment.boolean_operations(united_masks, mask_air, 'Unite')
    # Fill the combined masks
    smartfill_mask = mimics.segment.smart_fill_global(united_masks, 7)
    
    # TODO: replace these numbers with CONSTANTS or something
    # Define the extents of the region to analyse, based on which side is being analyzed.
    # This assumes a LPS coordinate system and that the midline is x == 0.
    if   orbit.center[1] < 0: # Right side
      x_lim = min_x - 10
      x_delta = (max_x - min_x) + 20
    elif orbit.center[1] > 0: # Left side
      x_lim = max_x + 10
      x_delta = (min_x - max_x) - 20
    else:                     # On midline
      print(f"ERROR: Can't tell which side to analyse (orbit at { orbit.x,  orbit.y,  orbit.z}")
    box_orig = (x_lim, min_y - 5, min_z -15)
    vector_x = [x_delta, 0, 0]
    vector_y = [0, 80, 0]
    vector_z = [0, 0, (max_z - min_z) + 30]

    orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    # Crop the filled mask and crrate a part
    smartfill_mask = mimics.segment.crop_mask(smartfill_mask, orbit_ROI)
    smartfill_part = mimics.segment.calculate_part(smartfill_mask, 'High')
    
    smartfill_wrap = mimics.tools.wrap(smartfill_part, 0.2, 10, False, True, True)
    wrapped_mask = mimics.segment.calculate_mask_from_part(part_wrap, None)
    
    orbit_vol = mimics.segment.boolean_operations(wrapped_mask, smartfill_mask, 'Minus')
    orbit_vol = mimics.segment.morphology_operations(orbit_vol, 'Erode', 1, 8, None, None)
    orbit_vol = mimics.segment.region_grow(orbit_vol, orbit_vol, sphere1.center, 'Axial', False, True, connectivity='6-connectivity') 
    orbit_vol.name = "Orbital Volume"
    
    #delete all masks
    for n in reversed(range(0,len(mimics.data.objects))):    
        if mimics.data.objects[n].name == "Bone":
            print("Keeping Bone Object")
        elif mimics.data.objects[n].name == "Spline 1":
            print("Keeping Spline Object")
        elif mimics.data.objects[n].name == "Sphere 1":
            print("Keeping Sphere Object")
        elif mimics.data.objects[n].name == "Orbital Volume":
            print("Keeping Orbital Volume Mask")
        else:
            mimics.data.objects.delete(mimics.data.objects[n])

## Have calcuated Orbital Volume, either this time or on a previous run

#convert sphere to mask, manually. 
mimics.dialogs.message_box("Convert the sphere object into a mask and rename the mask 'Globe'. Click okay to continue", title=None, ui_blocking=False)
globe_mask = mimics.data.objects.find("Globe", False)
if mimics.data.objects.find("Globe", False) != None:
    print("Globe Exists!")
    intersect_vol_mask = mimics.segment.boolean_operations(orbit_vol, globe_mask, 'Minus')
    intersect_vol_mask.name = "Intersect Mask"
else:
    print("ERROR: Cant Find Globe Mask")
 
 
#create parts
part_orbital_vol = mimics.segment.calculate_part(mask=intersect_vol_mask, quality='High')

# Put the orbital volume first in the parts list
parts = {}
parts['orbital'] = part_orbital_vol

# Create a list of masks and corresponding parts for each material in the materials dict
masks = {}
for matl in materials.keys():
  masks[matl] = material_mask(matl + ' mask', materials[matl])
  masks[matl] = mimics.segment.boolean_operations(masks[matl], intersect_vol_mask, 'Intersect')
  if masks[matl].number_of_pixels > 0:
    parts[matl] = mask_to_part(matl, masks[matl])


# Discrete object version
mask_air = mimics.segment.boolean_operations(material_mask("Air Mask", const.MATL_AIR), intersect_vol_mask, 'Intersect')
if mask_air.number_of_pixels > 0:
  part_air = mask_to_part("Air", mask_air)

mask_fat = mimics.segment.boolean_operations(material_mask("Fat Mask", const.MATL_FAT), intersect_vol_mask, 'Intersect')
if mask_fat.number_of_pixels > 0:
 part_fat = mask_to_part("Fat", mask_fat)

mask_muscle = mimics.segment.boolean_operations(material_mask("Muscle Mask", const.MATL_muscle), intersect_vol_mask, 'Intersect')
if mask_muscle.number_of_pixels > 0:
 part_muscle = mask_to_part("Muscle", mask_muscle)
print(f"Orbital Volume (excluding globe) = {part_orbital_vol.volume}mm^3")
if mask_air.number_of_pixels != 0:
    print(f"Air Volume = {part_air.volume}mm^3")
if mask_fat.number_of_pixels != 0:
    print(f"Fat Volume = {part_fat.volume}mm^3")
if mask_muscle.number_of_pixels != 0:
    print(f"Muscle Volume = {part_muscle.volume}mm^3")



# CSV logging below
result_log = Path('name_of_project_log.csv') # TODO: set this at top of program, or make it a configurable CONST

study_name = 'find the name of the project?' # TODO: get from the project name somehow
user = 'me' # TODO: pick from a list

from datetime import date # to get the current date
today = date.today().isoformat()

#results = [user, today, part_orbital_vol.volume, part_air.volume, part_fat.volume, part_muscle.volume]
# This method will automatically cope with changes to the measured volumes
headers = ["user", "date", "study"] + [n + '_volume' for n in parts.keys()]
results = [user,    today, study_name] + [p.volume for p in parts.values()]

log_to_file(result_log, headers, results)
