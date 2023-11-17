#Script: Perform automated segmentation of the orbit and analysis for estimation of enophthalmos in orbital reconstructions - Study with Dieter 2022. 

# Original	Ryan Collier 2022
# This version	Robert Day	2023-11-14

import const # useful 'constants'
# Axis names for points stored as [X, Y, Z]
X = 0
Y = 1
Z = 2

# Define a Material data structure and some helper functions
from utils import Material, material_mask, mask_to_part
# Materials to analyse
materials = {
  'air' : const.MATL_AIR, 
	'fat' : const.MATL_FAT,
	'muscle' : const.MATL_MUSCLE
}

# If there is no "Bone Mark" mask the make a bone mask and part
if mimics.data.objects.find("Bone Mask") is None:
	mask_bone = material_mask("Bone Mask", const.MATL_BONE)
	part_bone = mask_to_part("Bone", mask_bone)

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
		
# The DICOM co-ordinate system is defined (for BIPED) as patient LPS. That is:
#    X+ to the Left hand side of the patient
#    Y+ to the Posterior (towards the back)
#    Z+ to the Superior (towards the head)
# Mimics appears to use the DICOM LPS co-ordinate system.

# Which eye is being processed? The head is usually scanned in the 
# centre of the CT co-ordinate system, so X=0 is roughly the midline
# meaning the Right eye is X < 0 and the Left eye is X > 0
# TODO - a possibly more robust midline is the mean image X co-ordinates
if sphere1.x < 0:
	side = 'right'
elif sphere1.x > 0:
	side = 'left'
else:
	side = 'ambiguous'

#Check if already calculated orbital volume
if mimics.data.objects.find("Orbital Volume", False) == None:
    print(len(spline1.points))
	
	# Decompose the spline into straight lines between each pair of points in the spline.
    # Use these to intersect with planes, as mimics doesn't have a spline intersect plane function.
	# Use utils.circular_pairwise to join the last point back to the first point
	spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) for a, b in circular_pairwise(spline1.points)]
 
    # Find the x,y,z limits of the spline points (i.e. the max & min for each axis)
	# In the following code * unpacks the iterator (i.e. spline1_points), 
	# then zip() returns the unpacked lists in parallel (in this case the list of X, Y, Z co-ordinates).
	# The list comprehension then maps the max() or min() function across each of idx = x, y, z.
    # The comprehension retrns a tuple with the max of min for each axis, which are then assigned.
	(max_x, max_y, max_z) = [max(idx) for idx in list(zip(* spline1.points))]
	(min_x, min_y, min_z) = [min(idx) for idx in list(zip(* spline1.points))]
	
	print(f"X is {min_x} to {max_x}")
	print(f"Y is {min_y} to {max_y}")
	print(f"Z is {min_z} to {max_z}")
	
    # print(max_z)
    # print(min_z)
    delta_z = round(max_z - min_z, 0)
    print(delta_z)
    spacing_z = delta_z/const.NUM_PLANES
    
	# Create a list of plane origin points, each aligned to the orbit in the X,Y plane 
	# and spaced in Z from min_z to min_z + delta_z (which is close to max_z)
	orig_x = sphere1.center[0]
	orig_y = sphere1.center[1]
	plane_origins = [[orig_x, orig_y, min_z + (z + 1) * spacing_z] for z in range(NUM_PLANES - 1)]
	# Create planes using these origins and all in the X,Y plane (normal is Z+)
	norm_z = [0, 0, 1]
	z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) for orig in plane_origins]
	
	# Create a list of intersection points between each plane and the spline line segments.
	# Step through the planes from low Z to high Z, and for each plane find the two intersection points with the spline.
	# One intersection will be the spline going up through the plane, and the other going down.
  intersect_points = [] # List to store intersection points.
	for n, plane in enumerate(z_planes):
		# Go around the spline lines and for each line segment check if the y values span the plane.
		# If yes, then calculate the intercept. Going around the whole line gets both the upwards
		# and downwards intersections, one following the other so each pair of points is for the same plane
		for line in spline_lines:
			if ( (line.point1[Y] > plane.origin[Y] and line.point2[Y] < plane.origin[Y]) # line crosses from above to below 
				  or
				 (line.point1[Y] < plane.origin[Y] and line.point2[Y] > plane.origin[Y]) # line crosses from below to above
				):
				intersect_points.append(mimics.analyze.create_point_as_line_and_plane_intersection(line, plane))
	
    # intersect_points = {} # here len(intersect_points) == 0, so using len()+1 below starts indices from 1
    # for n in range(len(z_planes)):
        # for m in range(len(spline_lines)):
            # if spline_lines[m].point1[2] > z_planes[n].origin[2] and spline_lines[m].point2[2] < z_planes[n].origin[2]:          
                # intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n], name=None, color=None)
            # if spline_lines[m].point1[2] < z_planes[n].origin[2] and spline_lines[m].point2[2] > z_planes[n].origin[2]:
                # intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n], name=None, color=None) 
                
  crop_masks = []
  point_pairs = list(batched(intersect_points_RED, 2))
	for i, (pt_up, pt_down) in enumerate(point_pairs):
		vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0] # why 1.2???
		vector_y = [0, -20, 0] # -Y is anterior, so point towards front of face
		if i == 0: # first point
		  ofs = -1.5
		elif i == len(point_pairs):
		  ofs = +1.5
		else
		  ofs = 1
		vector_z = [0, 0, ofs * spacing_z]
		box_orig = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)  
    crop_box = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    # create the mask, threshold and crop it, then append to the list
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV))
    crop_masks_RED.append(mimics.segment.crop_mask(m, crop_box))

# crop_masks_RED now has len(intesect_points) // 2 masks
		
# Ryan version	
    crop_masks = {} 
    for n in range(len(intersect_points)):   
        if n % 2 == 1: 
			# If I read this correctly, this makes a line from points (1,2), (3, 4) etc.
			# What happened to point zero? ANS: we lost it in creating the list, as used [n+1] and n starts from 0
            vector_x = [(intersect_points[n].x-intersect_points[n+1].x)*1.2,(intersect_points[n].y-intersect_points[n+1].y)*1.2,0]
            vector_y = [0,-20,0]    
            if n == 1:
                crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
                mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
                box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z+spacing_z/2)         
                vector_z = [0,0,-1.5*spacing_z]        
                mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))       
            elif n == len(intersect_points)-1:
                crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
                mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
                box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)         
                vector_z = [0,0,1.5*spacing_z]        
                mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
            else:
                box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)          
                vector_z = [0,0,spacing_z] 
                crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
                mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
                mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
        
    #boolean masks
    for n in range(len(crop_masks)+1):
        if n == 0:
            boolean_masks = mimics.segment.create_mask()
        else:
            boolean_masks = mimics.segment.boolean_operations(boolean_masks, crop_masks[n], 'Unite')
    
    
    pause 
            
    #boolean w air
    boolean_masks = mimics.segment.boolean_operations(boolean_masks, mask_bone, 'Unite')
    mask_air = mimics.segment.create_mask()
    mask_air = mimics.segment.threshold(mask_air, mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
    boolean_masks = mimics.segment.boolean_operations(boolean_masks, mask_air, 'Unite')            
    smartfill_mask = mimics.segment.smart_fill_global(boolean_masks, 7)
    
        
    if max_x > 0 and min_x > 0: #right orbit
        box_orig = (min_x-10, min_y-5, min_z-15)          
        vector_x = [max_x - min_x+20,0,0]
        vector_y = [0,80,0]
        vector_z = [0,0,max_z - min_z + 30] 
    elif max_x < 0 and min_x < 0: #left orbit
        box_orig = (max_x+10, min_y-5, min_z-15)          
        vector_x = [min_x-max_x-20,0,0]
        vector_y = [0,80,0]
        vector_z = [0,0,max_z - min_z + 30]
    else:
        print("ERROR: Orbits not greater than or less than zero coord in x axis")    
    
    
    smartfill_mask = mimics.segment.crop_mask(smartfill_mask, mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
    part_temp = mimics.segment.calculate_part(smartfill_mask, 'High')
    part_wrap = mimics.tools.wrap(part_temp, 0.2, 10, False, True, True)
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

# Create a list of masks and corresponding parts for each material in the materials dict
masks = {}
parts = {}
for matl in materials.keys():
	masks[matl] = material_mask(matl + ' mask', materials[matl])
	masks[matl] = mimics.segment.boolean_operations(masks[matl], intersect_vol_mask, 'Intersect')
	if masks[matl].number_of_pixels > 0:
		parts[matl] = mask_to_part(matl, masks[matl])



# mask_air = mimics.segment.boolean_operations(material_mask("Air Mask", const.MATL_AIR), intersect_vol_mask, 'Intersect')
# if mask_air.number_of_pixels > 0:
	# part_air = mask_to_part("Air", mask_air)

# mask_fat = mimics.segment.boolean_operations(material_mask("Fat Mask", const.MATL_FAT), intersect_vol_mask, 'Intersect')
# if mask_fat.number_of_pixels > 0:
	# part_fat = mask_to_part("Fat", mask_fat)

# mask_muscle = mimics.segment.boolean_operations(material_mask("muscle Mask", const.MATL_muscle), intersect_vol_mask, 'Intersect')
# if mask_muscle.number_of_pixels > 0:
	# part_muscle = mask_to_part("muscle", mask_muscle)

	
mask_temp = mimics.segment.create_mask()
mimics.segment.threshold(mask_temp, mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
mask_air = mimics.segment.boolean_operations(mask_temp, intersect_vol_mask, 'Intersect')
mimics.data.objects.delete(mask_temp)
mask_air.name = "Air Mask"
if mask_air.number_of_pixels != 0: 
    part_air = mimics.segment.calculate_part(mask=mask_air, quality='High')
    part_air.name = "Air"

mask_temp = mimics.segment.create_mask() 
mimics.segment.threshold(mask_temp, mimics.segment.HU2GV(-200), mimics.segment.HU2GV(-24))
mask_fat = mimics.segment.boolean_operations(mask_temp, intersect_vol_mask, 'Intersect')
mimics.data.objects.delete(mask_temp)
mask_fat.name = "Fat Mask"
if mask_fat.number_of_pixels != 0: 
    part_fat = mimics.segment.calculate_part(mask=mask_fat, quality='High')
    part_fat.name = "Fat"
  
mask_temp = mimics.segment.create_mask() 
mimics.segment.threshold(mask_temp, mimics.segment.HU2GV(-24), mimics.segment.HU2GV(100))
mask_muscle = mimics.segment.boolean_operations(mask_temp, intersect_vol_mask, 'Intersect')
mimics.data.objects.delete(mask_temp)
mask_muscle.name = "Muscle Mask"
if mask_muscle.number_of_pixels != 0: 
    part_muscle = mimics.segment.calculate_part(mask=mask_muscle, quality='High')
    part_muscle.name = "Muscle"
    

print(f"Orbital Volume (excluding globe) = {part_orbital_vol.volume}mm^3")
if mask_air.number_of_pixels != 0:
    print(f"Air Volume = {part_air.volume}mm^3")
if mask_fat.number_of_pixels != 0:
    print(f"Fat Volume = {part_fat.volume}mm^3")
if mask_muscle.number_of_pixels != 0:
    print(f"Muscle Volume = {part_muscle.volume}mm^3")


### RED CSV logging below
# Get the user from a list at start of program
from datetime import date # to get the current date
today = date.today().isoformat()
results = [user, today, part_orbital_vol.volume, part_air.volume, part_fat.volume, part_muscle.volume]

# Below from https://stackoverflow.com/questions/71275961/append-new-line-to-csv-file
	
##Each time you use the write w mode, it overwrites the file, deleting the old contents in the process.
##
##It might make more sense to first check if the file exists, and write the header row only if it doesn't. 
##Here's one way to do that, and it shows the behaviour required in your question:

import csv
from pathlib import Path

FILE_PATH = Path('log.csv')

# If the log file does not exist create it with a header row
if not FILE_PATH.exists():
    with open(FILE_PATH, 'w', newline='') as log_csv:
        # Start a new blank log with column headings in the first row
        log_csv_write = csv.writer(log_csv)
		log_csv_write.writerow(["user", "date", "study", "Orbital Volume", "Air Volume", "Fat Volume", "Muscle Volume"])

# Append the current results (wtih autocloses at exit)
with open(FILE_PATH, 'a', newline='') as log_csv:
    log_csv_append = csv.writer(log_csv, dialect = 'excel')
    log_csv_append.writerow(results)


















    
