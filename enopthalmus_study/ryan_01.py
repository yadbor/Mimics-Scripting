


import mimics # needed for systax checker
import utils

#Generate a bone object
mask_bone = mimics.segment.create_mask()
mask_bone.name = "Bone Mask"
mimics.segment.threshold(mask_bone, mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071))
part_bone = mimics.segment.calculate_part(mask=mask_bone, quality='High')
part_bone.name = "Bone"

def ryan_version_1(sphere1, spline1, point1, mask_bone):
    # passed in sphere1, spline1, and point1
    # sphere1 = mimics.data.spheres[0]
    imag = mimics.data.images.get_active()
    globe_gv = imag.get_grey_value(point_coordinates= sphere1.center)
    print(f"Hu Value at Sphere Center {sphere1.center} is: {mimics.segment.GV2HU(globe_gv)} Hu")
    print()

    # Replace Ryan's code to create globe from sphere (but it doesn't get used?)
    mask_temp_globe = utils.sphere_to_mask(sphere1)


    #Indicate Plane from Landmarks
    #plane1 = mimics.analyze.indicate_plane_points(message="Indicate a plane in the orbit", show_message_box=True, confirm=True, title=None)
    spline1 = mimics.data.splines[0]
    plane1 = mimics.analyze.create_plane_fit_to_spline(spline=spline1)
    plane1.name = "Cut off Plane"

    ###### WHAT IS POINT 1?? ######
    #point1 = mimics.analyze.indicate_point(message='Please indicate point', show_message_box=True, confirm=True, title=None)
    #point1 = mimics.data.points[0]
    line1 = mimics.analyze.create_line(plane1.origin, point1, name=None, color=None)

    #Crop Bone to Orbit Volume based on Plane and line distance
    w = plane1.width
    h = plane1.height+10
    d = line1.length
    mimics.data.objects.delete(line1)
    corner_x = plane1.origin[0]+plane1.x_axis[0]*w/2+plane1.y_axis[0]*h/2
    corner_y = plane1.origin[1]+plane1.x_axis[1]*w/2+plane1.y_axis[1]*h/2
    corner_z = plane1.origin[2]+plane1.x_axis[2]*w/2+plane1.y_axis[2]*h/2
    x_vector = (plane1.x_axis[0]*(-w), plane1.x_axis[1]*(-w), plane1.x_axis[2]*(-w))
    y_vector = (plane1.y_axis[0]*(-h), plane1.y_axis[1]*(-h), plane1.y_axis[2]*(-h))
    z_vector = (plane1.z_axis[0]*d, plane1.z_axis[1]*d, plane1.z_axis[2]*d)
    corner = (corner_x, corner_y, corner_z)
    mimics.segment.crop_mask(mask_bone, bounding_box= mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector))

    #Repair Orbit Walls + Floor
    mask_temp_air = mimics.segment.create_mask()
    mimics.segment.threshold(mask_temp_air, mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
    mimics.segment.crop_mask(mask_temp_air, bounding_box= mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector))
    mimics.segment.keep_largest(mask_temp_air)
    mask_temp_morph = mimics.segment.morphology_operations(mask_temp_air, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=None, limited_to_mask=None)
    mask_bone_boolean = mimics.segment.boolean_operations(mask_temp_morph, mask_bone, 'Unite')
    mask_bone_filled = mimics.segment.smart_fill_global(mask=mask_bone_boolean, hole_closing_distance=7)
    part_bone_repaired = mimics.segment.calculate_part(mask_bone_filled, quality='High')
    mimics.data.objects.delete(mask_temp_air)
    mimics.data.objects.delete(mask_temp_morph)
    mimics.data.objects.delete(mask_bone_boolean)
    mimics.data.objects.delete(mask_bone_filled)
    part_bone_smoothed = mimics.tools.smooth(part_bone_repaired, 0.5, 5, False, False)
    part_bone_wrapped = mimics.tools.wrap(part_bone_smoothed, 0.2, 10, False, True, False)
    part_bone_wrapped.name = "Smoothed & Wrapped Orbit"
    subtraction_mask = mimics.segment.calculate_mask_from_part(part_bone_wrapped, None)
    subtraction_mask.name = "Subtraction Mask"


    #Create Orbit Volume Mask
    mask_temp_orbit = mimics.segment.create_mask()
    mimics.segment.threshold(mask_temp_orbit, mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))
    mimics.segment.crop_mask(mask_temp_orbit, bounding_box= mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector)) #change this to landmarks
    mask_intersect = mimics.segment.boolean_operations(mask_a= mask_temp_orbit, mask_b= subtraction_mask, operation='Minus')
    mask_orbit_volume = mimics.segment.region_grow(mask_intersect, None, sphere1.center, "Axial", keep_original_mask=False, multiple_layer=True, connectivity='6-connectivity')
    mask_orbit_volume.name = "Orbit Volume"
    mimics.data.objects.delete(subtraction_mask)
    mimics.data.objects.delete(mask_temp_orbit)
    mimics.data.objects.delete(mask_intersect)


    # Adjust mask_orbit_volume as needed

    
    #Create Air Mask for Calcs
    mask_temp_air = mimics.segment.create_mask()
    mimics.segment.threshold(mask_temp_air, mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
    mask_air = mimics.segment.boolean_operations(mask_a= mask_temp_air, mask_b= mask_orbit_volume, operation='Intersect')
    mimics.data.objects.delete(mask_temp_air)
    mask_air.name = "Air Mask"

    #Create Fat Mask for Calcs
    mask_temp_fat = mimics.segment.create_mask()
    mimics.segment.threshold(mask_temp_fat, mimics.segment.HU2GV(-200), mimics.segment.HU2GV(-24))
    mask_fat = mimics.segment.boolean_operations(mask_a= mask_temp_fat, mask_b= mask_orbit_volume, operation='Intersect')
    mimics.data.objects.delete(mask_temp_fat)
    mask_fat.name = "Fat Mask"

    #Create Air Mask for Calcs
    mask_temp_muscle = mimics.segment.create_mask()
    mimics.segment.threshold(mask_temp_muscle, mimics.segment.HU2GV(-24), mimics.segment.HU2GV(100))
    mask_muscle = mimics.segment.boolean_operations(mask_a= mask_temp_muscle, mask_b= mask_orbit_volume, operation='Intersect')
    mimics.data.objects.delete(mask_temp_muscle)
    mask_muscle.name = "Muscle Mask"


    return (mask_orbit_volume, mask_air, mask_fat, mask_muscle)

def ryan_version_2(spline1, sphere1):
    print(len(spline1.points))
    spline_lines = {}
    for n in range(len(spline1.points)):   
        if n == len(spline1.points)-1:
            spline_lines[n] = mimics.analyze.create_line(spline1.points[n], spline1.points[0], name=None, color=None)
        else:
            spline_lines[n] = mimics.analyze.create_line(spline1.points[n], spline1.points[n+1], name=None, color=None)   
    
    #loop to find max value of z in spline
    for n in range(len(spline1.points)) :
        if n == 0:
            max_x = spline1.points[n][0]  
            min_x = spline1.points[n][0]
            max_y = spline1.points[n][1]  
            min_y = spline1.points[n][1]
            max_z = spline1.points[n][2]  
            min_z = spline1.points[n][2]  
        else:
            if max_x < spline1. points[n][0]:       
                max_x = spline1.points[n][0]   
            if min_x > spline1. points[n][0]:       
                min_x = spline1.points[n][0] 
            if max_y < spline1. points[n][1]:       
                max_y = spline1.points[n][1]   
            if min_y > spline1. points[n][1]:       
                min_y = spline1.points[n][1] 
            if max_z < spline1. points[n][2]:       
                max_z = spline1.points[n][2]   
            if min_z > spline1. points[n][2]:       
                min_z = spline1.points[n][2]   
    
    print(max_z)
    print(min_z)
    delta_z = round(max_z - min_z, 0)
    print(delta_z)
    spacing_z = delta_z/10
    
    z_planes = {}
    for n in range(9):
        orig = [sphere1.center[0], sphere1.center[1], min_z + (n+1)*spacing_z]
        z_planes[n] = mimics.analyze.create_plane_origin_and_normal(orig, [0,0,1], name=None, color=None)
        
        
    intersect_points = {}
    for n in range(len(z_planes)):
        for m in range(len(spline_lines)):
            if spline_lines[m].point1[2] > z_planes[n].origin[2] and spline_lines[m].point2[2] < z_planes[n].origin[2]:          
                intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n], name=None, color=None)
            if spline_lines[m].point1[2] < z_planes[n].origin[2] and spline_lines[m].point2[2] > z_planes[n].origin[2]:
                intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n], name=None, color=None) 
    
    crop_masks = {} 
    for n in range(len(intersect_points)):   
        if n % 2 == 1: 
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
    
    
    # pause 
            
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

