#Script: Perform automated segmentation of the orbit and analysis for estimation of enophthalmos in orbital reconstructions - Study with Dieter 2022. 

#Generate a bone object if it does not already exist
if mimics.data.objects.find("Bone Mask", False) == None:
    mask_bone = mimics.segment.create_mask()
    mask_bone.name = "Bone Mask"
    mimics.segment.threshold(mask_bone, mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071))
    part_bone = mimics.segment.calculate_part(mask=mask_bone, quality='High')
    part_bone.name = "Bone"


#User Inputs
if mimics.data.objects.find("Spline 1", False) == None:
    spline1 = mimics.analyze.indicate_spline(message='Indicate spline on the orbital rim, passing through landmarks', show_message_box=True, confirm=False, title=None)
if mimics.data.objects.find("Sphere 1", False) == None:
    sphere1 = mimics.analyze.indicate_sphere(message="Indicate the globe using 3pts on the axial view", show_message_box=True, confirm=True, title=None)


#Check if already calculated orbital volume
if mimics.data.objects.find("Orbital Volume", False) == None:
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






















    
