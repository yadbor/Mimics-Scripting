#Script: Perform automated segmentation of the orbit and analysis for estimation of enophthalmos in orbital reconstructions - Study with Dieter 2022. 

#Generate a bone object
mask_bone = mimics.segment.create_mask()
mask_bone.name = "Bone Mask"
mimics.segment.threshold(mask_bone, mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071))
part_bone = mimics.segment.calculate_part(mask=mask_bone, quality='High')
part_bone.name = "Bone"


#user defined globe
#sphere1 = mimics.analyze.indicate_sphere(message="Indicate the globe using 3pts on the axial view", show_message_box=True, confirm=True, title=None)
imag = mimics.data.images.get_active()
globe_gv = imag.get_grey_value(point_coordinates= sphere1.center)
print(f"Hu Value at Sphere Center {sphere1.center} is: {mimics.segment.GV2HU(globe_gv)} Hu")
print()


#create globe - dynamic region grow works better. Use that until can solve this. 
#mask_temp_globe = mimics.segment.create_mask()
#mimics.segment.threshold(mask_temp_globe, globe_gv-5, globe_gv+5)
#globe_o = [sphere1.center[0]-sphere1.radius, sphere1.center[1]-sphere1.radius, sphere1.center[2]-sphere1.radius]
#globe_x = [sphere1.radius*2,0,0]
#globe_y = [0,sphere1.radius*2,0]
#globe_z = [0,0,sphere1.radius*2]
#mimics.segment.crop_mask(mask_temp_globe, bounding_box= mimics.BoundingBox3d(globe_o, globe_x, globe_y, globe_z))
#mimics.segment.keep_largest(mask_temp_globe)
#mimics.segment.smart_fill_global(mask_temp_globe, hole_closing_distance=2)
#mimics.segment.smooth_mask(mask_temp_globe)
#mimics.segment.morphology_operations(mask_temp_globe, operation='Erode', number_of_pixels=1,)
#mimics.segment.smooth_mask(mask_temp_globe)
#mimics.segment.keep_largest(mask_temp_globe)




#Indicate Plane from Landmarks
plane1 = mimics.analyze.indicate_plane_points(message="Indicate a plane in the orbit", show_message_box=True, confirm=True, title=None)
plane1.name = "Cut off Plane"
#point1 = mimics.analyze.indicate_point(message='Please indicate point', show_message_box=True, confirm=True, title=None)
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


mimics.segment.activate_edit_mask(mask_orbit_volume, edit_type=None, edit_mode=None)
mimics.dialogs.message_box("Adjust the orbital volume with the Edit Mask tool as needed", title=None, ui_blocking=False)



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

