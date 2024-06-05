plane1 = mimics.analyze.indicate_plane_points(message="Indicate a plane in the orbit", show_message_box=True, confirm=False, title=None)
plane1.name = "Cut off Plane"


mask_bone = mimics.segment.create_mask()
mask_bone.name = "Bone Mask"
mimics.segment.threshold(mask_bone, mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071))


print(plane1.origin)

w = plane1.width/2
h = plane1.height/2
corner_x = plane1.origin[0]+plane1.x_axis[0]*w+plane1.y_axis[0]*h
corner_y = plane1.origin[1]+plane1.x_axis[1]*w+plane1.y_axis[1]*h
corner_z = plane1.origin[2]+plane1.x_axis[2]*w+plane1.y_axis[2]*h
x_vector = (plane1.x_axis[0]*(-50), plane1.x_axis[1]*(-50), plane1.x_axis[2]*(-50))
y_vector = (plane1.y_axis[0]*(-50), plane1.y_axis[1]*(-50), plane1.y_axis[2]*(-50))
z_vector = (plane1.z_axis[0]*80, plane1.z_axis[1]*80, plane1.z_axis[2]*80)


print(x_vector)

corner = (corner_x, corner_y, corner_z)

mimics.segment.crop_mask(mask_bone, bounding_box= mimics.BoundingBox3d(corner, x_vector, y_vector, z_vector))
#part_bone = mimics.segment.calculate_part(mask=mask_bone, quality='High')
#part_bone.name = "Bone"