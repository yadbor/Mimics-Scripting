# Tutorial for views
# Open Heart.mcs project
input_dir=r'C:\MedData\DemoFiles\Heart.mcs'
mimics.file.open_project(input_dir)
#Show and select the masks
for m in mimics.data.masks:
    if not m.visible:
        m.visible = True
    m.selected = True
# Delete the parts
for p in mimics.data.parts:
    mimics.data.parts.delete(p)
# Activate 3D preview
mimics.view.enable_mask_3d_preview()
mimics.dialogs.question_box(message="Please inspect the heart segmentation",buttons='OK')
#Create 3D parts
for m in mimics.data.masks:
    p = mimics.segment.calculate_part(mask=m, quality='Optimal')
    p.name = m.name
    p.visible = True
# Disable 3D preview
mimics.view.disable_mask_3d_preview()
#Preparation for fluoroscopy
visualised_objects = []
contrast = 0.7
for p in mimics.data.parts:
	visualised_objects.append((p,contrast))
# Activate fluoroscopy
f = mimics.view.create_fluoroscopy_view_default()
# Activate simulation
sim_quality = "High"
f.simulate(objects_contrast=visualised_objects,quality=sim_quality)
