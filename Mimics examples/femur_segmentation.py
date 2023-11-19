#Tutorial: Perform segmentation and operations to the Hip.mcs. Calculate and export STL
# Open the project
mimics.file.open_project(r'C:\MedData\DemoFiles\Hip.mcs')
 # Create an empty mask        
mask_a = mimics.segment.create_mask()
mask_a.name = "Lower limb"
# Perform thresholding with selected min and max values
mimics.segment.threshold(mask=mask_a,threshold_min=1250,threshold_max=2650) # thresholds are set in gray values
#Fill holes in the segmentation mask
mimics.segment.fill_holes(mask_a)
# Create a point that will be used fot the region growing operation
point_1 = mimics.analyze.indicate_point(title="Region growing point",message= "Please indicate a point on the part of interest")
# Region growing. The original mask is preserved
mask_b = mimics.segment.region_grow(point=point_1,input_mask=mask_a,target_mask=None,slice_type="Axial",keep_original_mask=True)
#mimics.data.points.delete(point_1)
mask_b.name = "Segmented right femur"
# Perform the boolean operation ""Minus"" to take the anatomy of interest".
mask_c = mimics.segment.boolean_operations(mask_a=mimics.data.masks.find("Lower limb"), mask_b=mimics.data.masks.find("Segmented right femur"), operation="Minus")
mask_c.name ="Pelvis and left femur"
#Calculation of the 3D parts
part_a = mimics.segment.calculate_part(mask=mimics.data.masks.find("Segmented right femur"),quality="High")      
part_b = mimics.segment.calculate_part(mask=mimics.data.masks.find("Pelvis and left femur"),quality="High")      
# Smooth 3D parts 
objects = mimics.data.parts
for part in objects:
    part.visible = False
    smoothed_part = mimics.tools.smooth(object_to_smooth=part,smooth_factor=0.6,keep_originals=True)
    smoothed_part.visible = True
    # Export the STL
    mimics.file.export_part(object_to_convert=smoothed_part,file_name="C:\MedData\\" + smoothed_part.name + ".stl")
# Save the project and exit
mimics.file.save_project()
mimics.file.exit()
    
