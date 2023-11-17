#Tutorial: Perform segmentation, region growing and calculate 3D to Mimi.mcs
# Open the project
input_dir = r'C:\MedData\DemoFiles\Mimi.mcs'
mimics.file.open_project(input_dir)
# Create an empty mask        
mask_a = mimics.segment.create_mask()
mask_a.name = "Bone"
# Perform thresholding with selected min and max values
mimics.segment.threshold(mask=mask_a,threshold_min=1250,threshold_max=2800) # thresholds are set in gray values
# Create a point that will be used fot the region growing operation
point_1 = mimics.analyze.indicate_point(title="Region growing point",message= "Please indicate a point on the part of interest")
point_2 = point_1.coordinates
point_2 = tuple(point_2)
# Region growing. The original mask is preserved
mask_b = mimics.segment.region_grow(point=point_2,input_mask=mask_a,target_mask=None,slice_type="Axial",keep_original_mask=True)
mask_b.name = "Segmented skull"
#Calculation of the 3D part
part_a = mimics.segment.calculate_part(mask=mimics.data.masks.find("Segmented skull"),quality="High")
# Export the STL
mimics.file.export_part(object_to_convert=part_a,file_name=r"C:\MedData\skull_of_Mimi.stl")
# Save the project and exit
mimics.file.save_project()
mimics.file.exit()
    
            