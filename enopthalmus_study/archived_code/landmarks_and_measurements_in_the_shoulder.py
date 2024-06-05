#Tutorial: Landmarking the shoulder and perform measurements
# Open the project
mimics.file.open_project(r'C:\MedData\DemoFiles\Shoulder.mcs')
# Create an empty mask        
mask_a = mimics.segment.create_mask()
# Perform thresholding with selected min and max values
mimics.segment.threshold(mask=mask_a,threshold_min=1250,threshold_max=2800) # thresholds are set in gray values
mask_a.name = "Shoulder"
# Create a point that will be used fot the region growing operation
point_1 = mimics.analyze.indicate_point(title="Region growing point",message= "Please indicate a point on the part of interest")
# Region growing. The original mask is preserved
mask_b = mimics.segment.region_grow(point=point_1,input_mask=mask_a,target_mask=None,slice_type="Axial",keep_original_mask=True)
mimics.data.points.delete(point_1)
mask_b.name = "Segmented shoulder"
#Calculation of the 3D part
part = mimics.segment.calculate_part(mask=mimics.data.masks.find("Segmented shoulder"),quality="High")
# Set the anatomical landmarks of the shoulder        
anatomical_landmarks = ["Acromion","Coracoid process","Humerus"]
for point in anatomical_landmarks:
    p = mimics.analyze.indicate_point(title=point,message= "Please indicate a point on the {}".format(point))
    p.name = point
# Create distance measurement between coracoid & acromion and humerus
m = mimics.measure.create_distance_measurement(point1=mimics.data.points.find("Acromion").coordinates,point2=mimics.data.points.find("Humerus").coordinates)
m.name = "Acromion-Humerus"
m = mimics.measure.create_distance_measurement(point1=mimics.data.points.find("Coracoid process").coordinates,point2=mimics.data.points.find("Humerus").coordinates)
m.name = "Coracoid process-Humerus"
# Create Angle measurement between  the three landmarks in the shoulder area
mimics.measure.create_angle_measurement(point1=mimics.data.points.find("Acromion").coordinates,point2=mimics.data.points.find("Humerus").coordinates,point3=mimics.data.points.find("Coracoid process").coordinates)
# Save the project and exit
mimics.file.save_project()
mimics.file.exit()


