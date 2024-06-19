import mimics

mimics.file.open_project(filename=r'C:\MedData\DemoFiles\Skull_Cranial_Defect.mcs')

# Create the input mask
input_mask = mimics.segment.threshold(mask=mimics.segment.create_mask(), 
                                      threshold_min=mimics.segment.HU2GV(226), 
                                      threshold_max=mimics.segment.HU2GV(3071))
input_mask.name = 'input'

l_t = 1000
h_t = 3000
sd = 10

output_mask1 = mimics.segment.local_threshold(mask=input_mask, 
                                              threshold_min=l_t, 
                                              threshold_max=h_t, 
                                              search_distance=sd, 
                                              isolate=False, 
                                              bounding_box=None)
output_mask1.name = "auto"

output_mask2 = mimics.segment.activate_local_threshold(mask=input_mask, 
                                                       thresholds=None, 
                                                       search_distance=sd, 
                                                       clipping=None,
                                                       preview=None,
                                                       isolate=None,
                                                       is_modal=True)
output_mask2.name='manual'
