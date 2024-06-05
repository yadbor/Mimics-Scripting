# Import the numpy library that is useful for some operations
try:
    import numpy as np  # First need to install numpy package for Python. Type pip install numpy in your cmd
except ImportError as ie:
    print("================================================================")
    print("=== The 3rd party Python package 'numpy' is not installed! ===")
    print("=== To install it, use 'pip install numpy' in your cmd!    ===")
    print("================================================================")
    raise
# Define a shortcut
md = mimics.data 
# Constants declaration
MASK = "Threshold"
LANDMARKS = ("RA", "RA",
             "LA", "LA",
             "LV", "LV",
             "RV", "RV",
             "Aorta", "Aorta",
             "Pulmonary Artery", "Pulmonary Artery",
             )

SEED_RADIUS = dict(RA=10.0,
                   LA=10.0,
                   LV=10.0,
                   RV=10.0,
                   Aorta=8.0,
                   Pulmonary=8.0,
                   )

SEED_COLOR = dict(RA=(0, 255, 255),
                  LA=(255, 0, 255),
                  LV=(255, 205, 205),
                  RV=(145, 112, 255),
                  Aorta=(255, 0, 0),
                  Pulmonary=(0, 0, 255),
                  )

MASKS = ("RA", "LA", "LV", "RV", "Aorta", "Pulmonary Artery")

def activate_thresholding():
    m = mimics.segment.activate_thresholding()
    m.name  = MASK
    return

def indicate_landmark(pid: int):
    pdef = LANDMARKS[pid]
    name = pdef
    try:
        coords = mimics.indicate_coordinate(message="Indicate {} ".format(pdef),
                                                    confirm=False, show_message_box=True)
    except InterruptedError:
        return False

    mimics.view.navigate_to(coords)
    pnt = mimics.analyze.create_sphere_center_radius(coords, SEED_RADIUS[pdef.split()[0]])
    pnt.name = pdef
    pnt.color = tuple(np.array(SEED_COLOR[pdef.split()[0]]) / 255)
    return 

def calc_ct_heart():
    
    thres = md.masks.find(MASK)
    seeds = []
    for p in md.spheres:
       if p.name in LANDMARKS:
           seeds.append(p)        
    mimics.segment.calculate_ct_heart_from_mask(thres, seed_points=seeds)
    for p in md.masks:
       if p.name in MASKS:
                 p.color = tuple(np.array(SEED_COLOR[p.name.split()[0]])/255)
    return 

def create_3d_parts():
      for p in md.masks:
            if p.name in LANDMARKS:
                  par = mimics.segment.calculate_part(p,"Medium")
                  par.name = p.name
      return

# Main function
# Open Heart.mcs project
input_dir=r'C:\MedData\DemoFiles\Heart.mcs'
mimics.file.open_project(input_dir)
# Function for the thresholding 
activate_thresholding()
for l in LANDMARKS:
    # Function for the landmarking
      indicate_landmark(LANDMARKS.index(l))
# Function for the ct heart 
calc_ct_heart()
#Function for the creation of the 3d parts
create_3d_parts()
