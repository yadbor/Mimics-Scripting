# Check a chosen project file
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching

import mimics  # API to access mimics

# Define axis labels
X, Y, Z = 0, 1, 2

# Get the names of all projects analyzed by Sam
root = r'D:\Projects & Research\Enophthalmos Study'
# Just the chosen good patients - and only the processed files
good = ('Evas', 'Dabbs', 'Hayes', 'Madden', 'Umbulgarri', 'Rogers')
re_good = re.compile('('+'|'.join(good)+').*\\d{2}\\.processed\\.mcs$', flags=re.IGNORECASE)
projects = [f.path for f in os.scandir(root) if re.match(re_good, f.name)]

# Count the files we are processing
num_projects = len(projects)

def rename_vols():
    vols = mimics.data.masks.filter('Orbital Volume')
    for m in vols: m.visible = True
    for m in vols: m.selected = True
    mimics.view.enable_mask_3d_preview() # Show the vols
    
    # Work out which side is which, as the names are broken form the other code :-(
    origins = [mimics.measure.get_bounding_box(objects=[v]).origin[X] for v in vols]
    if len(origins) ==1:
        if origins[0] < 0:
            side_prefixes = ('right',)
        else:
            side_prefixes = ('left',)
    else:
        if origins[0] < origins[1]:
            side_prefixes = ('right', 'left')
        else:
            side_prefixes = ('left', 'right')

    for i, j in enumerate(vols):
        j.name = side_prefixes[i] + j.name


for i, p in enumerate(projects):
    # Process project
    mimics.file.open_project(filename=p, read_only_mode=True)
    project_path = mimics.file.get_project_information().project_path
    folder = os.path.dirname(project_path)
    filename = os.path.basename(project_path)
    
    # Make pictures of the extracted volumes
    # Get the masks to image, and turn everything else off
    for m in mimics.data.masks: m.visible = False
    for p in mimics.data.parts: p.visible = False
    
    vols = mimics.data.masks.filter('Orbital Volume')
    for m in vols: m.visible = True
    for m in vols: m.selected = True
    mimics.view.enable_mask_3d_preview() # Show the vols
    
    # # Work out which side is which, as the names are broken form the other code :-(
    # origins = [mimics.measure.get_bounding_box(objects=[v]).origin[X] for v in vols]

    # if len(origins) == 0:
    #     continue # Skip this project as no orbital volume mask was found

    # if len(origins) == 1:
    #     if origins[0] < 0:
    #         side_prefixes = ('right',)
    #     else:
    #         side_prefixes = ('left',)
    # else:
    #     if origins[0] < origins[1]:
    #         side_prefixes = ('right', 'left')
    #     else:
    #         side_prefixes = ('left', 'right')

    # for i, j in enumerate(vols):
    #     j.name = side_prefixes[i] + '_' + j.name
        
    if len(vols) == 0:
        mimics.file.close_project()
        continue # Skip as no orbital volumes to see

    # Get the combined bounding box to zoom to
    orbital_bb = mimics.measure.get_bounding_box(objects=vols)
    # Zoom each view to cover both masks
    for v in mimics.data.views: 
        view_cam = mimics.view.get_camera(v)
        view_settings = view_cam.get_settings()
        view_settings.zoom_to_bounding_box(orbital_bb, zoom_factor=0.8)
        view_cam.set_settings(view_settings)
    
    # Save a screenshot for this project
    # Use the 3D view
    settings = mimics.view.get_camera(view = mimics.data.views['3D']).get_settings()
    settings.zoom_to_bounding_box(orbital_bb, zoom_factor=1)
    pic_filename = filename.replace('processed.mcs', 'snapshot.jpg')
    mimics.file.export_view_by_type(filename=os.path.join(folder, pic_filename), view='3D', image_type = 'jpg', camera_settings=settings)

    mimics.file.close_project() # Close and move to next