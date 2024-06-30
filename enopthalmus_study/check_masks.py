import os  # for scandir() etc
import re  # for regexp matching

root = r'D:\Projects & Research\Enophthalmos Study'


def check_mask(side):
    # Find the Orbital Volume mask for this side. 
    # End the regex with '$' to skip any experimental or trial masks (e.g. with/without sinus)
    orbit_vol = mimics.data.masks.find(f'(?i){side}_Orbital Volume$', regex=True) # Use perl style ignore case flag

    if orbit_vol is None:
        # couldn't find the mask, so raise an insex error & bail
        raise (IndexError, ValueError)
        # Huston, we have a problem. Bail without returning results
        return

    # NEW - as we have changed the image set for the volume masks need to make sure that the image set they are linked to is active 
    mask_image = orbit_vol.image
    volume = orbit_vol.volume / 1000
    print(f'mask {orbit_vol.name} on {mask_image} is {volume} cc')



if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.
 
  # This version has one folder with the segmenting person int he file name
  
  root = r'D:\Projects & Research\Enophthalmos Study\re-do_DICOM'
  
  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]
  num_projects = len(projects)
  with mimics.disabled_gui_update():
    for i, p in enumerate(projects):
        mimics.file.open_project(filename=p, read_only_mode=True)
        print(f'project {i+1} of {num_projects} \t {p}')

        for side in ['left', 'right']:
            try:   
                check_mask(side)
            except (IndexError, ValueError):
                print(f'ERROR in {side} of {p} ({i}/{num_projects}')
                continue
        mimics.file.close_project()
