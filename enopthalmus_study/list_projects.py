# Measure the intersection volumes for various tissues, 
# given folders of projects with already segmented orbits.
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching

import numpy as np  # Was only used for np.array() in bounding box construction; now for several more

# CONSTANT definitions for this project (* is safe as only consts)
from const import *
import utils  # Utility functions to simplify code

import mimics  # API to access mimics


def list_a_project(p):
    project_info = mimics.file.get_project_information()
    path = project_info.path
    filename = os.path.basename(path)


def list_an_image_set(img):
    info = img.get_image_information()
    description = (
        f"algorithm:\t{info.algorithm}\n"
        f"obliqueness:\t{info.obliqueness}\n"
        f"orientation:\t{info.orientation}\n"
        f"gantry_tilt:\t{info.gantry_tilt}\n"
        f"slice_thickness:\t{info.slice_thickness}\n"
        f"slice_increment:\t{info.slice_increment}\n"
        f"pixel_size:\t{info.pixel_size}\n"
    )
    return description



