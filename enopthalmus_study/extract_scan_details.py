# Extract CT scan details from all Mimics projects in given folder
# Enophthalmous study 2024

import logging
import os  # for scandir() etc
import re  # for regexp matching
import time  # for timeing processes

import numpy as np  # Only used for np.array() in bounding box construction

# CONSTANT definitions for this project (* is safe as only consts)
from const import *
import utils  # Utility functions to simplify code

# from orbital_analysis import make_orbit_mask, make_orbit_ROI
import orbital_analysis
# Write measurements to a .CSV file
from result_logger import Path, log_to_file

import mimics  # API to access mimics

def extract_all_att(info):
    return (dict([(att, getattr(info, att)) for att in dir(info) if not re.match("^__", att)]))

# Get the names of all projects analyzed by Sam
root = r'D:\Projects & Research\Enophthalmos Study'
projects = [f.path for f in os.scandir(root) if re.match(r'.*SS\s+\d+\.mcs', f.name)]
num_projects = len(projects)

# Create a results file in the same folder as the project files
info_file = Path(os.path.join(os.path.dirname(projects[0]), 'project_info.csv'))

# Do the sanme for the DICOM tags
tags_file = Path(os.path.join(os.path.dirname(projects[0]), 'project_tags.csv'))
tag_headers = ['tag', 'VR', 'description', 'value']

for i, p in enumerate(projects):
    mimics.file.open_project(filename=p, read_only_mode=True)
    # Extract project information as a dict
    project_info = extract_all_att(mimics.file.get_project_information())
    # On first call this will create the file and write the headers, then the results.
    # Subsequent calls wil only write the results.
    log_to_file(info_file, project_info.keys(), project_info.values())

    # DICOM tags come as a dict already, but the parts need decoding
    # and different studies may have very different tags
    dicom_tags = mimics.get_dicom_tags()
    for key, tag in dicom_tags.items():
        tag_values = ['({0:04X} {1:04X})'.format(*key), tag.vr, tag.description, tag.value]
        log_to_file(tags_file, tag_headers ,tag_values)

