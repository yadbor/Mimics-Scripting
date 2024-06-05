# import required modules 
from collections import OrderedDict as od
import os
import sys
import subprocess
TEMPLATE = od([
            ("Patient" , ""),
            ("Study" , ""),
            ("Notes","")
            ])
TEMPLATE_3_MATIC = od([
            ("Processed" , "False")
            ])
################################
PATIENT_A = od([
            ("Patient" , "Mat patient"),
            ("Study" , "CT Heart scan"),
            ("Notes","")
            ])
################################
MIMICS_FILE_PATH = r"C:\MedData\DemoFiles"
MIMICS_FILE_NAME = "Heart.mcs"
PARTS_OF_INTEREST = ["LA", "LV", "Aorta"] 
################################

# One script is used for the metadata tutorial for both Mimics and 3-matic. 
# For that reason we have to check if we are in Mimics or in 3-matic 
in_mimics = False
try:
    import mimics
    if mimics.get_version():
        in_mimics = True
except:
    pass
# Mimics part
if in_mimics:
    parts = []
    # Open Mimics project
    mimics.file.open_project(os.path.join(MIMICS_FILE_PATH, MIMICS_FILE_NAME))
    # For this exercise we will remove all the metadata from the Mimics project
    for p in mimics.data.parts:
        for md in p.metadata:
            p.metadata.delete(md.name)
    # Group the required parts
    mdp = mimics.data.parts
    for p in PARTS_OF_INTEREST:
        parts.append(mdp[p])
    # Assign the template as metadata to all the parts of interest
    l = list(TEMPLATE.items())
    for p in parts:
        for i in range(len(TEMPLATE)):
            p.metadata.create(l[i][0],l[i][1])
    # Fill the metadata template
    patient_a = list(PATIENT_A.items())
    for p in parts:
        for i in range(len(PATIENT_A)):
            p.metadata[l[i][0]].value = patient_a[i][1]
    # Save Mimics project
    mimics.file.save_project()
    #Prepare to run 3-matic
    trimatic = mimics.file.get_path_to_3matic()
    command = trimatic 
    args = ("-run_script", __file__)
    process = subprocess.Popen((command,) + args, shell=False, stdout=subprocess.PIPE)

# 3-matic part
else:
    parts = []

    trimatic.import_project(os.path.join(MIMICS_FILE_PATH, MIMICS_FILE_NAME))
    # Group the required parts
    tp = trimatic.get_parts()
    for p in tp:
        if p.name in PARTS_OF_INTEREST:
            parts.append(p)
    # Assign the template as metadata elements to all the parts of interest
    l3m = list(TEMPLATE_3_MATIC.items())
    for p in parts:
        for i in range(len(TEMPLATE_3_MATIC)):
            mdata = p.get_metadata()
            mdata.create(l3m[i][0],l3m[i][1])
    # Smooth all the imported parts
    trimatic.smooth(entities = parts) 
    # Add the info that the parts are smoothed and processed
    l = list(TEMPLATE.items())
    for p in parts:
        mdata = p.get_metadata()
        notes = mdata.find(l[2][0],l[2][1])
        if notes:
            notes.value = "Part is smoothed with default values"
        processed = mdata.find(l3m[0][0],l3m[0][1])
        if processed:
            processed.value = "True"

                