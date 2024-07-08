# Check a chosen project file
# Enophthalmous study 2024

import os  # for scandir() etc
import re  # for regexp matching

import mimics  # API to access mimics

# Define axis labels
X, Y, Z = 0, 1, 2

# Get the names of all projects analyzed by Sam
root = r'D:\Projects & Research\Enophthalmos Study'
re_good = re.compile('\\.mcs$', flags=re.IGNORECASE)
projects = [f.path for f in os.scandir(root) if re.match(re_good, f.name)]

# Count the files we are processing
num_projects = len(projects)

for i, p in enumerate(projects):
  # process this file
  # for all image series
  for img in p.images:
    # get the name of the series
    # get the pixel size and slice spacing 
    # get the basis vector
    
    
    
  

