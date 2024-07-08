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

def get_basis_vector(img, get_origin = False):
  origin = img.get_voxel_center([0, 0, 0])
  dims = img.get_voxel_buffer().shape
  i = img.get_voxel_center([dims[0]-1, 0, 0])
  j = img.get_voxel_center([0, dims[1]-1, 0])
  k = img.get_voxel_center([0, 0, dims[2]-1])
  span = [np.asarray(v) - np.asarray(origin) for v in (i, j, k)]
  basis = [component / np.linalg.norm(component, ord=1) for component in span]
  if get_origin:
    return basis, origin
  else:
    return basis


for i, p in enumerate(projects):
  # process this file
  # for all image series
  for j, img in enumerate(mimics.data.images):
    # get the name of the series
    # get the pixel size and slice spacing 
    # get the basis vector