# Tutorial for access to nodes of an STL
# Import the required libraries
try:
    import numpy as np  # First need to install numpy package for Python. Type pip install numpy in your cmd
except ImportError as ie:
    print("================================================================")
    print("=== The 3rd party Python package 'numpy' is not installed! ===")
    print("=== To install it, use 'pip install numpy' in your cmd!    ===")
    print("================================================================")
    raise
# Open Heart.mcs project
input_dir = r'C:\MedData\DemoFiles\Heart.mcs'
mimics.file.open_project(input_dir)
# Get the LV part
p = mimics.data.parts.find("LV")
if p is not None:
# Get a copy of nodes and triangles
    nodes,triangles = p.get_triangles()
# Read them with numpy
    nodes = np.asarray(nodes)
    print(len(nodes))
    triangles = np.asarray(triangles)
# Find the point that is closest to the WCS
    mx = []
    for m in nodes:
        mx.append(np.linalg.norm(m))
    i_mx = mx.index(max(mx))
# Find the point that is furthest from the WCS
    mn = []
    for m in nodes:
        mn.append(np.linalg.norm(m))
    i_mn = mn.index(min(mn))
# Calculate the distace
    d = mimics.measure.create_distance_measurement(list(nodes[i_mx]),list(nodes[i_mn]))
else:
    print("The part LV could not be found.")

