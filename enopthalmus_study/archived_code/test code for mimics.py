# Test code for mimics


# import numpy as np
p = mimics.data.parts[0]
v,t = p.get_triangles()
v = np.array(v)
t = np.array(t)
for i in range(len(v)):
    v[i] = v[i]+100
mimics.segment.create_part(v,t)


cube_v = [[0,0,0], [1,0,0], [1,0,1], [0,0,1],
          [0,1,0], [1,1,0], [1,1,1], [0,1,1]]
cube_t = [[0,1,2], [2,3,0], # front
          [2,1,5], [5,6,3], # left
          [4,0,7], [0,3,7], # right
          [7,6,5], [5,4,7], # back
          [2,6,7], [7,3,2], # top
          [5,1,0], [0,4,5]  # bottom
          ]


# tested and works
sph_r = mimics.data.spheres[0].radius
sph_c = mimics.data.spheres[0].center
imp_sph = mimics.file.import_stl(filename=r"D:\Projects & Research\Enophthalmos Study\Mimics-Scripting\enopthalmus_study\unit_sphere_level_3.ml.stl")
u_v, u_t = imp_sph.get_triangles()
u_t = np.asarray(u_t)
u_v = np.asarray(u_v)
u_v = u_v * sph_r
u_v = u_v + sph_c
new_sph = mimics.segment.create_part(u_v, u_t)
new_sph.name = "created"
sph_mask = mimics.segment.calculate_mask_from_part(new_sph)
sph_mask.name = "created_from_imported_sphere"