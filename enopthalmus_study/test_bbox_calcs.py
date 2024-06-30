import numpy as np
def bbox_to_points(bbox):
    p1 = bbox.origin
    # Add the three vectors
    span = [(a + b + c) 
    for a, b, c 
    in zip(bbox.first_vector, bbox.second_vector, bbox.third_vector)]
    p2 = tuple((p + s) for p, s in zip(p1, span)) # Return a tuple so both points the same
    return (p1, p2)

#Find the unit basis vectors for the project.
# These may not be (1,0,0), (0,1,0), (0,0,1) if the project is tilted
def v_hat(v):
    '''Return a normalised unit vector.''' 
    mag = sum([i**2 for i in v])**0.5
    return([i/mag for i in v])

def mimics_image_vectors():
    '''return the three unit vectors that descrivbe this image volume.'''
    active_img = [i for i in mimics.data.images if i.active]#
    p0 = active_img.get_voxel_center([0, 0, 0])
    d = active_img.get_voxel_buffer().shape
    x = active_img.get_voxel_center(#-1, 0, 0])
    y = active_img.get_voxel_center([0, d#-1, 0])
    z = active_img.get_voxel_center([0, 0, d#)
    if 'numpy' in sys.modules:
    # use the faster neater version
    basis = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
    (i, j, k) = [b_i / np.linalg.norm(b_i, ord=1) for b_i in basis]
    else:
    basis = [tuple(b-a for a,b in zip(p0,v)) for v in (x, y, z)]
    (i, j, k) = [v_hat(v) for v in basis]
    
    return (i,j,k)

#Project edition information
#    Project created in: Materialise Mimics Medical 25.0
#    Project last modified in: Materialise Mimics Medical 25.0
#Open project
#    file name: D:\Projects & Research\Enophthalmos Study\Dabbs 2015.11.14 merged SS 01.processed.mcs

left_m_air = mimics.data.masks.find('left_m_air')
left_m_air.maximum_value, left_m_air.minimum_value
#(824, 0)
m_air_bb = mimics.measure.get_bounding_box(#)
m_air_bb
#<mimics.BoundingBox3d((-17.78324137441814, -415.53645670134574, -331.625), (96.73704539053142, -1.180451363325119e-05, 0.0), (1.0244548320770264e-07, 143.59886514861137, 0.0), (1.0244548320770264e-07, -1.180451363325119e-05, 56.03125))>
m_air_bb.origin
#(-17.78324137441814, -415.53645670134574, -331.625)
p1, p2 = bbox_to_points(m_air_bb)
p1, p2
#((-17.78324137441814, -415.53645670134574, -331.625), (78.95380422100425, -271.93761516176164, -275.59375))
v1, v2, v3 = m_air_bb.first_vector, m_air_bb.second_vector, m_air_bb.third_vector
dir(m_air_bb)
#['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'first_vector', 'origin', 'second_vector', 'third_vector']
np.array(v1) + np.array(v2) + np.array(v3)
#[ 96.7370456 143.59884154 56.03125 ]
v1, v2, v3
#((96.73704539053142, -1.180451363325119e-05, 0.0), (1.0244548320770264e-07, 143.59886514861137, 0.0), (1.0244548320770264e-07, -1.180451363325119e-05, 56.03125))
p1 = m_air_bb.origin
p2 = np.array(p1) + (np.array(v1) + np.array(v2) + np.array(v3))
p1, p2
#((-17.78324137441814, -415.53645670134574, -331.625), array([ 78.95380422, -271.93761516, -275.59375 ]))
p2 = np.array(p1) + np.array(v1) + np.array(v2) + np.array(v3)
p1, p2
#((-17.78324137441814, -415.53645670134574, -331.625), array([ 78.95380422, -271.93761516, -275.59375 ]))

p2 = np.array(p1) + np.array(v1) + np.array(v2) + np.array(v3)
# crop mask has p1=[-90.9875, -394.5027, -363.0781], p2=[95.4149, -272.9790, -275.0781]
m_air_bb.origin
#(-17.78324137441814, -415.53645670134574, -331.625)
#Crop Mask dialog is closed
# But crop mask was the whole volume...
i, j, k = mimics_image_vectors()
i,j,k
#(array([ 0.87810006, -0.12189994, 0. ]), array([0.12189993, 0.87810007, 0. ]), array([0., 0., 1.]))
delta = np.array(v1) + np.array(v2) + np.array(v3)
delta
#[ 96.7370456 143.59884154 56.03125 ]
np.array(p1) + np.array(delta)
#[ 78.95380422 -271.93761516 -275.59375 ]
p2 = np.array(p1) + np.array(delta)
v_align = [m * np.array(v) for m, v in zip(delta, #
v_align
#[array([ 84.94480508, -11.79224051, 0. ]), array([ 17.50468905, 126.09415249, 0. ]), array([ 0. , 0. , 56.03125])]
bb_align = mimics.BoundingBox3d(origin=p1, first_vector=v_align#, second_vector=v_align#, third_vector=v_align#)
m_align = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_align)

# match crop mask to mask and it reads p1 = [-16.5915, -404.8305, -331.4531], p2 = [84.4474, -300.7347, -275.0781]
p1, p2
#((-17.78324137441814, -415.53645670134574, -331.625), array([ 78.95380422, -271.93761516, -275.59375 ]))
v_align
#[array([ 84.94480508, -11.79224051, 0. ]), array([ 17.50468905, 126.09415249, 0. ]), array([ 0. , 0. , 56.03125])]
np.array(p1) + v_align
#[[ 6.71615637e+01 -4.27328697e+02 -3.31625000e+02]
#[-2.78552329e-01 -2.89442304e+02 -3.31625000e+02]
#[-1.77832414e+01 -4.15536457e+02 -2.75593750e+02]]
np.sum(v_align)
#272.7826561090635
np.sum(v_align, axis=0)
#[102.44949413 114.30191198 56.03125 ]

v_align
#[array([ 84.94480508, -11.79224051, 0. ]), array([ 17.50468905, 126.09415249, 0. ]), array([ 0. , 0. , 56.03125])]
np.array(p1) + np.sum(v_align, axis=0)
#[ 84.66625275 -301.23454472 -275.59375 ]
created_p2 = np.array(p1) + np.sum(v_align, axis=0)
created_bb = mimics.measure.get_bounding_box(mimics.data.masks.find('Turquoise'))
created_bb
#<mimics.BoundingBox3d((-17.78324137441814, -415.53645670134574, -331.625), (96.73704539053142, -1.180451363325119e-05, 0.0), (1.0244548320770264e-07, 143.59886514861137, 0.0), (1.0244548320770264e-07, -1.180451363325119e-05, 56.03125))>
#Crop mask
    Corner 1: [-16.5915, -404.8305, -331.4531]
    Corner 2: [84.4474, -300.7347, -275.0781]
#Crop Mask dialog is closed
pa = (-46.3041, -359.4182, -295.1016)
pb = (-7.6911, -370.3003, -295.1350)
line_ab = mimics.analyze.create_line(point1=pa, point2=pb, name='across_orbit')

v_ab = np.array(pb) - np.array(pa)
bb_angled = mimics.BoundingBox3d(origin=pa, first_vector=v_ab, second_vector=v_align#, third_vector=v_align#)
m_angled = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_angled)
v_ba = np.array(pa) - np.array(pb)
bb_angled2 = mimics.BoundingBox3d(origin=pa, first_vector=v_ba, second_vector=v_align#, third_vector=v_align#)
m_angled2 = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_angled2)
bb_angled3 = mimics.BoundingBox3d(origin=pb, first_vector=v_ba, second_vector=v_align#, third_vector=v_align#)
m_angled3 = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_angled3)

bb_angled4 = mimics.BoundingBox3d(origin=pa, first_vector=v_align#, second_vector=-v_align#, third_vector=v_align#)
m_angled4 = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_angled4)

bb_angled4 = mimics.BoundingBox3d(origin=pa, first_vector=v_ab, second_vector=-v_align#, third_vector=v_align#)
m_angled4 = mimics.segment.threshold(mask=mimics.segment.create_mask(), threshold_min=0, threshold_max=4095, bounding_box=bb_angled4)
