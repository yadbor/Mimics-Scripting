# The DICOM co-ordinate system is defined (for a BIPED) as patient LPS. That is:
#    X+ to the Left hand side of the patient
#    Y+ to the Posterior (towards the back)
#    Z+ to the Superior (towards the head)
#
# Mimics appears to use the DICOM LPS co-ordinate system.
#
# Which eye is being analysed? The head is usually scanned in the centre 
# of the CT co-ordinate system, so X=0 is roughly the midline, meaning 
# the Right eye is X < 0 and the Left eye is X > 0
# TODO - a possibly more robust midline is the mean image X co-ordinate
if sphere1.x < 0:
  side = 'right'
elif sphere1.x > 0:
  side = 'left'
else:
  side = 'ambiguous'



###############################
# Old test code to create dummy data
# Supesceded by creating dummy mimics.analyze.indicate_{spline|sphere} code

  if TESTING: ##########################################################
    # Fake mimics geometry objects
    class Sphere:
      def __init__(self, x, y, z, radius):
        self.radius = radius
        self.center = [x, y, z]
    
    class Spline:
      def __init__(self, points):
        self.points = points

    # Define the geometry that the user would normally enter
    orbit = Sphere(11.776, -29.9621, -36.1569, 76.0900)
    spline_pts = [[-22.31, -46.94, 94.00],
                  [-32.22, -46.54, 92.10],
                  [-41.65, -42.97, 87.38],
                  [-46.08, -38.57, 83.18],
                  [-45.93, -35.38, 75.76],
                  [-46.55, -35.92, 68.80],
                  [-44.12, -40.42, 61.95],
                  [-36.31, -44.23, 58.72],
                  [-26.31, -44.84, 60.14],
                  [-16.61, -46.71, 64.92],
                  [-11.09, -48.23, 70.07],
                  [-8.86, -46.85, 76.10],
                  [-10.93, -44.75, 82.40],
                  [-13.85, -45.99, 88.07],
                  [-17.94, -46.54, 92.42,]
                 ]
    spline1 = Spline(points = spline_pts)
  else: ################################################################



###############################

# Ryan version
  intersect_points = {} # here len(intersect_points) == 0, so using len()+1 below starts indices from 1
  for n in range(len(z_planes)):
    for m in range(len(spline_lines)):
      if spline_lines[m].point1[2] > z_planes[n].origin[2] and spline_lines[m].point2[2] < z_planes[n].origin[2]:
        intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n])
      if spline_lines[m].point1[2] < z_planes[n].origin[2] and spline_lines[m].point2[2] > z_planes[n].origin[2]:
        intersect_points[len(intersect_points)+1] = mimics.analyze.create_point_as_line_and_plane_intersection(spline_lines[m], z_planes[n])

# Rob version
  intersect_points_RED = [] # List to store intersection points.
  for n, plane in enumerate(z_planes):
    # Go around the spline lines and for each line segment check if the y values span the plane.
    # If yes, then calculate the intercept. Going around the whole line gets both the upwards
    # and downwards intersections, one following the other so each pair of points is for the same plane
    for line in spline_lines:
      pt1_z = line.point1[Z]
      pt2_z = line.point2[Z]
      plane_z = plane.origin[Z]
      from_above = (pt1_z > plane_z and pt2_z < plane_z)
      from_below = (pt1_z < plane_z and pt2_z > plane_z)
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z] )
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z] )
      if (from_above or from_below):
        intersect_points_RED.append(mimics_analyze_create_point_as_line_and_plane_intersection(line, plane))

# Compare versions. 
# Ryan's is a dict so need to unpack the values to a list first
  all([i.__dict__ == j.__dict__ for i, j in zip([v for v in intersect_points.values()], intersect_points_RED)])
## Identical (all True)

# Rob version - by planes
  intersect_points_RED2 = len(z_planes) * [None] # List to store intersection points.
  for n, plane in enumerate(z_planes):
    # Go around the spline lines and for each line segment check if the y values span the plane.
    # If yes, then calculate the intercept. Going around the whole line gets both the upwards
    # and downwards intersections, one following the other so each pair of points is for the same plane
    for line in spline_lines:
      # Does ths line cross the curent plane? If so, in which direction?
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z])
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z])
      # Add the corresponding point to the list
      if from_above:
        pt_down = intersect_points_RED.append(mimics.analyze.create_point_as_line_and_plane_intersection(line, plane))
      if from_below:
        pt_up = intersect_points_RED.append(mimics.analyze.create_point_as_line_and_plane_intersection(line, plane))


# Ryan version 
  crop_masks = {} 
  for n in range(len(intersect_points)):   
    if n % 2 == 1:
   # If I read this correctly, this makes a line from points (1,2), (3, 4) etc.
   # What happened to point zero? ANS: we lost it in creating the list, as used [n+1] and n starts from 0
      vector_x = [(intersect_points[n].x-intersect_points[n+1].x)*1.2,(intersect_points[n].y-intersect_points[n+1].y)*1.2,0]
      print(f"vector_x from ({intersect_points[n].x}, { intersect_points[n+1].x}) and ({intersect_points[n].y}, {intersect_points[n+1].y})")
      vector_y = [0,-20,0]    
      if n == 1:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z+spacing_z/2)         
        vector_z = [0,0,-1.5*spacing_z]        
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)) 
      elif n == len(intersect_points)-1:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)         
        vector_z = [0,0,1.5*spacing_z]        
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
      else:
        crop_masks[len(crop_masks)+1] = mimics.segment.create_mask()
        ##mimics.segment.threshold(crop_masks[len(crop_masks)], mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))      
        box_orig = (intersect_points[n+1].x, intersect_points[n+1].y, intersect_points[n+1].z-spacing_z/2)          
        vector_z = [0,0,spacing_z] 
        mimics.segment.crop_mask(crop_masks[len(crop_masks)], mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z))
  
  print("Ryan result")      
  print('\n'.join(map(str, [i.__dict__ for i in crop_masks.values()])))

# Rob version
  united_masks = mimics.segment.create_mask()
  crop_masks_RED = []
  point_pairs = list(batched(intersect_points_RED, 2))
  last_pair = len(point_pairs) - 1
  for i, (pt_up, pt_down) in enumerate(point_pairs):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0] # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # extend first and last masks in z beyond their planes to cover whole orbit
    if i == 0:
      vector_z = [0, 0, -1.5 * spacing_z] # extends down past start
    elif i == last_pair:
      vector_z = [0, 0, +1.5 * spacing_z] # extends up past end
    
    # Set the bounding box origin above the down intesection point
    # and midway between the z planes
    box_orig = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    crop_box = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    
    # create the mask, threshold and crop it, then append to the list
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, crop_box)
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    crop_masks_RED.append(m)
    
    
  print("Rob Result")
  print('\n'.join(map(str, [m.__dict__ for m in crop_masks_RED])))


  point_pairs = list(batched(intersect_points_RED, 2))
  def make_crop_box(pt_up, pt_down):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0]       # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # Put the origin above the down intesection point, midway between the z planes
    origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)
  
  boxes = [make_crop_box(pt_up, pt_down) for pt_up, pt_down in point_pairs]
  # Make the First and Last box a bit bigger to ensure coverage
  boxes[0].third_vector[2] = -1.5 * boxes[0].third_vector[2]
  boxes[-1].third_vector[2] = 1.5 * boxes[-1].third_vector[2]

  # Create the united mask
  united_masks = mimics.segment.create_mask()
  for bbox in boxes:
    # Create a mask, threshold to cover everything
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, bbox) # Crop it with BoundingBox for this plane
    # Then Union with the existing mask
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    
  
  last_pair = len(point_pairs) - 1
  for i, (pt_up, pt_down) in enumerate(point_pairs):
    # Calculate the extents of the cropping box.
    # Align the crop box along the XY line between the two intesection points, 
    # expanded slightly to ensure it covers the whole orbit
    vector_x = [ 1.2 * (pt_up.x - pt_down.x ), 1.2 * (pt_up.y - pt_down.y ), 0]
    vector_y = [0, -20, 0] # -Y is anterior, so point towards front of face
    vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
    # extend first and last masks in z beyond their planes to cover whole orbit
    if i == 0:
      vector_z = [0, 0, -1.5 * spacing_z] # extends down past start
    elif i == last_pair:
      vector_z = [0, 0, +1.5 * spacing_z] # extends up past end
    
    # Set the bounding box origin above the down intesection point
    # and midway between the z planes
    box_orig = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
    crop_box = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
    
    # create the mask, threshold and crop it, then append to the list
    m = mimics.segment.create_mask()
    mimics.segment.threshold(m, const.MIN_GV, const.MAX_GV)
    mimics.segment.crop_mask(m, crop_box)
    united_masks = mimics.segment.boolean_operations(united_masks, m, 'Unite')
    crop_masks_RED.append(m)
    


# Part from triangles
# in Mimics API help under mimcs.segment
import numpy as np
# Scale matrix for a 10 unit sphere
def scale_matrix(s):
    return np.array([[s, 0, 0],[0, s, 0], [0, 0, s]])

def scale_object(p, r):
    tx = scale_matrix(r/10)
    v,t = p.get_triangles()
    v = np.array(v)
    t = np.array(t) # Don't actually need this line
    vt = v.dot(tx.T)
    return mimics.segment.create_part(vt,t)



p = mimics.data.parts[0]
v,t = p.get_triangles()
v = np.array(v)
t = np.array(t)
for i in range(len(v)):
    v[i] = v[i]+100
mimics.segment.create_part(v,t)




def spline_intercept_plane(spline_lines, plane):
    for line in spline_lines:
		from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z]) # TODO: should this be >=
		from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z]) # TODO: should this be =<
		if (from_above):
		  pt_up = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)
		if (from_below):
		  pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)






def bbox_from_intersections(pt_a, pt_b, mult = MULT_XY, thickness = spacing, depth = ANT_EXT):
	k = (mult - 1)/2.0
	delta = [(a-b) for a, b in zip(pt_a, pt_b)]
	ofs = [ k * i for i in delta]

	p1 = [a + o for a, o in zip(pt_a, ofs)]
	p2 = [b - o for b, o in zip(pt_b, ofs)]

	origin = p1[X], p1[Y], p1[Z] - (thickness/2)
	first_vector=[p2[X]-p1[X], p2[Y]-p1[Y], 0]
	second_vector = [0, -depth, 0]
	third_vector = [0, 0, thickness]
	bbox_ab = mimics.BoundingBox3d(origin, first_vector, second_vector, third_vector)
	#mask_ab = mimics.segment.threshold(mask=mimics.segment.create_mask(select_new_mask=False), threshold_min=materials.MIN_GV, threshold_max=materials.MAX_GV, bounding_box=bbox_ab)
	
	return bbox_ab


def bbox_thicken(bbox, extra):
	if extra < 0:
		bbox.origin = [bbox.origin[X], bbox.origin[Y], bbox.origin[Z] + extra]
		bbox.third_vector = [bbox.third_vector[X], bbox.third_vector[Y], bbox.third_vector[Z] - extra]
	else:
		bbox.third_vector = [bbox.third_vector[X], bbox.third_vector[Y], bbox.third_vector[Z] + extra]
	return bbox




src_part = mimics.data.parts.find('m_smooth')
for res in (0.2, 0.5): 
  for dist in range(10, 12): 
    name = f"wrap_{dist}x{res}"
    p_new = mimics.tools.wrap(src_part, smallest_detail=res, gap_closing_distance=dist,
                                dilate_result=False, protect_thin_walls=True, keep_originals=True)
    mimics.data.parts[-1].name = name

    print(name)



#################### Test bone segmentation again 2025-06-03
# Use the mimics Spongiosal Bone low threshold to cover all thin and cancellous bone parts.
m_bone_narrow = mimics.segment.threshold(mimics.segment.create_mask(), threshold_min=mimics.segment.HU2GV(226), threshold_max=mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
m_bone_narrow.name = 'm_bone_narrow'
m_bone = mimics.segment.threshold(mimics.segment.create_mask(), threshold_min=mimics.segment.HU2GV(148), threshold_max=mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
m_bone.name = 'm_bone'
m_bone = mimics.segment.smooth_mask(m_bone)

# Get air by thresholding, erode to seperate into internal and external regions
m_air = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200), bounding_box=orbit_ROI)
m_air.name = 'm_air'
m_air_erode = mimics.segment.morphology_operations(m_air, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name='m_air_erode', limited_to_mask=None)

# Get a point on the antero-lateral corner of the air mask to start region growing in the external part.
m_air_external_1 = mimics.segment.region_grow(input_mask=m_air_erode, target_mask=None, point=lateral_pt,
                                                slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')

# Dilate only the external part back to original size, then subtract that from the total to get the internal part
m_air_external = mimics.segment.morphology_operations(m_air_external_1, operation='Dilate', number_of_pixels=10, connectivity=26, target_mask_name='m_air_external', limited_to_mask=m_air)
# Subtract external from total air to get internal
m_air_internal = mimics.segment.boolean_operations(m_air, m_air_external, 'Minus')
m_air_internal.name = 'm_air_internal'

# Open internal air to remove small regions caused by low density tissue overlapping with our definition of air.
m_air_int_open = mimics.segment.morphology_operations(m_air_internal, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_open', limited_to_mask=None)
# Dilate what remains to extend to cover thin bone regions, then unite with bone
m_air_int_dilate = mimics.segment.morphology_operations(m_air_int_open, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_dilate', limited_to_mask=None)
m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
m_unite_bone_air.name = 'm_bone+air'
# Plan was to use local threshold to add more bone around the edges of the selected mask, but sadly, 
# local_threshold doesn't seem to work when called from the API, although it does using the GUI.
# Work around that by just using the lower threshold in the first place. This may get extra non-bone soft tissue, but that's hard to fix.

# Fill holes in the bone mask, then smooth it. Use keep_largest as some scans have air pockets in front of globe that breaks wrap
m_fill = mimics.segment.smart_fill_global(mask=m_bone, hole_closing_distance=3) # Is this distance too big?
m_fill = mimics.segment.smooth_mask(m_fill)
m_fill = mimics.segment.keep_largest(m_fill)
m_fill.name = 'm_fill'
# Fill holes in the bone+air mask, then smooth it. Use keep_largest as some scans have air pockets in front of globe that breaks wrap later.                                         
m_fill_ba = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=3) # Is this distance too big?
m_fill_ba = mimics.segment.smooth_mask(m_fill_ba)
m_fill_ba = mimics.segment.keep_largest(m_fill_ba)
m_fill_ba.name = 'm_fill_ba'
                                           
# Next turn the mask into a part, smooth and wrap it to get a hopefully watertight orbital border all around.
p_fill = mimics.segment.calculate_part(mask=m_fill, quality='High')
p_fill.name = 'p_fill'
# Note keep_originals=True for debugging - not used later so could be deleted here
p_smooth = mimics.tools.smooth(p_fill, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
p_smooth.name = 'p_smooth'
p_wrap = mimics.tools.wrap(p_smooth, smallest_detail=0.2, gap_closing_distance=6, dilate_result=False, protect_thin_walls=True, keep_originals=True)
p_wrap.name = 'p_wrap'

p_fill_ba = mimics.segment.calculate_part(mask=m_fill_ba, quality='High')
p_fill_ba.name = 'p_fill_ba'
# Note keep_originals=True for debugging - not used later so could be deleted here
p_smooth_ba = mimics.tools.smooth(p_fill_ba, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
p_smooth_ba.name = 'p_smooth_ba'
p_wrap_ba = mimics.tools.wrap(p_smooth_ba, smallest_detail=0.2, gap_closing_distance=6, dilate_result=False, protect_thin_walls=True, keep_originals=True)
p_wrap_ba.name = 'p_wrap_ba'

# Make the smoothed, wrapped part back into a mask. It should now define the orbit boundary
m_wrap = mimics.segment.calculate_mask_from_part(part=p_wrap, target_mask=None)
m_wrap.name = 'm_wrap'
m_wrap_ba = mimics.segment.calculate_mask_from_part(part=p_wrap_ba, target_mask=None)
m_wrap_ba.name = 'm_wrap_ba'

# Create a block of everything (this is ALL values - should be soft tissue + air? or all less than bone?)
m_all_tissue = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
m_all_tissue.name = 'm_all_tissue'

# Then subtract the boundary mask to separate the orbial contents from everything else.
m_temp_orbit = mimics.segment.boolean_operations(m_all_tissue, m_boundary, 'Minus')
m_temp_orbit_ba = mimics.segment.boolean_operations(m_all_tissue, m_boundary_ba, 'Minus')
# Use region growing to (hopefully) select only the internal volume
# This may fail, so erode first to try to seperate regions joined by thin necks or small holes
# To make sure there is enough left in the orbit to pick a seed point first combine with the globe
m_temp_orbit_1 = mimics.segment.boolean_operations(m_temp_orbit, m_globe, 'Minus')
m_temp_eroded = mimics.segment.morphology_operations(m_temp_orbit_1, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_open', limited_to_mask=None)

# Turn the wrapped part into a mask
m_wrap_ba = mimics.segment.calculate_mask_from_part(part=p_wrap_ba, target_mask=None)
m_wrap_ba.name = 'm_wrap_ba'

# Get the volume contained within the bony orbit

# Add the wrapped bone, original bone (to restore any corners), anterior block and external air
b1 = mimics.segment.boolean_operations(m_wrap_ba, m_fill_ba, operation='Unite')
b2 = mimics.segment.boolean_operations(m_ant_big, m_air_external, operation='Unite')
b3 = mimics.segment.boolean_operations(b1, b2, operation='Unite')
# Make a temporary 'all tissues' mask and subtract the combined boundary masks
temp_tissue = mimics.segment.threshold(mimics.segment.create_mask(), threshold_min=mimics.segment.HU2GV(-1024), threshold_max=mimics.segment.HU2GV(3071), bounding_box=mimics.measure.get_bounding_box((b3)))
t2 = mimics.segment.boolean_operations(temp_tissue, b3, operation='Minus')
# Temporarily add in the globe, so that the globe centre point will alwaye be in the mask
m_globe = mimics.data.masks.find('globe_mask') ## FIX THIS TO CORRECT OBJECT
t3 = mimics.segment.boolean_operations(t2, m_globe, operation='Unite')
p_globe = mimics.data.spheres[1]  ## FIX THIS TO CORRECT OBJECT
# Region grow from the centre of the globe, then subtract the globe again
t4 = mimics.segment.region_grow(input_mask=t3, target_mask=None, point=p_globe.center, slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
t5 = mimics.segment.boolean_operations(t4, m_globe, operation='Minus')
# The resulting mask may be a bit smaller than the bony orbit, so expand it a bit, then remove the filled bony orbit
t6 = mimics.segment.morphology_operations(input_mask=t5, operation='Dilate', number_of_pixels=2, target_mask_name='t6', connectivity=8,limited_to_mask=None)
t7 = mimics.segment.boolean_operations(t6, m_fill_ba, operation='Minus')
# Finally subtract the globe, as the dilate will have expanded into that region as well
t8 = mimics.segment.boolean_operations(t7, m_globe, operation='Minus')


# Failures 2024-06-05
failures = ((11, 424),
            (12, 424),
            (17, 499),
)


def mimics_image_orientation():
  active_img = [i for i in mimics.data.images if i.active][0]
  p0 = active_img.get_voxel_center([0, 0, 0])
  d = active_img.get_voxel_buffer().shape
  x = active_img.get_voxel_center([d[0]-1, 0, 0])
  y = active_img.get_voxel_center([0, d[1]-1, 0])
  z = active_img.get_voxel_center([0, 0, d[2]-1])
  return((p0, (x, y, z)))

# Usage:
p0, (x,y,z) = mimics_image_orientation()
# Calc spans for each axis
[tuple(b-a for a,b in zip(p0,v)) for v in (x, y, z)]

# Calc unit vectors for each axis
import numpy as np
basis = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
n_basis = basis / np.linalg.norm(basis, ord=1)

# Create a hoizontal plane in an oblique scan
# [00:38:58] ››› n_basis = basis / np.linalg.norm(basis, ord=1)
# [00:39:05] ››› n_basis[2]
# [00:39:05] [ 0.0365846 -0.27650366 0.69707959]
# [00:39:55] ››› new_p = mimics.analyze.create_plane_origin_and_normal(globe.center, n_basis[2])
# [00:39:55] CAD object created
#     type: Plane
#     name: Plane 25
#     color: [196, 196, 196, 128]
#     point1: [-38.6146, -356.1226, 687.8764]
#     point2: [-40.0221, -356.2587, 687.8963]
#     point3: [-39.4214, -355.0581, 688.3410]
def v_hat(v):
  mag = sum([i**2 for i in v])**0.5
  return([i/mag for i in v])

def mimics_image_vectors():
  active_img = [i for i in mimics.data.images if i.active][0]
  p0 = active_img.get_voxel_center([0, 0, 0])
  d = active_img.get_voxel_buffer().shape
  x = active_img.get_voxel_center([d[0]-1, 0, 0])
  y = active_img.get_voxel_center([0, d[1]-1, 0])
  z = active_img.get_voxel_center([0, 0, d[2]-1])
  if 'numpy' in sys.modules:
    # use the faster neater version
    basis = [np.asarray(v) - np.asarray(p0) for v in (x, y, z)]
    (i, j, k) = [b_i / np.linalg.norm(b_i, ord=1) for b_i in basis]
  else:
    basis = [tuple(b-a for a,b in zip(p0,v)) for v in (x, y, z)]
    (i, j, k) = [v_hat(v) for v in basis]
  return(i,j,k)

def mimics_z():
  active_img = [i for i in mimics.data.images if i.active][0]
  p0 = active_img.get_voxel_center([0, 0, 0])
  d = active_img.get_voxel_buffer().shape
  z = active_img.get_voxel_center([0, 0, d[2]-1])
   


# itertools versions of finding which side each component is on
import itertools

# globes

temp = [c.center[X] for c in itertools.islice(mimics.data.spheres, 2)]
if len(temp) > 0: 
    rg = which_min(temp)
    lg = 1 - rg
temp = [utils.spline_center(c)[X] for c in itertools.islice(mimics.data.spheres, 2)]
if len(temp) > 0: 
    rr = which_min(temp)
    lr = 1 - rr

# Consise & clever versions - maybe don't do this...
        # rg = which_min([s.center[X] for s in [mimics.data.spheres[i] for i in (0,1)]]); lg = 1 - rg
        # rr = which_min([utils.spline_center(s)[X] for s in [mimics.data.splines[i] for i in (0,1)]]); lr = 1 - rr

       

