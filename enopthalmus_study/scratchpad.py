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

import mimics # for sytax checker and dummy routines

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
p = mimics.data.parts[0]
v,t = p.get_triangles()
v = np.array(v)
t = np.array(t)
for i in range(len(v)):
    v[i] = v[i]+100
mimics.segment.create_part(v,t)

import numpy as np
# Create a Scale matrix
def scale_matrix(s):
    return np.array([[s, 0, 0],[0, s, 0], [0, 0, s]])
# Scale an object by s (in our case, r/10 for a 10 unit sphere)
def scale_object(p, s):
    tx = scale_matrix(s)
    v,t = p.get_triangles()
    v = np.array(v)
    t = np.array(t) # Don't actually need this line 
    vt = v.dot(tx.T)
    return mimics.segment.create_part(vt,t)


# Text crop box calculation
# Make a tuple of spans (min_x, max_x)

spans = ((10, 20), (-20, -10))
for mn, mx in spans:
  box_orig_pos = mn - 10
  vector_x_pos = mx - mn + 20
  box_orig_neg = mx + 10
  vector_x_neg = mn - mx - 20
  print(f"box ({mn}, {mx}) neg from {box_orig_neg} by {vector_x_neg} to {box_orig_neg + vector_x_neg}, pos from {box_orig_pos} by {vector_x_pos} to {box_orig_pos + vector_x_pos}")
