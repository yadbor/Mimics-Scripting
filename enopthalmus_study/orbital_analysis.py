
from const import * # CONSTANT definitions for this project (* is safe as only consts in file)
import utils # Utility functions to simplify code
import materials # Material definitions (used to make masks)

import mimics # API to access mimics (not needed inside mimics)
import mimics.analyze # API to access mimics (not needed inside mimics)

# Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
def make_crop_box(pt_up, pt_down, spacing_z):
  """Create a mimics.BoundingBox3D in the X,Y plane, given two plane intersection points.
  Input:  two points for the posterior edge of the box and a thickness.
          The posterior edge will be extended by MULT_XY in both directions and
          the box will extend SIZE_Y anteriorly.
  Output: A mimcs.BoundingBox3d to sue for cropping a mask."""

  # Align the crop box along the XY line between the two intesection points, 
  # expanded slightly to ensure it covers the whole orbit
  vector_x = [ MULT_XY * (pt_up.x - pt_down.x), MULT_XY * (pt_up.y - pt_down.y), 0]
  vector_y = [0, SIZE_Y, 0]    # -Y is anterior, so point towards front of face
  vector_z = [0, 0, spacing_z] # thickness is distance bwtween planes
  # Put the origin above the down intesection point, midway between the z planes
  origin = (pt_down.x, pt_down.y, pt_down.z + spacing_z / 2)
  return mimics.BoundingBox3d(origin, vector_x, vector_y, vector_z)

def make_orbit_mask(rim, globe):
  """Calculate Orbital contents volumes.
  Input:  a mimics.Spline delmiting the rim, and a mimics.Sphere for the globe.
  Output: a Mask from the 'surface' of the rim, extending forward, for cropping."""

  # Decompose the spline into straight lines between each pair of points.
  # Use these to intersect with planes, as mimics doesn't have a spline intersect plane function.
  # Use utils.looped_pairwise to join the last point back to the first point.
  spline_lines = [mimics.analyze.create_line(point1 = a, point2 = b) 
                  for a, b 
                  in utils.looped_pairwise(rim.points)]
  rim_geom = utils.spline_geometry(rim)
  # How tall is the orbital rim?
  delta_z = round(rim_geom.span[Z], 0) # Round the nearest whole number
  spacing = delta_z/NUM_PLANES
  
  # Make a list of plane origin points, matching the globe in the X,Y plane 
  # and spaced in Z from min_z to (min_z + delta_z) (approximately max_z)
  globe_x = globe.center[X]
  globe_y = globe.center[Y]
  min_z = rim_geom.min[Z]
 
  plane_origins = [(globe_x, globe_y, min_z + (z + 1) * spacing)
                    for z
                    in range(NUM_PLANES - 1)]

  # Create planes using these origins paralle to the X,Y plane (normal is Z+)
  norm_z = [0, 0, 1]
  z_planes = [mimics.analyze.create_plane_origin_and_normal(orig, norm_z) 
              for orig 
              in plane_origins]
  
  # For each plane, find intersections with the spline lines 
  # and create a bounding box based on the intersection points.
  boxes = [] # Start with an empty list for the boxes
  for plane in z_planes:
    pt_up, pt_down = None, None
    for line in spline_lines:
      # If a line crosses this plane it will have one endpoint above the plane and one below. 
      # Most lines will not cross this plane, in which case both these will be False
      from_above = (line.point1[Z] > plane.origin[Z] and line.point2[Z] < plane.origin[Z]) # TODO: should this be >=
      from_below = (line.point1[Z] < plane.origin[Z] and line.point2[Z] > plane.origin[Z]) # TODO: should this be =<
      if (from_above):
        pt_up = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)
      if (from_below):
        pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line, plane)
    
      # Each line will only cross any given plane once, so can stop
      # when there is an intersection in each direction.
      if (pt_up is not None) and (pt_down is not None):
        boxes.extend(make_crop_box(pt_up, pt_down)) # add the current crop box to the list
        break # We have found that line for this plane   
    else:
      # Fell through the loop without breaking
      print(f"did not find intersection for plane {plane}")

  # Adjust the first and last bounding box to ensure full overlap.
  boxes[0].third_vector[2] = -MULT_Z * boxes[0].third_vector[2]  # first goes down
  boxes[-1].third_vector[2] = MULT_Z * boxes[-1].third_vector[2] # and last goes up

  # Now use the bounding boxes to create a Region Of Interest
  mask_union = mimics.segment.create_mask() # Start with an empty mask
  # Loop through each box, create a mask, crop it with the box, and Union it together.
  for bbox in boxes:
    # Create a temporary mask, thresholding it to cover everything and
    # crop it with BoundingBox for this plane.
    m = mimics.segment.threshold(mask = mimics.segment.create_mask(), 
                                 threshold_min = materials.MIN_GV, 
                                 threshold_max = materials.MAX_GV,
                                 bounding_box = bbox)
    # Union with the existing mask
    mask_union = utils.unite(mask_union, m)
    # then delete it
    mimics.data.masks.delete(m)

  # The union of all the masks is now everything in front of the orbital rim
  return(mask_union)

def make_orbit_ROI(rim):
  # Get the rim extents 
  rim_geometry = utils.spline_geometry(rim)
  min_pt = rim_geometry.min_point
  max_pt = rim_geometry.max_point

  # The DICOM co-ordinate system is defined (for a BIPED) as patient LPS:
  #    X+ to the patient's Left
  #    Y+ to the patient's Posterior
  #    Z+ to the patient's Superior

  # Make a bounding box that extends beyond the rim by +/- 10 X and +/- 15 Z
  # with Y extending ant by -5 and posterior 80. Could go to max_pt[Y] + 80 ?
  expand = (10, 5, 15) # I think Y could be more?

  box_orig = (min_pt[X] - expand[X], min_pt[Y] - expand[Y], min_pt[Z] - expand[Z])
  vector_x = (max_pt[X] - min_pt[X] + (2 * expand[X]), 0, 0)
  vector_y = (0, 80, 0)
  vector_z = (0, 0, max_pt[Z] - min_pt[Z] + (2 * expand[Z]))
  
  orbit_ROI = mimics.BoundingBox3d(box_orig, vector_x, vector_y, vector_z)
  return orbit_ROI
