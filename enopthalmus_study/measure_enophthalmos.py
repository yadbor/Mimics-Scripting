import os
import re
import csv
import numpy as np

import mimics

import result_logger

from const import * # Safe to use * as only contains CONST variables
import utils # My utility functions


def find_spline_plane_intersections(spline, plane):
  pt_up, pt_down = None, None
  norm = np.array(plane.normal)
  origin = np.array(plane.origin)
  # The spline is a continous loop. Check if each pair of points intersects this plane
  for p1, p2 in utils.looped_pairwise(spline.geometry_points):
    # If this segment of spline crosses this plane then one endpoint will be above the plane and one below. 
    p1_delta = np.array(p1) - origin
    p2_delta = np.array(p2) - origin
    p1_side = np.dot(p1_delta, norm)  #>0 = same side as normal pointing
    p2_side = np.dot(p2_delta, norm)  #<0 = opposite side to normal 
    # Most lines will not cross this plane, in which case both these will be False, so this is fast.
    if (p1_side >= 0 and p2_side <= 0): # crosses from above
      line_int = mimics.analyze.create_line(p1, p2) # temp line to get intersection point
      pt_up = mimics.analyze.create_point_as_line_and_plane_intersection(line_int, plane)
      mimics.data.lines.delete(line_int) # remove the temp line
  
    if (p1_side <= 0 and p2_side >= 0): # crosses from below
      line_int = mimics.analyze.create_line(p1, p2) # temp line to get intersection point
      pt_down = mimics.analyze.create_point_as_line_and_plane_intersection(line_int, plane)
      mimics.data.lines.delete(line_int) # remove the temp line
      
    # Each segment will only cross a given plane once. Stop when have found one in each direction.
    if (pt_up is not None) and (pt_down is not None):
      break # We have found that line for this plane   
  else:
    # Fell through the loop without breaking
    print(f"WARNING: did not find intersection for plane {plane}")

  return (pt_up, pt_down)

def find_extrema_points(splines, plane):
  intersection_points = [pts for s in splines for pts in find_spline_plane_intersections(spline=s, plane=plane)]
  intersection_array = np.array(intersection_points)
  max_idx = np.argmax(intersection_array, axis = 0)[0] # test the X axis and only return the x index
  min_idx = np.argmin(intersection_array, axis = 0)[0] # test the X axis and only return the x index
  return (intersection_points[min_idx], intersection_points[max_idx])

# face_pts = find_extrema_points(mimics.data.splines.filter("rim", regex=True), measure_plane)

def dist_line_to_points(p, q, rs):
    x = p-q
    return np.linalg.norm(
        np.outer(np.dot(rs-q, x)/np.dot(x, x), x)+q-rs,
        axis=1)

def measure_enophthalmos(globes, rims):
  # Create a plane through the centre of the globes and perpendicular to the frontal plane
  # This will be used to get the intersections of the two orbital rims.
  # The extreme points of those intersections wil define face_line, 
  # which is the baseline for measuring the globe positions.

  # OR - create the midpoint of the line between globes, then move it anteriorly in the Y direction, 
  # then create the plane with those three points. This is the same as being perpendicular to the frontal plane.

  pt1 = globes[0].center
  pt2 = globes[1].center
  pt_mid = mimics.analyze.create_midpoint(point1=pt1, point2=pt2)

  # Get the mimics basis vectors to establish what "anterior" is, as scan could be skewed or tilted.
  basis = utils.basis_vectors()
  pt_ant = np.array(pt_mid) + np.array((0, -20, 0)) * basis[Y]

  measure_plane = mimics.analyze.create_plane_points(point1=pt1, point2=pt2, point3=pt_ant)

  face_pts = find_extrema_points(rims, measure_plane)
  for p in face_pts:
    p.visible = True
    p.color = [(255/255), (68/255), (255/255)]
  line_1, line_2 = [np.array(p) for p in face_pts]
 
  distances = dist_line_to_points(line_1, line_2, np.array([pt1, pt2])) + np.array([g.radius for g in globes])

  return distances


if __name__ == '__main__':
  # Execute when the module is not initialized from an import statement.

  root = r'D:\Projects & Research\Enophthalmos Study\to_measure'

  results_file = result_logger.Path(os.path.join(root, 'enophthalmos.csv'))
  # Start with a clean empty file so that the headers are written correctly
  if os.path.exists(results_file):
    os.remove(results_file)

  # Get a list of all the .mcs files in root
  projects = [f.path for f in os.scandir(root) if re.match(r'.*.mcs', f.name)]
  num_projects = len(projects) 
 
  for i, p in enumerate(projects):
      project_name = os.path.basename(p)
      print(f"process {i+1}/{num_projects}\t{project_name}")

      if re.match(".*_measured.mcs", p):
        print("skip already measured file")
        continue

      mimics.file.open_project(filename=p, read_only_mode=True)

      headers = "project,left,right".split(",")
      try:
        # Use this method to get the eye order correct
        globes = [mimics.data.spheres[f"{side}_globe"] for side in ["left", "right"]]
        rims   = [mimics.data.splines[f"{side}_rim"]   for side in ["left", "right"]]

        en_left, en_right = measure_enophthalmos(globes, rims)

        results = [project_name[:-4], en_left, en_right]
      except KeyboardInterrupt:
        exit
      except Exception as e:
        results = [project_name[:-4], None, None]
        print("Error measuring distances: ", e)
      
      mimics.file.close_project()
      result_logger.log_to_file(results_file, headers, results)
  
  print("Finished.")