
import dummy_data

def box_corners(orig, x, y, z): 
  near = orig
  pts = [orig, x, y, z]
  far = [sum(idx) for idx in list(zip(* pts))]
  return (near, far)

s_points = dummy_data.left_eye.spline_points

(max_x, max_y, max_z) = [max(idx) for idx in list(zip(* s_points))]
(min_x, min_y, min_z) = [min(idx) for idx in list(zip(* s_points))]

print(f"Eye at {dummy_data.major.sphere.center}")
# Right side
x_lim = min_x - 10
x_delta = (max_x - min_x) + 20

box_orig = (x_lim, min_y - 5, min_z -15)
vector_x = [x_delta, 0, 0]
vector_y = [0, 80, 0]
vector_z = [0, 0, (max_z - min_z) + 30]

print("\nRight Ryan")
print(f"x_lim {x_lim} x_delta {x_delta}")
print(f"{box_orig} -> {vector_x} {vector_y} {vector_z}")

corners = box_corners(box_orig, vector_x, vector_y, vector_z)
print(corners)

# Left side
x_lim = max_x + 10
x_delta = (min_x - max_x) - 20

box_orig = (x_lim, min_y - 5, min_z -15)
vector_x = [x_delta, 0, 0]
vector_y = [0, 80, 0]
vector_z = [0, 0, (max_z - min_z) + 30]

print("\nLeft Ryan")
print(f"x_lim {x_lim} x_delta {x_delta}")
print(f"{box_orig} -> {vector_x} {vector_y} {vector_z}")

corners = box_corners(box_orig, vector_x, vector_y, vector_z)
print(corners)

print('------------------------------------------------')

s_points = dummy_data.right_eye.spline_points

(max_x, max_y, max_z) = [max(idx) for idx in list(zip(* s_points))]
(min_x, min_y, min_z) = [min(idx) for idx in list(zip(* s_points))]

print(f"Eye at {dummy_data.right_eye.sphere.center}")
# Right side
x_lim = min_x - 10
x_delta = (max_x - min_x) + 20

box_orig = (x_lim, min_y - 5, min_z -15)
vector_x = [x_delta, 0, 0]
vector_y = [0, 80, 0]
vector_z = [0, 0, (max_z - min_z) + 30]

print("\nRight Ryan")
print(f"x_lim {x_lim} x_delta {x_delta}")
print(f"{box_orig} -> {vector_x} {vector_y} {vector_z}")

corners = box_corners(box_orig, vector_x, vector_y, vector_z)
print(corners)

# Left side
x_lim = max_x + 10
x_delta = (min_x - max_x) - 20

box_orig = (x_lim, min_y - 5, min_z -15)
vector_x = [x_delta, 0, 0]
vector_y = [0, 80, 0]
vector_z = [0, 0, (max_z - min_z) + 30]

print("\nLeft Ryan")
print(f"x_lim {x_lim} x_delta {x_delta}")
print(f"{box_orig} -> {vector_x} {vector_y} {vector_z}")

corners = box_corners(box_orig, vector_x, vector_y, vector_z)
print(corners)
