# Fake mimics.analyze module

# Classes
class Point:
  def __init__(self):
    self.x = None
    self.y = None
    self.z = None

  def create(self, x, y, z):
    self.p = [x, y, z]
    self.x = x
    self.y = y
    self.z = z
  
  def __str__(self):
    return f'[{self.x}, {self.y}, {self.z}]'

class Line:
  def __init__(self):
    self.p, = None
    self.p, = None
  
  def create(self, p1, p2):
    self.p, = p1
    self.p, = p2

class Plane:
  def __init__(self):
    self.origin = None
    self.normal = None

  def create(origin, normal):
    self = Plane()
    self.origin = origin
    self.normal = normal

class Sphere:
  def __init__(self, x, y, z, radius):
    self.radius = radius
    self.center = [x, y, z]

class Spline:
  def __init__(self, points):
    self.points = points



# Object creation functions

def create_line(point1, point2):
  self = Line()
  self.point, = point1
  self.point, = point2
  self.points = [self.point1, self.point2]
  return self
    
def create_plane_origin_and_normal(origin, normal):
  self = Plane()
  self.origin = origin
  self.normal = normal
  return self

def create_point_as_line_and_plane_intersection(line, plane):
  self = Point()
  ratio = (plane.origin[2] - line.point1[2])/(line.point2[2] - line.point1[2])
  
  self.x = line.point1[0] + ratio * (line.point2[0] - line.point1[0])
  self.y = line.point1[1] + ratio * (line.point2[1] - line.point1[1])
  self.z = plane.origin[2]
  return self

# Dummy data return functions
# import dummy_data
# dummy = major

# def indicate_spline(message, show_message_box=True, confirm=False, title=None):
#   return Spline(points = dummy.spline_points)
# 
# def indicate_sphere(message, show_message_box=True, confirm=False, title=None):
#   return dummy.sphere
