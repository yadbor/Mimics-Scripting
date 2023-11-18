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
    self.p1 = None
    self.p2 = None
  
  def create(self, p1, p2):
    self.p1 = p1
    self.p2 = p2

class Plane:
  def __init__(self):
    self.origin = None
    self.normal = None

  def create(origin, normal):
    self = Plane()
    self.origin = origin
    self.normal = normal

# Object creation functions

def create_line(point1, point2):
  self = Line()
  self.point1 = point1
  self.point2 = point2
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
