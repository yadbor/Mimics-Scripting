# Dummy mimics module
# for top-level objects and methods ie.e mimics.method()

from collections import namedtuple

BoundingBox3d = namedtuple('bbox', 'origin first_vector second_vector third_vector')

#def BoundingBox3d(orig, vector_x, vector_y, vector_z):
#  return [orig, vector_x, vector_y, vector_z]

class Object:
  def __init__(self):
    self.name = ''
    return
  
