import sys
import numpy as np

from collections import namedtuple


Triangle = namedtuple("Triangle", "a,b,c")
Point = namedtuple("Point", "x,y,z")


def normalize(p):
    s = sum(u*u for u in p) ** 0.5
    return Point(*(u/s for u in p))


def midpoint(u, v):
    return Point(*((a+b)/2 for a, b in zip(u, v)))

 
def subdivide_edge(tri, depth):
    if depth == 0:
        yield tri
        return
    #       p0
    #      /  \
    # m01 /....\ m02
    #    / \  / \
    #   /___\/___\
    # p1    m12   p2
    p0, p1, p2 = tri
    m01 = normalize(midpoint(p0, p1))
    m02 = normalize(midpoint(p0, p2))
    m12 = normalize(midpoint(p1, p2))
    triangles = [
        Triangle(p0,  m01, m02),
        Triangle(m01, p1,  m12),
        Triangle(m02, m12, p2),
        Triangle(m01, m02, m12),
    ]
    for t in triangles:
        yield from subdivide_edge(t, depth-1)


def subdivide(faces, depth, method):
    for tri in faces:
        yield from method(tri, depth)


def simple_sphere(radius = 1, centre = [0, 0, 0], depth = 0):

    method = subdivide_edge
    # octahedron
    p = 2**0.5 / 2
    faces = [
        # top half
        Triangle(Point(0, 1, 0), Point(-p, 0, p), Point( p, 0, p)),
        Triangle(Point(0, 1, 0), Point( p, 0, p), Point( p, 0,-p)),
        Triangle(Point(0, 1, 0), Point( p, 0,-p), Point(-p, 0,-p)),
        Triangle(Point(0, 1, 0), Point(-p, 0,-p), Point(-p, 0, p)),

        # bottom half
        Triangle(Point(0,-1, 0), Point( p, 0, p), Point(-p, 0, p)),
        Triangle(Point(0,-1, 0), Point( p, 0,-p), Point( p, 0, p)),
        Triangle(Point(0,-1, 0), Point(-p, 0,-p), Point( p, 0,-p)),
        Triangle(Point(0,-1, 0), Point(-p, 0, p), Point(-p, 0,-p)),
    ]

    X = []
    Y = []
    Z = []
    T = []

    for i, tri in enumerate(subdivide(faces, depth, method)):
        X.extend([p.x for p in tri])
        Y.extend([p.y for p in tri])
        Z.extend([p.z for p in tri])
        T.append([3*i, 3*i+1, 3*i+2])

    # Scale unit sphere by radis and translate to centre
    X = np.array(X) * radius + centre[0]
    Y = np.array(Y) * radius + centre[1]
    Z = np.array(Z) * radius + centre[2]

    return (T, X, Y, Z)

def mimics_sphere(radius = 1, centre = [0, 0, 0], depth = 0):

    method = subdivide_edge
    # octahedron
    p = 2**0.5 / 2
    faces = [
        # top half
        Triangle(Point(0, 1, 0), Point(-p, 0, p), Point( p, 0, p)),
        Triangle(Point(0, 1, 0), Point( p, 0, p), Point( p, 0,-p)),
        Triangle(Point(0, 1, 0), Point( p, 0,-p), Point(-p, 0,-p)),
        Triangle(Point(0, 1, 0), Point(-p, 0,-p), Point(-p, 0, p)),

        # bottom half
        Triangle(Point(0,-1, 0), Point( p, 0, p), Point(-p, 0, p)),
        Triangle(Point(0,-1, 0), Point( p, 0,-p), Point( p, 0, p)),
        Triangle(Point(0,-1, 0), Point(-p, 0,-p), Point( p, 0,-p)),
        Triangle(Point(0,-1, 0), Point(-p, 0, p), Point(-p, 0,-p)),
    ]

    X = []
    Y = []
    Z = []
    T = []

    for i, tri in enumerate(subdivide(faces, depth, method)):
        X.extend([p.x for p in tri])
        Y.extend([p.y for p in tri])
        Z.extend([p.z for p in tri])
        T.append([3*i, 3*i+1, 3*i+2])

    # Scale unit sphere by radis and translate to centre
    X = np.array(X) * radius + centre[0]
    Y = np.array(Y) * radius + centre[1]
    Z = np.array(Z) * radius + centre[2]

    V = [[x, y, z] for x, y, z in list(zip(X, Y, Z))]

    return (np.asarray(T), np.column_stack((X, Y, Z)))

if __name__ == '__main__':

    (T, X, Y, Z) = simple_sphere(depth = 4)
    
    # Combine the X, Y, & Z components into a single array of vertices V
    #mimics.segment.create_part probably wants V & T (see below)
    V = [[x, y, z] for x, y, z in list(zip(X, Y, Z))]

    print(T)
    print('\n')
    
    # import numpy as np
    # p = mimics.data.parts[0]
    # v,t = p.get_triangles()
    # v = np.array(v)
    # t = np.array(t)
    # for i in range(len(v)):
    #     v[i] = v[i]+100
    # mimics.segment.create_part(v,t)
