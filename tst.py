
#import numpy as np
import itertools
# spline points
foo = list([[10,11,12], [12,11,13], [5,6,9], [2,5,10], [11,4,3], [8,6, 9]])
# make line segments from pairs of spline points
[print(a,b) for a, b in itertools.pairwise(foo)]
            
# for each plane, find point pair that crosses plane 
# going up & going down and calculate intersection

# for each pair of intersections,(i.e. per plane)
# make a bounding box aligned with points on x,y plane (the plane)
# and as thick as the plane slice in z
# then subtract from mask to get front surface
# of orbital volume
