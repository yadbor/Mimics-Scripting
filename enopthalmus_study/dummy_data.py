# Dummy data to be returned in indicate_spline() and indicate_sphere()

from mimics.analyze import Sphere # get the Sphere class to return below

# A dataclass is like a record or c struct, with named attributes
from dataclasses import dataclass
# Define a class to hold the thresholds for each material
@dataclass
class dummy_data:
  sphere: Sphere
  spline_points: []

left_eye = dummy_data(
  sphere = Sphere(23.1270, 37.7933, -30.2142, 12.239),
  spline_points = [
            [11.04, 20.62,  -24.66],
            [25.89, 20.46,  -24.88],
            [37.23, 26.44,  -27.01],
            [41.37, 38.46,  -30.73],
            [39.41, 44.17,  -44.81],
            [28.93, 39.68,  -50.08],
            [16.13, 40.32,  -52.34],
            [7.46,  36.78,  -47.44],
            [3.52,  33.11,  -45.33],
            [1.55,  29.56,  -39.30],
            [1.03,  29.09,  -33.55],
            [6.15,  22.50,  -28.77]
           ]
)

right_eye = dummy_data(
  sphere = Sphere(-29.9621, -36.1569, 76.0900, 11.7760),
  spline_points = [
            [-22.3, -46.9, 94.00],
            [-32.2, -46.5, 92.10],
            [-41.6, -42.9, 87.38],
            [-46.0, -38.5, 83.18],
            [-45.9, -35.3, 75.76],
            [-46.5, -35.9, 68.80],
            [-44.1, -40.4, 61.95],
            [-36.3, -44.2, 58.72],
            [-26.3, -44.8, 60.14],
            [-16.6, -46.7, 64.92],
            [-11.0, -48.2, 70.07],
            [-8.8,  -46.8, 76.10],
            [-10.9, -44.7, 82.40],
            [-13.8, -45.9, 88.07],
            [-17.9, -46.5, 92.42]
           ]
)

dummy = right_eye # change this to use a different data set
