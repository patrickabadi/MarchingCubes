#include "MarchingCubes.h"
#include <pcl/geometry/mesh_conversion.h>
#include <pcl/common/transforms.h>
#include <pcl/common/time.h>
#include <set>
#include "Common.h"

using namespace RTS;
using namespace pcl;
using namespace Eigen;
using namespace std;

/*
   * Tables, and functions, derived from Paul Bourke's Marching Cubes implementation:
   * http://paulbourke.net/geometry/polygonise/
   * Cube vertex indices:
   *   y_dir 4 ________ 5
   *         /|       /|
   *       /  |     /  |
   *   7 /_______ /    |
   *    |     |  |6    |
   *    |    0|__|_____|1 x_dir
   *    |    /   |    /
   *    |  /     |  /
   z_dir|/_______|/
   *   3          2
   */
const unsigned int edgeTable[256] = {
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};
const int triTable[256][16] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

struct VoxelInfo
{
public:
  Vector3i voxel_index;
  PointXYZRGBNormal point_transformed;
  VoxelInfo ()
  {
    voxel_index << -1, -1, -1;
  }

};

RTS::MarchingCubes::MarchingCubes ( float leaf_size )
  : _total_faces ( 0 )
  , _mesh_points ( new PointCloudXYZRGBNormal )
{
  _leaf_size[0] = leaf_size; _leaf_size[1] = leaf_size; _leaf_size[2] = leaf_size;

  // Use multiplications instead of divisions
  _inverse_leaf_size = Eigen::Array3f::Ones () / _leaf_size.array ();

  _size = (int)( 20.0f / leaf_size ) + 1;

  // setting up virtual voxel bounds
  _voxel_min_bounds << -10.0f, -10.0f, -10.0f;
  _voxel_max_bounds << 10.0f, 10.0f, 10.0f;

  _mesh_points->width = 0;
  _mesh_points->height = 1;
  _mesh_points->is_dense = false;
}

void RTS::MarchingCubes::Reset ()
{
  _voxel_grid.clear ();
  _cube_grid.clear ();
  _mesh_points->clear ();
  _mesh_points->width = 0;
  _total_faces = 0;
}

void RTS::MarchingCubes::Merge ( const PointCloudXYZRGBNormalPtr &input, const Eigen::Matrix4f &transformation )
{
  int original_cube_count = _cube_grid.size ();
  std::vector<Eigen::Vector3i> added_voxels;
  {
    pcl::ScopeTime t( __FUNCTION__ );
    added_voxels = Voxelize ( input, transformation );
    if ( added_voxels.size () == 0 )
      return;

    Cubelize ( added_voxels );
  }

  std::cerr << Format ( "%s - PointCloud Points: %d", __FUNCTION__, input->size() ) << endl;
  std::cerr << Format ( "%s - Added %d voxels, cubes: %d", __FUNCTION__, added_voxels.size (), _cube_grid.size() - original_cube_count ) << endl;
  std::cerr << Format ( "%s - Mesh vertices: %d, faces: %d", __FUNCTION__, _mesh_points->size(), _total_faces  ) << endl;
}



PointCloudXYZRGBNormalPtr RTS::MarchingCubes::BuildPointCloud ()
{
  pcl::ScopeTime t( __FUNCTION__ );
  auto output = PointCloudXYZRGBNormalPtr ( new PointCloudXYZRGBNormal );

  pcl::PointXYZRGBNormal pNan;
  pNan.x = pNan.y = pNan.z = std::numeric_limits<float>::quiet_NaN ();
  pNan.normal_x = pNan.normal_y = pNan.normal_z = std::numeric_limits<float>::quiet_NaN ();

  output->width = _voxel_grid.size ();
  output->height = 1;
  output->is_dense = false;
  output->points.resize ( _voxel_grid.size (), pNan );

  int index = 0;
  for ( const auto &p : _voxel_grid )
  {
    output->points[index++] = p.second.point;
  }

  return output;
}

pcl::PolygonMesh RTS::MarchingCubes::BuildMesh ()
{
  pcl::ScopeTime t( __FUNCTION__ );
  std::vector<pcl::Vertices> polys;
  pcl::PolygonMesh poly;


  pcl::Vertices p;
  p.vertices.emplace_back ( 0 );
  p.vertices.emplace_back ( 0 );
  p.vertices.emplace_back ( 0 );
  polys.resize ( _total_faces, p );

  int index = 0;
  for ( const auto &cube : _cube_grid )
  {
    for ( const auto &face : cube.second->face_indices )
    {
      auto &v = polys[index++];

      v.vertices[0] = face.vertices[0];
      v.vertices[1] = face.vertices[1];
      v.vertices[2] = face.vertices[2];
    }
  }

  poly.polygons = polys;

  pcl::PCLPointCloud2::Ptr cloud_blob ( new pcl::PCLPointCloud2 );
  pcl::toPCLPointCloud2 ( *_mesh_points, *cloud_blob );

  poly.cloud = *cloud_blob;

  return poly;
}

std::vector<Vector3i> RTS::MarchingCubes::Voxelize ( const PointCloudXYZRGBNormalPtr &input, const Eigen::Matrix4f &transformation )
{
  std::vector<Vector3i> added_voxels;

  if ( !input || input->size () == 0 )
    return added_voxels;

  added_voxels.reserve ( input->size () );

  std::vector< VoxelInfo > voxel_info ( input->size (), VoxelInfo () );

  Eigen::Affine3f affine_transformation;
  affine_transformation.matrix () = transformation.eval ();

#pragma omp parallel for
  for ( int i = 0; i < (int)input->size (); i++ )
  {
    const PointXYZRGBNormal &point = input->points[i];
    auto &voxel_item = voxel_info[i];

    if ( isnan ( point.x ) ||
      isnan ( point.y ) ||
      isnan ( point.z ) ||
      isnan ( point.normal_x ) ||
      isnan ( point.normal_y ) ||
      isnan ( point.normal_z ) )
      continue;

    PointXYZRGBNormal point_transformed = pcl::transformPointWithNormal ( point, affine_transformation );

    Vector3f p3 = point_transformed.getArray3fMap ();

    voxel_item.point_transformed = point_transformed;
    auto voxel_coord = ( p3 - _voxel_min_bounds );
    voxel_item.voxel_index <<
      int ( floor ( voxel_coord[0] * _inverse_leaf_size[0] ) ),
      int ( floor ( voxel_coord[1] * _inverse_leaf_size[1] ) ),
      int ( floor ( voxel_coord[2] * _inverse_leaf_size[2] ) );

  }

  for ( int i = 0; i < (int)input->size (); i++ )
  {
    const auto &voxel_item = voxel_info[i];
    if ( voxel_item.voxel_index[0] < 0 )
      continue;

    auto iter = _voxel_grid.find ( voxel_item.voxel_index );
    if ( iter != _voxel_grid.end () )
    {
      // we already have a point defined in that grid
      continue;

    }
    else
    {
      // add a new point
      _voxel_grid.insert ( make_pair ( voxel_item.voxel_index, VoxelItem ( voxel_item.point_transformed ) ) );
      added_voxels.emplace_back ( voxel_item.voxel_index );
    }
  }
  return added_voxels;
}

void RTS::MarchingCubes::Cubelize ( const std::vector<Vector3i> &added_voxels )
{
  std::set<std::tuple<int, int, int>> cube_set;

  for ( int i = 0; i < added_voxels.size (); ++i )
  {
    auto voxel_index = added_voxels[i];

    // there are 8 points in a cube, so eight potential cubes that could be affected by this point.
    // Using and std::set to avoid duplicates
    cube_set.insert ( std::make_tuple ( voxel_index[0], voxel_index[1], voxel_index[2] ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0] - 1, voxel_index[1], voxel_index[2] ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0], voxel_index[1] - 1, voxel_index[2] ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0] - 1, voxel_index[1] - 1, voxel_index[2] ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0], voxel_index[1], voxel_index[2] - 1 ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0] - 1, voxel_index[1], voxel_index[2] - 1 ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0], voxel_index[1] - 1, voxel_index[2] - 1 ) );
    cube_set.insert ( std::make_tuple ( voxel_index[0] - 1, voxel_index[1] - 1, voxel_index[2] - 1 ) );

  }

  for ( const auto &c : cube_set )
  {
    AddCube ( std::get<0> ( c ), std::get<1> ( c ), std::get<2> ( c ) );
  }

}

bool RTS::MarchingCubes::AddCube ( int cube_index_x, int cube_index_y, int cube_index_z )
{
  if ( cube_index_x < 0 || cube_index_x >= _size - 1 ||
    cube_index_y < 0 || cube_index_y >= _size - 1 ||
    cube_index_z < 0 || cube_index_z >= _size - 1 )
  {
    // out of bounds
    return false;
  }

  auto v0 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z ) );
  auto v1 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x + 1, cube_index_y, cube_index_z ) );
  auto v2 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x + 1, cube_index_y, cube_index_z + 1 ) );
  auto v3 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z + 1 ) );
  auto v4 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y + 1, cube_index_z ) );
  auto v5 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x + 1, cube_index_y + 1, cube_index_z ) );
  auto v6 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x + 1, cube_index_y + 1, cube_index_z + 1 ) );
  auto v7 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y + 1, cube_index_z + 1 ) );

  int cubeindex = 0;
  if ( v0 != _voxel_grid.end () ) cubeindex |= 1;
  if ( v1 != _voxel_grid.end () ) cubeindex |= 2;
  if ( v2 != _voxel_grid.end () ) cubeindex |= 4;
  if ( v3 != _voxel_grid.end () ) cubeindex |= 8;
  if ( v4 != _voxel_grid.end () ) cubeindex |= 16;
  if ( v5 != _voxel_grid.end () ) cubeindex |= 32;
  if ( v6 != _voxel_grid.end () ) cubeindex |= 64;
  if ( v7 != _voxel_grid.end () ) cubeindex |= 128;

  // Cube is entirely in/out of the surface
  if ( edgeTable[cubeindex] == 0 )
    return false;

  std::vector<PointXYZRGBNormal> edge_midpoints;
  auto defaultPoint = PointXYZRGBNormal ();
  defaultPoint.r = 255;
  defaultPoint.normal_x = 1.0f;
  defaultPoint.normal_y = 1.0f;
  defaultPoint.normal_z = 1.0f;
  edge_midpoints.resize ( 12, defaultPoint );


  std::shared_ptr<CubeItem> cube = GetCube ( cube_index_x, cube_index_y, cube_index_z, true );

  _total_faces -= cube->face_indices.size ();

  // remove any faces in this cube because we're going to create some new ones
  cube->face_indices.clear ();

  for ( int i = 0; triTable[cubeindex][i] != -1; i += 3 )
  {
    int cube_edge_index_0 = triTable[cubeindex][i];
    int vertex_index_0 = GetEdgePointIndex ( cube_index_x, cube_index_y, cube_index_z, cube_edge_index_0, cube );

    int cube_edge_index_1 = triTable[cubeindex][i + 1];
    int vertex_index_1 = GetEdgePointIndex ( cube_index_x, cube_index_y, cube_index_z, cube_edge_index_1, cube );

    int cube_edge_index_2 = triTable[cubeindex][i + 2];
    int vertex_index_2 = GetEdgePointIndex ( cube_index_x, cube_index_y, cube_index_z, cube_edge_index_2, cube );

    pcl::Vertices v;
    v.vertices.emplace_back ( vertex_index_0 );
    v.vertices.emplace_back ( vertex_index_1 );
    v.vertices.emplace_back ( vertex_index_2 );

    cube->face_indices.emplace_back ( v );
    _total_faces++;
  }

  Eigen::Vector3i cube_index ( cube_index_x, cube_index_y, cube_index_z );
  _cube_grid.insert_or_assign ( cube_index, cube );

  return true;
}

int RTS::MarchingCubes::GetEdgePointIndex ( int cube_index_x, int cube_index_y, int cube_index_z, int edge_index, std::shared_ptr<CubeItem> cube )
{
  // look at the edge indices at this link: http://paulbourke.net/geometry/polygonise/
  // each cube in our cube_grid only tracks edge index 0,8,3.  The other ones we retrieve from surrounding cubes

  switch ( edge_index )
  {
  case 1:
    return GetEdgePointIndex ( cube_index_x + 1, cube_index_y, cube_index_z, 3 );
  case 2:
    return GetEdgePointIndex ( cube_index_x, cube_index_y, cube_index_z + 1, 0 );
  case 4:
    return GetEdgePointIndex ( cube_index_x, cube_index_y + 1, cube_index_z, 0 );
  case 5:
    return GetEdgePointIndex ( cube_index_x + 1, cube_index_y + 1, cube_index_z, 3 );
  case 6:
    return GetEdgePointIndex ( cube_index_x, cube_index_y + 1, cube_index_z + 1, 0 );
  case 7:
    return GetEdgePointIndex ( cube_index_x, cube_index_y + 1, cube_index_z, 3 );
  case 9:
    return GetEdgePointIndex ( cube_index_x + 1, cube_index_y, cube_index_z, 8 );
  case 10:
    return GetEdgePointIndex ( cube_index_x + 1, cube_index_y, cube_index_z + 1, 8 );
  case 11:
    return GetEdgePointIndex ( cube_index_x, cube_index_y, cube_index_z + 1, 8 );
  case 0:
  case 3:
  case 8:
    // will be handled below
    break;
  default:
    // we got problems
    break;
  }

  if ( !cube )
  {
    cube = GetCube ( cube_index_x, cube_index_y, cube_index_z, true );
  }

  std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator v0;
  std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator v1;

  int cube_bit = 0;
  int mesh_index = 0;
  switch ( edge_index )
  {
  case 0:
    cube_bit = cube->edge_bit_0;
    mesh_index = cube->mesh_vertex_index_edge_0;
    v0 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z ) );
    v1 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x + 1, cube_index_y, cube_index_z ) );
    break;
  case 8:
    cube_bit = cube->edge_bit_0;
    mesh_index = cube->mesh_vertex_index_edge_0;
    v0 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z ) );
    v1 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y + 1, cube_index_z ) );
    break;
  case 3:
    cube_bit = cube->edge_bit_3;
    mesh_index = cube->mesh_vertex_index_edge_3;
    v0 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z ) );
    v1 = _voxel_grid.find ( Eigen::Vector3i ( cube_index_x, cube_index_y, cube_index_z + 1 ) );
    break;
  }

  int check_bit = 0;
  if ( v0 != _voxel_grid.end () )
    check_bit |= 2;
  if ( v1 != _voxel_grid.end () )
    check_bit |= 1;

  if ( check_bit != cube_bit || mesh_index < 0 )
  {
    PointXYZRGBNormal point;

    int vertex_index_0 = 0;
    int vertex_index_1 = 0;

    switch ( edge_index )
    {
    case 0:
      vertex_index_0 = 0;
      vertex_index_1 = 1;
      cube->edge_bit_0 = check_bit;
      break;
    case 8:
      vertex_index_0 = 0;
      vertex_index_1 = 4;
      cube->edge_bit_8 = check_bit;
      break;
    case 3:
      vertex_index_0 = 0;
      vertex_index_1 = 3;
      cube->edge_bit_3 = check_bit;
    }

    InterpolatePoint ( v0, v1, cube_index_x, cube_index_y, cube_index_z, vertex_index_0, vertex_index_1, point );

    if ( mesh_index < 0 )
    {
      _mesh_points->points.emplace_back ( point );
      _mesh_points->width++;
      mesh_index = _mesh_points->size () - 1;

      // save the mesh index back into the cube
      switch ( edge_index )
      {
      case 0:
        cube->mesh_vertex_index_edge_0 = mesh_index;
        break;
      case 8:
        cube->mesh_vertex_index_edge_8 = mesh_index;
        break;
      case 3:
        cube->mesh_vertex_index_edge_3 = mesh_index;
      }
    }
    else
    {
      _mesh_points->points[mesh_index] = point;
    }

  }

  return mesh_index;

}

std::shared_ptr<RTS::MarchingCubes::CubeItem> RTS::MarchingCubes::GetCube ( int cube_index_x, int cube_index_y, int cube_index_z, bool create_if_none )
{
  std::shared_ptr<CubeItem> cube = nullptr;
  Eigen::Vector3i cube_index ( cube_index_x, cube_index_y, cube_index_z );
  auto cube_check = _cube_grid.find ( cube_index );
  if ( cube_check == _cube_grid.end () )
  {
    if ( create_if_none )
    {
      cube = std::make_shared<CubeItem> ( CubeItem () );

      _cube_grid.insert ( make_pair ( cube_index, cube ) );

      return cube;
    }
    else
    {
      return nullptr;
    }
  }

  return cube_check->second;
}

void RTS::MarchingCubes::InterpolatePoint (
  std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator &v0,
  std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator &v1,
  int cube_x, int cube_y, int cube_z,
  int vertex_index_0, int vertex_index_1,
  PointXYZRGBNormal &output )
{
  if ( v0 == _voxel_grid.end () && v1 == _voxel_grid.end () )
  {
    // nothing to do
    return;
  }

  if ( v0 != _voxel_grid.end () && v1 != _voxel_grid.end () )
  {
    auto p0 = v0->second.point;
    auto p1 = v1->second.point;

    output.getArray3fMap () = ( p1.getArray3fMap () + p0.getArray3fMap () ) / 2.0f;
    output.getNormalVector3fMap () = ( p1.getNormalVector3fMap () + p0.getNormalVector3fMap () ) / 2.0f;
    output.rgba = p1.rgba;

  }
  else if ( v0 != _voxel_grid.end () )
  {
    auto p0 = v0->second.point;
    auto p1_map = GetCubePoint ( cube_x, cube_y, cube_z, vertex_index_1 );

    PointXYZ p1;
    p1.x = p1_map[0];
    p1.y = p1_map[1];
    p1.z = p1_map[2];

    output.getArray3fMap () = ( p1.getArray3fMap () + p0.getArray3fMap () * 3.f ) / 4.0f;
    output.getNormalVector3fMap () = p0.getNormalVector3fMap ();
    output.rgba = p0.rgba;
  }
  else
  {
    auto p0_map = GetCubePoint ( cube_x, cube_y, cube_z, vertex_index_0 );
    auto p1 = v1->second.point;

    PointXYZ p0;
    p0.x = p0_map[0];
    p0.y = p0_map[1];
    p0.z = p0_map[2];

    output.getArray3fMap () = ( p1.getArray3fMap () * 3.f + p0.getArray3fMap () ) / 4.0f;
    output.getNormalVector3fMap () = p1.getNormalVector3fMap ();
    output.rgba = p1.rgba;
  }
}

Eigen::Vector3f RTS::MarchingCubes::GetCubePoint ( int cube_x, int cube_y, int cube_z, int vertex_index )
{
  Eigen::Vector3f baseCorner ( cube_x, cube_y, cube_z );
  baseCorner[0] *= _leaf_size[0];
  baseCorner[1] *= _leaf_size[1];
  baseCorner[2] *= _leaf_size[2];
  baseCorner += _leaf_size / 2.0f;
  baseCorner += _voxel_min_bounds;

  switch ( vertex_index )
  {
  case 1:
    baseCorner[0] += _leaf_size[0];
    break;
  case 2:
    baseCorner[0] += _leaf_size[0];
    baseCorner[2] += _leaf_size[2];
    break;
  case 3:
    baseCorner[2] += _leaf_size[2];
    break;
  case 4:
    baseCorner[1] += _leaf_size[1];
    break;
  case 5:
    baseCorner[0] += _leaf_size[0];
    baseCorner[1] += _leaf_size[1];
    break;
  case 6:
    baseCorner[0] += _leaf_size[0];
    baseCorner[1] += _leaf_size[1];
    baseCorner[2] += _leaf_size[2];
    break;
  case 7:
    baseCorner[1] += _leaf_size[1];
    baseCorner[2] += _leaf_size[2];
    break;
  default:
    break;
  }

  return baseCorner;
}
