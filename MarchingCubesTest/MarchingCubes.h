#pragma once
#include "PclDefines.h"
#include <pcl/geometry/triangle_mesh.h>

namespace RTS
{ 
  template<class ArgumentType, class ResultType>
  struct unary_function
  {
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  };

  // https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
  // Hash function for Eigen matrix and vector.
  // The code is from `hash_combine` function of the Boost library. See
  // http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
  template<typename T>
  struct matrix_hash : unary_function<T, size_t>
  {
    std::size_t operator()( T const &matrix ) const
    {
      // Note that it is oblivious to the storage order of Eigen matrix (column- or
      // row-major). It will give you the same hash value for two different matrices if they
      // are the transpose of each other in different storage order.
      size_t seed = 0;
      for ( size_t i = 0; i < matrix.size (); ++i )
      {
        auto elem = *( matrix.data () + i );
        seed ^= std::hash<typename T::Scalar> ()( elem ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
      }
      return seed;
    }
  };

  class MarchingCubes 
  {

    class VoxelItem
    {
    public:
      pcl::PointXYZRGBNormal point;
      VoxelItem (const pcl::PointXYZRGBNormal& p)
        :point(p)
      {}

    };

    class CubeItem
    {
    public:
      std::vector<int> vertex_indices;
      std::vector<pcl::Vertices> face_indices;
      int mesh_vertex_index_edge_0;
      int mesh_vertex_index_edge_8;
      int mesh_vertex_index_edge_3;

      int edge_bit_0;
      int edge_bit_8;
      int edge_bit_3;

      CubeItem ()
        : mesh_vertex_index_edge_0 (-1)
        , mesh_vertex_index_edge_8 (-1)
        , mesh_vertex_index_edge_3 (-1)
        , edge_bit_0 (0)
        , edge_bit_8 (0)
        , edge_bit_3 (0)
      {
        vertex_indices.resize ( 12, -1 );
      }
    };

  public:
    MarchingCubes (float leaf_size = 0.05f);

    void Reset ();

    void Merge ( const PointCloudXYZRGBNormalPtr &input, const Eigen::Matrix4f &transformation );

    PointCloudXYZRGBNormalPtr BuildPointCloud ();
    pcl::PolygonMesh BuildMesh ();

  protected:
    int _size;
    Eigen::Vector3f _leaf_size;
    Eigen::Vector3f _inverse_leaf_size;
    Eigen::Vector3f _voxel_min_bounds;
    Eigen::Vector3f _voxel_max_bounds;
    PointCloudXYZRGBNormalPtr _mesh_points;

    int _total_faces;

    std::unordered_map<Eigen::Vector3i, VoxelItem, matrix_hash<Eigen::Vector3i>> _voxel_grid;
    std::unordered_map<Eigen::Vector3i, std::shared_ptr<CubeItem>, matrix_hash<Eigen::Vector3i>> _cube_grid;

    std::vector<Eigen::Vector3i> Voxelize ( const PointCloudXYZRGBNormalPtr &input, const Eigen::Matrix4f &transformation );
    void Cubelize ( const std::vector<Eigen::Vector3i>& added_voxels );
    bool AddCube ( int cube_index_x, int cube_index_y, int cube_index_z );
    int GetEdgePointIndex ( int cube_index_x, int cube_index_y, int cube_index_z, int edge_index, std::shared_ptr<CubeItem> cube = nullptr );
    void InterpolatePoint (
      std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator &v0,
      std::unordered_map<Eigen::Vector3i, VoxelItem>::iterator &v1,
      int cube_x, int cube_y, int cube_z,
      int vertex_index_0, int vertex_index_1,
      pcl::PointXYZRGBNormal &output );

    Eigen::Vector3f GetCubePoint ( int cube_x, int cube_y, int cube_z, int vertex_index );
    std::shared_ptr<CubeItem> GetCube ( int cube_index_x, int cube_index_y, int cube_index_z, bool create_if_none = false );

  };
}

