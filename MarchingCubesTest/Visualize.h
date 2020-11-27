#pragma once

#include <pcl/visualization/pcl_visualizer.h>
#include "PclDefines.h"

class Visualize
{
public:
  Visualize ();
  ~Visualize ();

  void Initialize ();

  void AddPointCloud ( PointCloudXYZRGBNormalPtr &cloud, std::string name, int r = -1, int g = -1, int b = -1 );
  void AddMesh ( PclMesh &cloud, std::string name, int r = -1, int g = -1, int b = -1 );
  void AddMesh ( PclMeshPtr &cloud, std::string name, int r = -1, int g = -1, int b = -1 );
  void AddPointCloudNormals(PointCloudXYZRGBNormalPtr& cloud, std::string name, int level = 10, int r = -1, int g = -1, int b = -1);
  void AddSphere ( std::string name, pcl::PointXYZ &center, float radius, int r, int g, int b );
  void AddLine ( std::string name, pcl::PointXYZ &from, pcl::PointXYZ &to, int r, int g, int b );
  void AddText ( std::string text );
  bool Display ( );
  std::string Display ( std::vector<std::string> keys );
  void Clear ( );


private:
  pcl::visualization::PCLVisualizer* _viewer;
};

