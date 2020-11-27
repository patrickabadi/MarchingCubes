#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>

extern template class pcl::PointCloud<pcl::PointXYZRGBNormal>;
extern template class pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr;
extern template class boost::shared_ptr< pcl::PolygonMesh >;

typedef pcl::PointCloud<pcl::PointXYZRGBNormal> PointCloudXYZRGBNormal;
typedef pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr PointCloudXYZRGBNormalPtr;

typedef pcl::PolygonMesh PclMesh;
typedef boost::shared_ptr< pcl::PolygonMesh > PclMeshPtr;