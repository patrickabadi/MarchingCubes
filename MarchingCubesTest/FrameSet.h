#pragma once
#include "PclDefines.h"

class FrameSet
{
public:
  struct Item
  {
    PointCloudXYZRGBNormalPtr cloud;
    Eigen::Matrix4f transform;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
  std::vector<Item> items;

  bool Load ( std::string directory );
  size_t Size () { return items.size (); }

  inline Item &operator[] ( size_t n ) { return ( items[n] ); }

};

