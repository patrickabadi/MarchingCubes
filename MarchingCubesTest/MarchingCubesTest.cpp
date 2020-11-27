#include <iostream>
#include <string>
#include "FrameSet.h"
#include "Visualize.h"
#include "MarchingCubes.h"

int main (int argc, char* argv[]) 
{
  pcl::console::setVerbosityLevel ( pcl::console::L_ALWAYS );

  if (argc < 2)
  {
    printf ("Usage :\n");
    printf ("\t\t%s c:/test_data/\n", argv[0]);
    PCL_ERROR ("Provide a directory for your test data.\n");
    return (-1);
  }

  std::string directory = argv[1];

  FrameSet frame_set;
  RTS::MarchingCubes marching_cubes(0.05f); // 5cm voxels. Set your desired resolution
  Visualize viewer;
  viewer.Initialize ();

  frame_set.Load ( directory );
  if ( frame_set.Size () == 0 )
  {
    PCL_ERROR ( "No Files found in test data directory.\n" );
    return ( -1 );
  }

  bool display = true;
  Eigen::Matrix4f global_transform = Eigen::Matrix4f::Identity ();

  for ( int i = 0; i < frame_set.Size (); ++i )
  {
    FrameSet::Item &frame = frame_set[i];

    global_transform = global_transform * frame.transform;

    marching_cubes.Merge ( frame.cloud, global_transform );

    auto polygonMesh = marching_cubes.BuildMesh ();

    viewer.AddMesh ( polygonMesh, "mesh" );
    display = viewer.Display ();
    if ( !display )
      break;
  }
 
  return (0);
}
