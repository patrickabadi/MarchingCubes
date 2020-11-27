#include "FrameSet.h"
#include "Common.h"
#include <filesystem>
#include <pcl/io/pcd_io.h>

const std::string lcFilename = "lp_frame_%06d.";

bool FrameSet::Load ( std::string directory )
{
  items.clear ();

  int index = 0;

  std::string file;

  while ( true )
  {
    FrameSet::Item item;

    file = directory + Format ( lcFilename, index ) + "pcd";
    if ( !std::filesystem::exists ( file ) )
      break;

    item.cloud = PointCloudXYZRGBNormalPtr ( new PointCloudXYZRGBNormal );

    if ( pcl::io::loadPCDFile ( file, *item.cloud ) == -1 )
      break;

    file = directory + Format ( lcFilename, index ) + "txt";
    if ( !std::filesystem::exists ( file ) )
      break;

    std::ifstream ff ( file.c_str () );

    std::string line;
    std::getline ( ff, line, '\n' );

    std::istringstream ss ( line );

    std::string fval;
    float fv;
    for ( int y = 0; y < 4; y++ )
    {
      for ( int x = 0; x < 4; x++ )
      {
        std::getline ( ss, fval, ',' );
        sscanf_s ( fval.c_str (), "%f", &fv );
        item.transform ( y, x ) = fv;
      }
    }

    items.emplace_back ( item );

    index++;
  }

  return true;
}
