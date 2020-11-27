#include "Common.h"
#include "Visualize.h"
#include <string>

using namespace std;

bool next_iteration = false;
bool skip_visualize = false;

std::string last_key_press;

void keyboardEventOccurred ( const pcl::visualization::KeyboardEvent &event, void *nothing )
{
  if ( event.keyDown () )
  {
    last_key_press = event.getKeySym ();

    if ( last_key_press == "space" )
      next_iteration = true;
    else if ( last_key_press == "s" )
      skip_visualize = true;
  }

}

Visualize::Visualize ()
  :_viewer ( nullptr )
{

}

Visualize::~Visualize ()
{
  DEL ( _viewer )
}

void Visualize::Initialize ()
{
  DEL ( _viewer );

  _viewer = new pcl::visualization::PCLVisualizer ( "Virtual Scanner" );
  _viewer->setSize ( 800, 800 );  // Visualiser window size
  _viewer->registerKeyboardCallback ( &keyboardEventOccurred, (void *)NULL );
  _viewer->addCoordinateSystem ();
}

void Visualize::AddPointCloud ( PointCloudXYZRGBNormalPtr &cloud, std::string name, int r, int g, int b )
{
  if ( !_viewer )
    return;

  if ( !cloud )
  {
    if ( _viewer->contains ( name ) )
      _viewer->removePointCloud ( name );
    return;
  }

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGBNormal> one_colored_cloud ( cloud, (int)r, (int)g, (int)b );
  pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal> colored_cloud ( cloud );

  if ( _viewer->contains ( name ) )
  {
    // NOTE: I used to call viewer.updatePointCloud but if you call it too many times you end up with an error
    _viewer->removePointCloud ( name );
  }

  if ( r < 0 )
  {
    _viewer->addPointCloud<pcl::PointXYZRGBNormal> ( cloud, colored_cloud, name );
  }
  else
  {
    _viewer->addPointCloud<pcl::PointXYZRGBNormal> ( cloud, one_colored_cloud, name );
  }
}

void Visualize::AddMesh ( PclMeshPtr &cloud, std::string name, int r, int g, int b )
{
  if ( !cloud )
  {
    if ( _viewer->contains ( name ) )
      _viewer->removePolygonMesh ( name );
    return;
  }

  AddMesh ( *cloud, name, r, g, b );
}

void Visualize::AddMesh ( PclMesh & cloud, std::string name, int r, int g, int b )
{
  if ( !_viewer )
    return;

  if ( _viewer->contains ( name ) )
  {
    // NOTE: I used to call viewer.updatePointCloud but if you call it too many times you end up with an error
    _viewer->removePolygonMesh ( name );
  }

  _viewer->addPolygonMesh ( cloud, name, 0 );
}

void Visualize::AddPointCloudNormals ( PointCloudXYZRGBNormalPtr &cloud, std::string name, int level, int r, int g, int b )
{
  if ( !_viewer )
    return;

  if ( !cloud )
  {
    if ( _viewer->contains ( name ) )
      _viewer->removePointCloud ( name );
    return;
  }

  if ( _viewer->contains ( name ) )
  {
    // NOTE: I used to call viewer.updatePointCloud but if you call it too many times you end up with an error
    _viewer->removePointCloud ( name );
  }

  _viewer->addPointCloudNormals<pcl::PointXYZRGBNormal> ( cloud, level, 0.05, name );

  if ( !( r < 0 ) )
  {
    _viewer->setPointCloudRenderingProperties ( pcl::visualization::PCL_VISUALIZER_COLOR, (double)r / 255.0, (double)g / 255.0, (double)b / 255.0, name );
  }

}

void Visualize::AddSphere ( std::string name, pcl::PointXYZ &center, float radius, int r, int g, int b )
{
  _viewer->addSphere ( center, radius, r, g, b, name );
}

void Visualize::AddLine ( std::string name, pcl::PointXYZ &from, pcl::PointXYZ &to, int r, int g, int b )
{
  if ( !_viewer )
    return;

  if ( _viewer->contains ( name ) )
  {
    _viewer->removeShape ( name );
  }

  _viewer->addLine ( from, to, (double)r / 255.0, (double)g / 255.0, (double)b / 255.0, name );
}

void Visualize::AddText ( std::string text )
{
  if ( !_viewer )
    return;

  if ( _viewer->contains ( "text" ) )
  {
    _viewer->updateText ( text, 20, 20, "text" );
  }
  else
  {
    _viewer->addText ( text, 20, 20, "text" );
  }
}

bool Visualize::Display ()
{
  if ( !_viewer )
    return false;

  while ( !_viewer->wasStopped () )
  {
    _viewer->spinOnce ();

    if ( next_iteration )
    {
      next_iteration = false;
      return true;
    }

    if ( skip_visualize )
    {
      skip_visualize = false;
      return false;
    }
  }

  return true;
}

string Visualize::Display ( vector<string> keys )
{
  if ( !_viewer )
    return "";

  last_key_press.clear ();

  while ( !_viewer->wasStopped () )
  {
    _viewer->spinOnce ();

    if ( last_key_press.empty () )
      continue;

    for ( const auto &key : keys )
    {
      if ( last_key_press == key )
        return key;
    }

  }

  return last_key_press;
}

void Visualize::Clear ()
{
  if ( !_viewer )
    return;

  _viewer->removeAllShapes ();
  _viewer->removeAllPointClouds ();
  _viewer->removeAllCoordinateSystems ();
  _viewer->removeCorrespondences ();
  next_iteration = false;
  skip_visualize = false;
}
