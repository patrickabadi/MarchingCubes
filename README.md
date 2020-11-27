# MarchingCubes
Dynamic optimized mesh-generating Marching Cubes algorithm.  The goal of this project is to be able to give MarchingCubes successive point clouds along with a transformation (in world space) and it will continually update the mesh in the quickest amount of time.

## Requirements
- PointCloud Library [https://pointclouds.org](https://pointclouds.org)
- C++ 17

## Building the Project
- Use cmake to generate the project
- Steps I took to build on my Windows machine:

```
> mkdir vcpkg
> cd vcpkg
> git clone https://github.com/microsoft/vcpkg
> .\vcpkg\bootstrap-vcpkg.bat

> .\vcpkg install pcl:x64-windows

/* This command is used so other projects can more easily find and link to VCPKG libraries */
> .\vcpkg integrate install
> cd..

/* Now create your MarchingCubes folder */
> mkdir MarchingCubes
> cd MarchingCubes
> git clone https://github.com/patrickabadi/MarchingCubes
> mkdir build
> cd build
> cmake .. "-DCMAKE_TOOLCHAIN_FILE=D:\src\vcpkg\scripts\buildsystems\vcpkg.cmake"
```
## Running the Project
> Note: Debug timing will always be slower because of the debug iterator slowdowns

> Note: The debug build will build but fail to run unless you copy libpng16.dll and zlib1.dll from ./MarchingCubes/build/bin/Release to ./MarchingCubes/build/bin/Debug

- Run the app and on the visualization window rotate the scene until you see the model in front of you.  Zoom in/out with the mouse wheel
- Press spacebar to move to the next frame
- Press 's' to end early

<img src="./Screenshots/screenshot.png"/>

## Algorithm explained
- Incoming points are placed into a virtualized voxel grid.  This saves on space and allows it to be updated on subsequent incoming points, without having to recreate the voxel grid.
- Points added to each voxel grid now become vertices in a cube grid.  Every added point in the voxel grid can affect up to 8 cubes in the cube grids so they need to be checked
- Followed the following for the efficient generation of triangle faces [http://paulbourke.net/geometry/polygonise/](http://paulbourke.net/geometry/polygonise/)

## Contributers
- Patrick Abadi [https://github.com/patrickabadi](https://github.com/patrickabadi)
- Daniel Packard [https://github.com/daniel-packard](https://github.com/daniel-packard)