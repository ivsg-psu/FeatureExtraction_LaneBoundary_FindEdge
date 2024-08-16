
# FeatureExtraction_LaneBoundary_FindEdge

<!--
The following template is based on:
Best-README-Template
Search for this, and you will find!
>
<!-- PROJECT LOGO -->
<br />
<p align="center">
  <!-- <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h2 align="center"> FeatureExtraction_LaneBoundary_FindEdge
  </h2>

  <pre align="center">
    <img src=".\Images\RaceTrack.jpg" alt="main laps picture" width="960" height="540">
    <!--figcaption>Fig.1 - The typical progression of map generation.</figcaption -->
    <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

  <p align="center">
    This library finds the road edge from LIDAR data and XYZ trajectory data of a mapping vehicle. The road edge is the location where the pavement stops being a drivable surface, usually the edge of the pavement and the vegetation next to the road.
    <br />
    <!-- a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents">View Demo</a>
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Report Bug</a>
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Request Feature</a -->
  </p>
</p>

***

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About the Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="structure">Repo Structure</a>
      <ul>
        <li><a href="#directories">Top-Level Directories</li>
        <li><a href="#dependencies">Dependencies</li>
      </ul>
    </li>
    <li><a href="#functions">Functions</li>
      <ul>
        <li><a href="#data-preparation-functions">Data Preparation Functions</li>
        <ul>
          <li><a href="#fcn_findedge_loadlidardata">fcn_findEdge_loadLIDARData - Loads LiDAR and vehicle pose data</li>
          <li><a href="#fcn_findedge_pointsatrangeoflidarfromstation">fcn_findEdge_pointsAtRangeOfLiDARFromStation - Finds the scan lines located "range of LiDAR" meters from station 1 and station 2. </li>
          <li><a href="#fcn_findedge_finddrivenpathboundarypoints">fcn_findEdge_findDrivenPathBoundaryPoints - Find the boundary points of the driven path to create a bounding box for identifying the driven path grids.</li>
          <li><a href="#fcn_findedge_extractscanlines">fcn_findEdge_extractScanLines - Extracts vehicle pose ENU, vehicle pose unit orthogonal vectors, LiDAR ENU scans, LiDAR Intensity, LiDAR scan line, and Ring IDs of the scan line range</li>
          <li><a href="#fcn_findedge_finddrivablesurface">fcn_findEdge_findDrivableSurface - Finds the driven path points in LIDAR scans</li>
          <li><a href="#fcn_findedge_findpointsindomain">fcn_findEdge_findPointsInDomain - Finds the points in the domain from LiDAR ENU scans of the scan line range</li>
          <li><a href="#fcn_findedge_findmaxminofxyz">fcn_findEdge_findMaxMinOfXYZ - Finds the grid boundaries</li>
          <li><a href="#fcn_findedge_findgridswithpoints">fcn_findEdge_findGridsWithPoints - Separates the data into grids to find the empty and non-empty grids</li>
          <li><a href="#fcn_findedge_finddrivenpathgrids">fcn_findEdge_findDrivenPathGrids - Finds the driven path grids from the non-empty grids </li>
        </ul>
         <li><a href="#drivability-of-grids-functions">Drivability of Grids Functions</li>
        <ul>
          <li><a href="#fcn_findedge_determinegridpointdensity">fcn_findEdge_determineGridPointDensity - Determines the suitable "point density" for the analysis by comparing the point densities of driven grids with those of neighboring grids</li>
          <li><a href="#fcn_findedge_classifygridsbasedondensity">fcn_findEdge_classifyGridsBasedOnDensity - Classifies the grids based on the chosen point density</li>
          <li><a href="#fcn_findedge_calcnumberofgridscanlines">fcn_findEdge_calcNumberOfGridScanLines - Determines the number of LiDAR scan lines in each grid</li>
          <li><a href="#fcn_findedge_classifygridsbasedonscanlines">fcn_findEdge_classifyGridsBasedOnScanLines - Classifies the grids based on number of scan lines in each grid</li>
          <li><a href="#fcn_findedge_determinetransversespanthreshold">fcn_findEdge_determineTransverseSpanThreshold - Finds the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring. </li>
          <li><a href="#fcn_findedge_classifygridsbasedontransversespan">fcn_findEdge_classifyGridsBasedOnTransverseSpan - Classifies the grids based on chosen transverse span</li>
          <li><a href="#fcn_findedge_classifyqualifiedgrids">fcn_findEdge_classifyQualifiedGrids - Classifies the grids as qualified and unqualified</li>
          <li><a href="#fcn_findedge_finddrivenpathgrids">fcn_findEdge_findDrivenPathGrids - Finds the driven paths in the qualified grids</li>
        </ul>
        <li><a href="#grid-voting-functions">Grid Voting Functions</li>
        <ul>
          <li><a href="#fcn_findedge_determinestdinzerror">fcn_findEdge_determineSTDInZError - Finds the standard deviations of the z-fit errors for all qualified grids</li>
          <li><a href="#fcn_findedge_histogramstdinzerror">fcn_findEdge_histogramSTDinZError - Plots the histogram of standard deviations of qualified grids and driven path grids</li>
          <li><a href="#fcn_findedge_determineangledeviation">fcn_findEdge_determineAngleDeviation - Finds the angles between the unit normal vectors of the fitted planes and vertical ([0 0 1]) for all qualified grids</li>
          <li><a href="#fcn_findedge_histogramangledeviation">fcn_findEdge_histogramAngleDeviation - Plots the histogram of angle deviations of qualified grids and driven path grids</li>
          <li><a href="#fcn_findedge_classifygridsasdrivable">fcn_findEdge_classifyGridsAsDrivable - Classifies the qualified grids as drivable, non-drivable and uncertain</li>
        </ul>
        <li><a href="#post-processing-functions">Post Processing Functions</li>
        <ul>
          <li><a href="#fcn_findedge_prepgridcentersforboundarydetection">fcn_findEdge_prepGridCentersForBoundaryDetection - Prepares the grid centers for boundary detection</li>
          <li><a href="#fcn_findedge_findboundarypoints">fcn_findEdge_findBoundaryPoints - Find the boundary points of the grids based on the given grid centers</li>
          <li><a href="#fcn_findedge_findtrueboundarypoints">fcn_findEdge_findTrueBoundaryPoints - Finds the true boundary points by removing the boundary points of qualified and unqualified from drivable and non-drivable boundary points</li>
          <li><a href="#fcn_findedge_findnearestboundarypoints">fcn_findEdge_findNearestBoundaryPoints - Finds the nearest boundary points to the driven path </li>
          <li><a href="#fcn_findedge_seperateleftrightboundaries">fcn_findEdge_seperateLeftRightBoundaries - Separates the left and right boundaries from the nearest boundaries</li>
        </ul>
      </ul>
    <li><a href="#usage">Usage</a></li>
     <ul>
     <li><a href="#general-usage">General Usage</li>
     <li><a href="#examples">Examples</li>
     <li><a href="#definition-of-endpoints">Definition of Endpoints</li>
     </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

***

<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

The most common location of our testing is the Larson Test Track, and we regularly use “laps around the track” as replicates, hence the name of the library. And when not on the test track and on public roads, data often needs to be segmented from one keypoint to another. For example, it is a common task to seek a subset of path data that resides only from one intersection to the next. While one could segment this data during data collection by simply stopping the vehicle recordings at each segment, it is impractical and dangerous to stop data collection at each and every possible intersection or feature point. Rather, vehicle or robot data is often collected by repeated driving of an area over/over without stopping. So, the final data set may contain many replicates of the area of interest.

This "Laps" code assists in breaking recorded path data into paths by defining specific start and end locations, for example from intersection "A" to stop sign "B". Specifically, the purpose of this code is to break data into "laps", e.g. segments of data that are defined by a clear start condition and end condition. The code finds when a given path meets the "start" condition, then meets the "end" condition, and returns every portion of the path that is inside both conditions. There are many advanced features as well including the ability to define excursion points, the number of points that must be within a condition for it to activate, and the ability to extract the portions of the paths before and after each lap, in addition to the data for each lap.

* Inputs:
  * either a "traversals" type, as explained in the Path library, or a path of XY points in N x 2 format
  * the start, end, and optional excursions can be entered as either a line segment or a point and radius.  
* Outputs
  * Separate arrays of XY points, or of indices for the lap, with one array for each lap
  * The function also can return the points that were not used for laps, e.g. the points before the first start and after the last end

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo

   ```sh
   git clone https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps
   ```

3. Run the main code in the root of the folder (script_demo_Laps.m), this will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory and clear all global variables in MATLAB (type: "clear global *").
4. Confirm it works! Run script_demo_Laps. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->
### Directories

The following are the top level directories within the repository:
<ul>
 <li>/Documents folder: Descriptions of the functionality and usage of the various MATLAB functions and scripts in the repository.</li>
 <li>/Functions folder: The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
 <li>/Utilities folder: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>
</ul>

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

* [PathPlanning_PathTools_PathClassLibrary](https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary) - the PathClassLibrary contains tools used to find intersections of the data with particular line segments, which is used to find start/end/excursion locations in the functions. The repo can be found at: <https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary>

    Each should be installed in a folder called "Utilities" under the root folder, namely ./Utilities/DebugTools/ , ./Utilities/PathClassLibrary/ . If you wish to put these codes in different directories, the main call stack in script_demo_Laps can be easily modified with strings specifying the different location, but the user will have to make these edits directly.

    For ease of getting started, the zip files of the directories used - without the .git repo information, to keep them small - are included in this repo.

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- FUNCTION DEFINITIONS -->
## Functions

### Data Preparation Functions

#### fcn_findEdge_loadLIDARData

Loads LIDAR data from a specified file. Peforms loading of vehicle pose
and LIDAR data with many conditional checks to confirm the file is
present, the variables were not already loaded, the files are up-to-date,
etc.

```MATLAB
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num), (fig_num2));
```
<pre align="center">
  <img src=".\Images\fcn_findEdge_loadLIDARData_1.jpg" alt="fcn_findEdge_loadLIDARData picture" width="500" height="400">
  <figcaption>Fig.1 - Vehicle Trajectory</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_loadLIDARData_2.jpg" alt="fcn_findEdge_loadLIDARData picture" width="500" height="400">
  <figcaption>Fig.2 - Vehicle Trajectory in 2D</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_pointsAtRangeOfLiDARFromStation

This code finds the index of a point that is a certain range before
station 1 and the index of a point that is a certain range after
station 2. 

```MATLAB
[station1_minus_range_index, station2_plus_range_index]...
= fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,...
starting_index,ending_index,(range))
```

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findDrivenPathBoundaryPoints

Find the boundary points of the driven path to create a bounding box for finding the driven path grids 

```MATLAB
fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose,...
 scanLineRange, Nscans, shift, (fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathBoundaryPoints_1.jpg" alt="fcn_findEdge_findDrivenPathBoundaryPoints picture" width="500" height="400">
  <figcaption>Fig.3 - Boundary Points of Driven Path</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>


<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathBoundaryPoints_2.jpg" alt="fcn_findEdge_findDrivenPathBoundaryPoints picture" width="500" height="400">
  <figcaption>Fig.4 - Left and Right Boundary Lane</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_extractScanLines

 Extracts vehicle pose ENU, vehicle pose unit orthogonal vectors,
 LiDAR ENU scans, LiDAR Intensity, LiDAR scan line, and Ring IDs of the scan line range

```MATLAB
     [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
     LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
     fcn_findEdge_extractScanLines(VehiclePose,
     LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (fig_num), (fig_num2))
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_extractScanLines_1.jpg" alt="fcn_findEdge_extractScanLines picture" width="500" height="400">
  <figcaption>Fig.5 - LIDAR Data in 3D ENU </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_extractScanLines_2.jpg" alt="fcn_findEdge_extractScanLines picture" width="500" height="400">
  <figcaption>Fig.6 - Vehicle Trace in LLA </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findDrivableSurface

 Find the driven path points in LIDAR scans

```MATLAB
     [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface ...
     (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,(fig_num),(fig_num2))
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivableSurface_1.jpg" alt="fcn_findEdge_findDrivableSurface picture" width="500" height="400">
  <figcaption>Fig.7 - LIDAR Data Underneath the Vehicle in ENU</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivableSurface_2.jpg" alt="fcn_findEdge_findDrivableSurface picture" width="500" height="400">
  <figcaption>Fig.7 - LIDAR Data Underneath the Vehicle in LLA</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>
<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findPointsInDomain

Find the points in the domain from LiDAR ENU scans of the scan line range

```MATLAB
     [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,(fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_findPointsInDomain.jpg" alt="fcn_findEdge_findPointsInDomain picture" width="500" height="400">
  <figcaption>Fig.8 - LIDAR Points, Vehicle Pose and Boundary Points in LLA</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findMaxMinOfXYZ

Find the maximum and minimum numbers of x, y, and z of the given points 

```MATLAB
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(N_points,(fig_num))
```   
<pre align="center">
  <img src=".\Images\fcn_findEdge_findMaxMinOfXYZ.jpg" alt="fcn_findEdge_findMaxMinOfXYZ picture" width="500" height="400">
  <figcaption>Fig.9 - ENU Trace Geometry </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findGridsWithPoints

Find the grid boundaries and separate the data into grids to find the empty and non-empty grids

```MATLAB
    [gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_point,
     grids_greater_than_zero_points, gridCenters_zero_point_density,
     gridCenters_greater_than_zero_point_density] = fcn_findEdge_findGridsWithPoints(input_points,
     grid_size,grid_boundaries,(fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_findGridsWithPoints.jpg" alt="fcn_findEdge_findGridsWithPoints picture" width="500" height="400">
  <figcaption>Fig.10 - Grids of LIDAR Points </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findDrivenPathGrids

Find the driven path grids from the non-empty grids 

```MATLAB
    [total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,grids_greater_than_zero_points, (fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathGrids_1.jpg" alt="fcn_findEdge_findDrivenPathGrids picture" width="500" height="400">
  <figcaption>Fig.11 - ENU XY Plot of Vehicle Trajectory </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathGrids_2.jpg" alt="fcn_findEdge_findDrivenPathGrids picture" width="500" height="400">
  <figcaption>Fig.12 - Boundary Points in LLA  </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

### Drivability of Grids Functions

#### fcn_findEdge_determineGridPointDensity

Determine the suitable "point density" for the analysis by comparing the point densities of driven grids with those of neighboring grids

```MATLAB
    [point_density] = fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero ...
    ,total_points_in_each_grid_in_the_driven_path,grid_size,(N_bins_grid_with_points_greater_than_zero)...
    (N_bins_grid_in_the_driven_path),(fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_determineGridPointDensity.jpg" alt="fcn_findEdge_determineGridPointDensity picture" width="500" height="400">
  <figcaption>Fig.13 - Grids Density</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>


<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_classifyGridsBasedOnDensity

Determine the suitable "point density" for the analysis by comparing the point densities of driven grids with those of neighboring grids

 ```MATLAB
     [grid_indices_with_required_point_density, gridCenters_low_point_density] =  fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points...
     ,total_N_points_in_each_grid,point_density, gridCenters, (format_1),(format_2),
     (fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_classifyGridsBasedOnDensity.jpg" alt="fcn_findEdge_classifyGridsBasedOnDensity picture" width="500" height="400">
  <figcaption>Fig.14 - </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_calcNumberOfGridScanLines

Determine the number of LiDAR scan lines in each grid

 ```MATLAB
     [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points)
```  


<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_classifyGridsBasedOnScanLines

Determine the suitable "point density" for the analysis by comparing the point densities of driven grids with those of neighboring grids

 ```MATLAB
       [grid_indices_with_more_than_one_scan_line, gridCenters_with_one_scan_line] = fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points...
       ,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters, (format_1),(format_2), (fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_classifyGridsBasedOnScanLines.jpg" alt="fcn_findEdge_classifyGridsBasedOnScanLines picture" width="500" height="400">
  <figcaption>Fig.15 - Grid Center with One Scan Line</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_determineTransverseSpanThreshold

Find the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring. 

 ```MATLAB
      fcn_findEdge_determineTransverseSpanThreshold:
      [transverse_span_threshold,transverse_span_each_grid]= 
      fcn_findEdge_determineTransverseSpanThreshold(grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines, fig_num, fig_num2, fig_num3)
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_determineTransverseSpanThreshold_1.jpg" alt="fcn_findEdge_determineTransverseSpanThreshold picture" width="500" height="400">
  <figcaption>Fig.16 - Gridlines and Points of the Grids Greater than Zero Point Density </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineTransverseSpanThreshold_2.jpg" alt="fcn_findEdge_determineTransverseSpanThreshold picture" width="500" height="400">
  <figcaption>Fig.17 - Grid Lines of the Grids Greater than Zero </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineTransverseSpanThreshold_3.jpg" alt="fcn_findEdge_determineTransverseSpanThreshold picture" width="500" height="400">
  <figcaption>Fig.18 - Grid Lines of the Grids Greater than Zero </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_classifyGridsBasedOnTransverseSpan

Find the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring. 

 ```MATLAB
      [grid_indices_with_more_than_transverse_span_threshold, gridCenters_with_less_than_transverse_span_threshold] =
      fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points, gridCenters, (format_1),(format_2),(fig_num))
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_classifyGridsBasedOnTransverseSpan.jpg" alt="fcn_findEdge_classifyGridsBasedOnTransverseSpan picture" width="500" height="400">
  <figcaption>Fig.19 - Grid with Threshold </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_classifyQualifiedGrids

If either of the conditions fail, the non-empty grids are classified as UNQUALIFIED. If all the conditions pass, the non-empty grids are classified as QUALIFIED.

 ```MATLAB
      [original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
      fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,...
      grid_indices_with_more_than_transverse_span_threshold, grids_greater_than_zero_points,gridCenters,(format_unqualified),(format_qualified),(fig_num))
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_classifyQualifiedGrids.jpg" alt="fcn_findEdge_classifyQualifiedGrids picture" width="500" height="400">
  <figcaption>Fig.20 - Unqualified and Qualified Grids </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findDrivenPathGrids

Recalculate the driven path grids among qualified grids

 ```MATLAB
      [total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
     = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
     grids_greater_than_zero_points, (fig_num))
```
<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathGrids_3.jpg" alt="fcn_findEdge_findDrivenPathGrids picture" width="500" height="400">
  <figcaption>Fig.21 - New ENU XY Plot of Vehicle Trajectory </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_findDrivenPathGrids_4.jpg" alt="fcn_findEdge_findDrivenPathGrids picture" width="500" height="400">
  <figcaption>Fig.21 - New Boundary Points in LLA </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

### Grid Voting Functions

#### fcn_findEdge_determineSTDInZError

Find the standard deviations of the z-fit errors for all qualified grids, and compare these standard deviations with those of the grids in the driven path to determine the suitable standard deviation threshold

 ```MATLAB
      [input_points,original_mapped_gridIndices_cell,total_mapped_grids,total_points_in_mapped_grids,standard_deviation_in_z,gridlines_mapped_grids,...
 driven_path_grid_indices_in_current_mapped_grids,std_in_z_driven_path,std_in_z_other_mapped_grids,mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path]...
              = fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,...
              current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,(fig_num_gridCenters),(fig_num_gridSTD),(fig_num))
```
<pre align="center">
  <img src=".\Images\fcn_findEdge_determineSTDInZError_1.jpg" alt="fcn_findEdge_determineSTDInZError picture" width="500" height="400">
  <figcaption>Fig.22 - Grid centers mapped grids in ENU </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineSTDInZError_2.jpg" alt="fcn_findEdge_determineSTDInZError picture" width="500" height="400">
  <figcaption>Fig.23 - Standard deviation of Z of mapped grids </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineSTDInZError_3.jpg" alt="fcn_findEdge_determineSTDInZError picture" width="500" height="400">
  <figcaption>Fig.24 - Qualified grid centers vs standard deviation in Z </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_histogramSTDinZError

Histogram of standard deviations of the z-fit errors for all qualified grids

 ```MATLAB
      fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))
```
<pre align="center">
  <img src=".\Images\fcn_findEdge_histogramSTDinZError.jpg" alt="fcn_findEdge_histogramSTDinZError picture" width="500" height="400">
  <figcaption>Fig.25 - Standard Deviation in Z </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_determineAngleDeviation

Find the angles between the unit normal vectors of the fitted planes and vertical ([0 0 1]) for all qualified grids, and compare these angles with those of the grids in the driven path to determine the suitable theta threshold

 ```MATLAB
      [angle_btw_unit_normals_and_vertical, mean_angle_btw_unit_normals_and_vertical_driven_path] = ...
      fcn_findEdge_determineAngleDeviation(LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
      gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path,(fig_num),(fig_num2),(fig_num3) )
```
<pre align="center">
  <img src=".\Images\fcn_findEdge_determineAngleDeviation_1.jpg" alt="fcn_findEdge_determineAngleDeviation picture" width="500" height="400">
  <figcaption>Fig.26 - Grid Centers Mapped Grids in ENU </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineAngleDeviation_2.jpg" alt="fcn_findEdge_determineAngleDeviation picture" width="500" height="400">
  <figcaption>Fig.27 - Angle Deviation of Mapped Grids</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_findEdge_determineAngleDeviation_3.jpg" alt="fcn_findEdge_determineAngleDeviation picture" width="500" height="400">
  <figcaption>Fig.28 - Mapped Grid Centers vs Standard Deviation in Z</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_histogramAngleDeviation

Find the angles between the unit normal vectors of the fitted planes and vertical ([0 0 1]) for all qualified grids, and compare these angles with those of the grids in the driven path to determine the suitable theta threshold

 ```MATLAB
      theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path, fig_num);
``` 
<pre align="center">
  <img src=".\Images\fcn_findEdge_histogramAngleDeviation.jpg" alt="fcn_findEdge_histogramAngleDeviation picture" width="500" height="400">
  <figcaption>Fig.29 - Histogram of angle deviation </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_classifyGridsAsDrivable

If all conditions pass, DRIVABLE
If all conditions fail, NON_DRIVABLE
If only a few fail, UNCERTAIN

 ```MATLAB
      [grid_indices_with_required_point_density, gridCenters_low_point_density] =  fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points...
     ,total_N_points_in_each_grid,point_density, gridCenters, (format_1),(format_2),
     (fig_num))
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_classifyGridsAsDrivable.jpg" alt="fcn_findEdge_classifyGridsAsDrivable picture" width="500" height="400">
  <figcaption>Fig.30 - Points Density</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

### Post Processing Functions

#### fcn_findEdge_prepGridCentersForBoundaryDetection

Prepare the grid centers of qualified and unqualified for boundary detection

 ```MATLAB
      [X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);
``` 

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findBoundaryPoints

Find the boundary points of qualified and unqualified grids

```MATLAB
      boundary_points_qualified_unqualified = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_qualified_unqualified);;
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_findBoundaryPoints.jpg" alt="fcn_findEdge_findBoundaryPoints picture" width="500" height="400">
  <figcaption>Fig.31 - Boundary Points </figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findTrueBoundaryPoints

Find the true boundary points by removing the boundary points of qualified and unqualified from drivable and non-drivable boundary points

```MATLAB
      [true_boundary_points] = fcn_findEdge_findTrueBoundaryPoints(boundary_points,boundary_points_qualified_unqualified,fig_num_bd_pts_ENU);
```  
<pre align="center">
  <img src=".\Images\fcn_findEdge_findTrueBoundaryPoints.jpg" alt="fcn_findEdge_findTrueBoundaryPoints picture" width="500" height="400">
  <figcaption>Fig.32 - Computed Boundary Points</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_findNearestBoundaryPoints

Find the nearest boundary points to the driven path 

```MATLAB
    [true_borders] = fcn_findEdge_findNearestBoundaryPoints(boundaryPointsXY,
     gridCenters_non_drivable_grids, gridCenters_driven_path, (fig _num))
```   
<pre align="center">
  <img src=".\Images\fcn_findEdge_findNearestBoundaryPoints.jpg" alt="fcn_findEdge_findNearestBoundaryPoints picture" width="500" height="400">
  <figcaption>Fig.33 - Nearest Boundary Points</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***

#### fcn_findEdge_seperateLeftRightBoundaries

Separate the left and right boundaries from the nearest boundaries

```MATLAB
  [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
   (VehiclePose,station_1,station_2,nearest_boundary_points, grid_size,
   transverse_shift, (fig_num)).
``` 
 
<pre align="center">
  <img src=".\Images\fcn_findEdge_seperateLeftRightBoundaries.jpg" alt="fcn_findEdge_seperateLeftRightBoundaries picture" width="500" height="400">
  <figcaption>Fig.34 - Left and Right Boundary</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#featureextraction_laneboundary_findedge">Back to top</a>

***


<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

### General Usage

Each of the functions has an associated test script, using the convention

```sh
script_test_fcn_fcnname
```

where fcnname is the function name as listed above.

As well, each of the functions includes a well-documented header that explains inputs and outputs. These are supported by MATLAB's help style so that one can type:

```sh
help fcn_fcnname
```

for any function to view function details.

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

### Examples

1. Run the main script to set up the workspace and demonstrate main outputs, including the figures included here:

   ```sh
   script_demo_Laps
   ```

    This exercises the main function of this code: fcn_Laps_breakDataIntoLaps

2. After running the main script to define the included directories for utility functions, one can then navigate to the Functions directory and run any of the functions or scripts there as well. All functions for this library are found in the Functions sub-folder, and each has an associated test script. Run any of the various test scripts, such as:

   ```sh
   script_test_fcn_Laps_breakDataIntoLapIndices
   ```

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

### Definition of Endpoints

The codeset uses two types of zone definitions:

1. A point location defined by the center and radius of the zone, and number of points that must be within this zone. An example of this would be "travel from home" or "to grandma's house". The point "zone" specification is given by an X,Y center location and a radius in the form of [X Y radius], as a 3x1 matrix. Whenever the path passes within the radius with a specified number of points within that radius, the minimum distance point then "triggers" the zone.

    <img src=".\Images\point_zone_definition.png" alt="point_zone_definition picture" width="200" height="200">

2. A line segment. An example is the start line or finish line of a race. A runner has not started or ended the race without crossing these lines. For line segment conditions, the inputs are condition formatted as: [X_start Y_start; X_end Y_end] wherein start denotes the starting coordinate of the line segment, end denotes the ending coordinate of the line segment. The direction of start/end lines of the segment are defined such that a correct crossing of the line is in the positive cross-product direction defined from the vector from start to end of the segment.

    <img src=".\Images\linesegment_zone_definition.png" alt="linesegment_zone_definition picture" width="200" height="200">

These two conditions can be mixed and matched, so that one could, for example, find every lap of data where someone went from a race start line (defined by a line segment) to a specific mountain peak defined by a point and radius.

The two zone types above can be used to define three types of conditions:

1. A start condition - where a lap starts. The lap does not end until and end condition is met.
2. An end condition - where a lap ends. The lap cannot end until this condition is met.
3. An excursion condition (optional) - a condition that must be met after the start point, and before the end point. The excursion condition must be met before the end point is counted.

Why is an excursion point needed? Consider an example: it is common for the start line of a marathon to be quite close to the start line, sometimes even just a few hundred feet after the start line. This setup is for the practical reason that runners do not want to make long walks to/from starting locations to finish location either before, and definitely not after, such a race. As a consequence, it is common that, immediately after the start of the race, a runner will cross the finish line before actually finishing the race. This happens in field data collection when one accidentally passes a start/end station, and then backs up the vehicle to reset. In using these data recordings, we would not want these small segment to count as a complete laps, for example the 100-ish meter distance to be counted as a marathon run. Rather, one would require that the recorded data enter some excursion zone far away from the starting line for such a "lap" to count. Thus, this laps code allows one to define an excursion point as a location far out into the course that one must "hit" before the finish line is counted as the actual "finish" of the lap.

* For each lap when there are repeats, the resulting laps of data include the lead-in and fade-out data, namely the datapoint immediately before the start condition was met, and the datapoint after the end condition is met. THIS CREATES REPLICATE DATA. However, this allows better merging of data for repeated laps, for example averaging data exactly from start to finish, or to more exactly calculate velocities on entry and exit of a lap by using windowed averages or filters.

* Points inside the lap can be set for the point-type zones. These occur as optional input arguments in fcn_Laps_findPointZoneStartStopAndMinimum and in the core definition of a point zone as the 2nd argument. For example, the following code:

  ```Matlab
  start_definition = [10 3 0 0]; % Radius 10, 3 points must pass near [0 0]
  ```

  requires 3 points to occur within the start zone area.

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

## Major release versions

This code is still in development (alpha testing)

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [hhttps://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps](https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps)

<a href="#featureextraction_dataclean_breakdataintolaps">Back to top</a>

***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
