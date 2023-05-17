# PathPlanning_PathTools_PathClassLibrary

This repo supports mathematical operations for paths and trajectories.
<!-- PROJECT LOGO -->
<br>
<p align="center">
  <h2 align="center">
    PathPlanning_PathTools_PathClassLibrary
  </h2>
  <pre align="center">
        <img src=".\Images\TwoRoadsDivergeInForest_small.jpg" alt="main pathtools picture" width="960" height="540">
        <!-- figcaption>Fig.1 - The typical progression of map generation.</figcaption -->
        <font size="-2">Photo by <a href="https://unsplash.com/@madebyjens?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Jens Lelie</a> on <a href="https://unsplash.com/photos/u0vgcIOQG08?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a>
    </font>
  </pre>
</p>

<p align="left">
    This is the Path class library repo for MATLAB and listing of all functions within this class, whose purpose is to abstract very common "path" operations. What's a "path"? Why do we need them? How do we use them? Read more below for some exciting tools!
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About the Project</a>
      <ul>
        <li><a href="#before-you-start">Before you start</a></li>
      </ul>
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
    <li><a href="#functions">Functions</li>
        <ul>
        <li><a href="#motivation-and-definitions">Motivation and Definitions</li>
          <ul>
            <li><a href="#definition-of-station">Definition of station</li>
            <li><a href="#definition-of-stations">Definition of stations</li>
            <li><a href="#definition-of-path">Definition of path</li>
            <li><a href="#definition-of-paths">Definition of paths</li>
            <li><a href="#definition-of-traversal">Definition of traversal</li>
            <li><a href="#definition-of-traversals">Definition of traversals</li>
          </ul>
        <li><a href="#basic-getting-started-functions">Basic Getting Started Functions</li>
          <ul>
            <li><a href="#fcn_path_checkinputstofunctions">fcn_Path_checkInputsToFunctions - a general-use input checking function (deprecated by the Debug library, so replaced soon).</li>
            <li><a href="#fcn_path_fillsamplepaths">fcn_Path_fillSamplePaths - Fills test data and includes a click-to-draw.</li>
            <li><a href="#fcn_path_calcdiffanglesbetweenpathsegments">fcn_Path_calcDiffAnglesBetweenPathSegments - Calculates the change in angles between path segments.</li>
            <li><a href="#fcn_path_calcyawfrompathsegments">fcn_Path_calcYawFromPathSegments - Calculates yaw angle for each path segment.</li>
            <li><a href="#fcn_path_convertpathtotraversalstructure">fcn_Path_convertPathToTraversalStructure - Converts paths into traversals.</li>
            <li><a href="#fcn_path_plottraversalsxy">fcn_Path_plotTraversalsXY - Plots the XY value of the traversals.</li>
            <li><a href="#fcn_path_plottraversalsyaw">fcn_Path_plotTraversalsYaw - Plots the yaw values of the traversals.</li>
          </ul>
        <li><a href="#basic-path-operations">Basic Path Operations</li>
          <ul>
            <li><a href="#fcn_path_findprojectionhitontopath">fcn_Path_findProjectionHitOntoPath - Finds intersection between a sensor and a path.</li>
            <li><a href="#fcn_path_findintersectionsbetweentraversals">fcn_Path_findIntersectionsBetweenTraversals - Finds where two traversals intersect.</li>
            <li><a href="#fcn_path_snappointontonearestpath">fcn_Path_snapPointOntoNearestPath - snaps a point onto the nearest path segment.</li>
            <li><a href="#fcn_path_snappointontonearesttraversal">fcn_Path_snapPointOntoNearestTraversal - snaps point to nearest traversal.</li>
          </ul>
        <li><a href="#functions-that-trim-paths">Functions That Trim Paths</li>
          <ul>
            <li><a href="#fcn_path_findtraversalstationsegment">fcn_Path_findTraversalStationSegment - crops traversal by given station interval.</li>
            <li><a href="#fcn_path_removepinchpointintraversal">fcn_Path_removePinchPointInTraversal - Eliminates self-crossings in traversals.</li>
          </ul>
        <li><a href="#functions-that-project-paths">Functions That Project Paths</li>
          <ul>
            <li><a href="#fcn_path_findorthogonaltraversalvectorsatstations">fcn_Path_findOrthogonalTraversalVectorsAtStations - calculates orthogonal vectors to a traersal at given stations.</li>
            <li><a href="#fcn_path_findorthogonalhitfromtraversaltotraversal">fcn_Path_findOrthogonalHitFromTraversalToTraversal - finds which traversals are hit at ortho projections from one traversal to another.</li>
            <li><a href="#fcn_path_findorthoscatterfromtraversaltotraversals">fcn_Path_findOrthoScatterFromTraversalToTraversals - finds closest points on many traversals to a given central traversal. </li>
            <li><a href="#fcn_path_filloffsettraversalsabouttraversal">fcn_Path_fillOffsetTraversalsAboutTraversal - fills in an array of traversals about a reference traversal at user-defined offset distances.</li>
            <li><a href="#fcn_path_converttraversalxytosy">fcn_Path_convertTraversalXYtoSy - calculates the SY (e.g ST or "station") coordinates for a given traversal, given a reference traversal and another query traversal.</li>
            <li><a href="#fcn_path_fillrandomtraversalsabouttraversal">fcn_Path_fillRandomTraversalsAboutTraversal - generates random traversals about a given traversal.</li>
          </ul>
        <li><a href="#statistical-analysis-of-paths">Statistical Analysis of Paths</li>
          <ul>
            <li><a href="#fcn_path_plottraversalxywithupperlowerbands">fcn_Path_plotTraversalXYWithUpperLowerBands - plots a central traversal with a band around it defined by an upper and lower traversal. All traversals must have same length.</li>
            <li><a href="#fcn_path_calcsingletraversalstandarddeviation">fcn_Path_calcSingleTraversalStandardDeviation - estimates the standard deviation in the offsets of a single traversal using variance in path angles and average path segment length.</li>
            <li><a href="#fcn_path_plottraversalxywithvariancebands">fcn_Path_plotTraversalXYWithVarianceBands - plots a traversal with standard deviation boundaries around it.</li>
          </ul>
        <li><a href="#path-averaging-methods">Path Averaging Methods</li>
          <ul>
            <li><a href="#fcn_path_findaveragetraversalviastationalignment">fcn_Path_findAverageTraversalViaStationAlignment - estimates an average traversal based on aligning traversals by station.</li>
            <li><a href="#fcn_path_findaveragetraversalviaclosestpoint">fcn_Path_findAverageTraversalViaClosestPoint - estimates an average traversal based on closest point.</li>
            <li><a href="#fcn_path_newtraversalbystationresampling">fcn_Path_newTraversalByStationResampling - redecimates a traversal at given station points.</li>
            <li><a href="#fcn_path_findaveragetraversalviaorthoprojection">fcn_Path_findAverageTraversalViaOrthoProjection - estimates average traversal using orthogonal projection.</li>
            <li><a href="#fcn_path_findtraversalwithmostdata">fcn_Path_findTraversalWithMostData - finds the traversal in a traversals array with the most data (most elements in X array).</li>
            <li><a href="#fcn_path_cleanpathfromforwardbackwardjogs">fcn_Path_cleanPathFromForwardBackwardJogs - removes back/forth jogs from data.</li>
            <li><a href="#fcn_path_findclosestpointstotraversal">fcn_Path_findClosestPointsToTraversal - finds closest points on a set of traversals to reference traversal.</li>
          </ul>
        <li><a href="#3d-and-elevated-paths">3D and Elevated Paths</li>
          <ul>
            <li><a href="#fcn_path_addelevationtopath">fcn_Path_addElevationToPath - fills in the elevation in a Path's XY data using a reference traversal and interpolation.</li>
          </ul>
        <li><a href="#miscellaneous-functions">Miscellaneous Functions</li>
        <ul>
        </ul>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
     <ul>
     <li><a href="#general-usage">General Usage</li>
     <li><a href="#examples">Examples</li>
     </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- ABOUT THE PROJECT -->
## About The Project

This "path" class library is one of the first and most commonly used libraries developed within the Intelligent Vehicles and Systems Group, as it contains the set of functions for processing directionality operations related to the "paths" that vehicles follow. This includes calculating intersections, offsets, projections, etc. For example: do you want to know the closest point a vehicle travels to an point in space (e.g. path to point)? Use this library. Or: what is the nearest road centerline to the current vehicle location (e.g. point to path)? Use this library. Or, figure out when a vehicle collides with another vehicle (intersection)? Use this library! Or even, what is the average trajectory of many different vehicles going down the road? You guessed it - use this library.

Most of the tools used herein are abstracted from video gaming algorithms, and thus are highly optimized for speed in an algorithm context, e.g. they are parallelizable, support vectorization of inputs, minimize internal looping, and are highly robust to user input errors and fault cases. Where possible, we try to use matrix forms. All codes are in MATLAB to support generalization into other languages (Python, C, C++, etc.) and some students have started this conversion as these codes are FAR faster. However, MATLAB still seems to be the best "scratchpad" to test algorithms against each other (with Python being great also) - and if an algorithm is poorly chosen, a fast programming language will rarely overcome this flaw.

### Before you start

To learn more about "paths" before starting, please watch the excellent YouTube video, "The Continuity of Splines", by Freya Holmer (her videos are amazing, among the best out there). Her introduction to splines can be found at this link: <https://www.youtube.com/watch?v=jvPPXbo87ds&ab_channel=FreyaHolm%C3%A9r>.

Chapter 1 of her video introduces the blending "t" notation of interpolation (used throughout the functions in "Paths"), the "lerp" operation (e.g. linear interpolation), Bezier curve derivations as a function of lerp in quadratic and cubic forms via DeCasteljau's algorithm. It also introduces variants such as the Bernstein form using basis polynomials (which follows Pascal's triangle, see <https://math.stackexchange.com/questions/447718/algorithm-for-bezier-curve-with-eleven-control-points>) and is based on convex combinations. It also shows how to pre-calculate polynomial coefficients for Bezier curves, and the characteristic matrix form of the Bezier curve (which, by the way, is invertable!). It also introduces the problem of local versus global control of a path and why global control is problematic.

Chapter 2 (at 10 minutes in) introduces Bezier splines, definition of u-values and curve indices, the relationships between parameter space, control points, splines, joins and knots, knot values and knot intervals. It also discussed broken tangents, aligned tangents, and mirrored tangents.

Chapter 3 introduces parametric continuity and measures of connectivity including C0 continuity, C1 continuity, etc. and how to obtain these via derivative functions (velocity, acceleration, jerk, snap, crackle, pop). It also shows how the addition of constraints drastically affects the sensitivity of a trajectory to perturbations.

Chapter 4 introduces geometric continuity, the tangent, and normal vectors. It introduces G1 and higher continuity. The issue with G1 seam continuity shows up in curvature calculations for vehicle dynamics. If a path has a G1 discontinuity, then the steering angle has to change instantly (which is not possible) for a vehicle to stay on the road. Namely, the osculating circle must change continuously. It also introduces the calculation of curvature, $\kappa$ as:

$$ \kappa = \frac{
    \begin{vmatrix}
    P^{'}_x & P^{''}_x \\
    \\
    P^{'}_y & P^{''}_y \\
    \end{vmatrix}
    }{||P^{'}||^{3}} $$

where $P^{'}_x$ is the velocity in x, etc. and the determinant is the 2D cross product (wedge product producing a bivector). Note: curvature is increadibly important in vehicle dynamics because it directly relates to steering angle! This chapter also introduces the idea of a curvature comb, e.g. the plot of the curvature normal to the curve itself. Note that all roads are G2 or higher continuous, e.g. all roads are class-A surfaces. Further, roads that are G3 continous will - for non-slip conditions of the tire-road boundary - have smooth steering inputs to stay on the road. Most modern roads are G2 or higher continuous. (Note: this makes a HARD problem of collecting road data, and fitting the G2 or G3 curves to the data!). It also introduces the regularity condition, which shows up in vehicle dynamics as well where the bicycle model representation of vehicle dynamics also requires regularity.

Chapter 5 introduces the problem of splines as a boundary value problem, which is useful for calculating trajectories for vehicle behaviors. It also introduces the Hermite spline as the extension of Bezier spline. Most of the Path library deals with the Linear spline representation, which is also introduced in this chapter. We do not use cardinal splines but they do appear in racing lines. It also introduces the Catmull-Rom spline. It then derives the B-spline representation from constriants on the Catmull-Rom, and wraps up the video there.

Most of the issues discussed in this video sneak into the code in this library (as well as MANY more).

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo

   ```sh
   git clone https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
   ```

3. Run the main code in the root of the folder (script_demo_AlignCoordinates.m), this will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory and clear all global variables in MATLAB (type: "clear global *").

4. Confirm it works! Run script_demo_AlignCoordinates. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->
### Directories

The following are the top level directories within the repository:
<ul>
 <li>/Data folder: contains any example datasets used in the code (not applicable to this repo).</li>
 <li>/Documents folder: contains reference documents used for the creation of the code.</li>
 <li>/Functions folder: The majority of the codes are found within functions in this directory. All functions as well as test scripts are provided.</li>
 <li>/Images folder: images used for this README.md file are found in this directory.</li>
 <li>/Releases folder: contains current and historic releases of the code in zip file archives.</li>
 <li>/Utilities folder: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>
</ul>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

The dependencies are automatically installed by running the root master script (see below).

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- FUNCTION DEFINITIONS -->
## Functions

### Motivation and Definitions

Many applications require a reference path, but particularly driving. For example, roads are a form of path. The task of driving is to stay on the road, which means that a vehicle's position relative to the path - the offset to the path, travel along the path, angle relative to path, etc. -  are all basically path calculations.

<pre align="center">
  <img src=".\Images\RoadsAreAFormOfPath.png" alt="roads are a form of path picture" width="400" height="300">
  <figcaption>Roads are a type of path.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Definitions of terms are important so that we do not get confused later. These key definitions are defined in the following sections.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Station

The station of a path is the location on that path as measured by the s-coordinate of that path. Physically, it represents the distance traveled along that path from the path's origin.

For example, mile markers on a road are a station type of measurement; they do not represent distance to a specific point, but rather represent how far along a road someone has traveled relative to some start location and independent of the curvature of the road. In this case, the station is "39 miles". For our work, we nearly always measure station in meters.

<pre align="center">
  <img src=".\Images\DefinitionOfStation.jpg" alt="definition of station" width="400" height="300">
  <figcaption>Mile markers are a type of station measurement.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The important aspect of the station is that it does NOT depend on Euclidian distance from a location. The station distance is distance ALONG a path, so if a path is curvy, then the station distance from a point can be very large even if the point is relatively close in straight-line distance.In fact, for curvy paths, the Euclidean distance may be rise, fall, and then rise again as vehicles drive "away" from a point!

<pre align="center">
  <img src=".\Images\DefinitionOfStation2_small.jpg" alt="definition of station" width="900" height="300">
  <figcaption>Station coordinates are not Euclidian and so true distances from a point, for example the bottom left of the image above, may increase and then decrease and increase, even while station distance is steadily increasing.</figcaption>
  <font size="-2">Photo by <a href="https://unsplash.com/de/@umityildirim?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Umit Yildirim</a> on <a href="https://unsplash.com/photos/mY8-8Fsn3qw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a>
  </font>
</pre>

There are several key issues when using station coordinates. The first is that, for real-world paths, the station length of a path will change depending on how finely decimated the path is. For very curvy paths, the station distance can tend toward infinity with finer and finer scales of measurement as some paths follow fractacl curves. This is called the "Coastline Paradox"; see, for example, <a href="https://www.youtube.com/watch?v=7dcDuVyzb8Y"> Measuring Coastline - Numberphile </a>. A solution to this is to seek Allan Variance minima in measurements of time or space, a process covered in other repos. The idea of a distance that is not quite distance but not quite area falls into the class of problems called "intermediate asymptotics", "self similarity", and "fractal dimension" and is a difficult challenge.

<pre align="center">
  <img src=".\Images\CoastlineParadox.png" alt="definition of the coastline paradox" width="600" height="800">
  <figcaption>An illustration of the Coastline Paradox</figcaption>
  <font size="-2">Photo from <a href="https://sketchplanations.com/the-coastline-paradox">Sketchplanations.com</a>
  </font>
</pre>

A practical solution, and one used for most vehicle calculations, is to choose a decimation length that is fine enough that the measurements of station distance converge to the desired accuracy, for example 1 meter, 10 cm, or 1 cm depending on the required accuracy. An easy method to approximate the error of a path definition is to take the change in angle of a path segment to another, and multiply by the distance of the path segment. Over the whole path, this can be done as maximum angle times maximum length, or average angle times average length. The statistical function listed below does this calculation.

Even with a well-defined reference path, the length of other paths along the reference path depends very strongly on how much weaving is occurring and where the other paths are relative to curves - this is why, for example, in track and field races, runners are not started at the same starting line if there are races that have curves. The same issue arises when vehicles are driving along a reference path. For example, in comparing similar paths plotted at same station points - in this case every 40 meters - each path's station may be very distant from its neighbors. The following MATLAB code illustrates this issue (Example 2 in the main script). Most of the challenges with vehicle guidance require paths, but unless a clear reference path is used, the station distance on one side of the road or even between adjacent lanes may be very different. This illustrates a need a COMMON station for the road or lane.

```MATLAB
%% Example 2: Show station discrepancies
% Plot station markers
 
fig_num = 2;
 
figure(fig_num);
clf;
hold on;
grid on;
grid minor;
 
Station_step = 40;
 
fcn_Path_plotTraversalsXY(data,fig_num);
xlabel('X [m]');
ylabel('Y [m]');
 
for i_traveral = 1:length(paths)
    traversal_stations = data.traversal{i_traveral}.Station;
    for i_station = Station_step:Station_step:traversal_stations(end)
        index = find(traversal_stations >= i_station,1);
        plot(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            'k.','Markersize',15);
        text(data.traversal{i_traveral}.X(index),...
            data.traversal{i_traveral}.Y(index),...
            sprintf('Station: %.2f',...
            data.traversal{i_traveral}.Station(index)));
    end
end

```

<pre align="center">
  <img src=".\Images\SimilarPathsMayHaveDifferentStationCoordinates.jpg" alt="definition of the coastline paradox" width="800" height="600">
  <figcaption>Similar paths can have different station coordinates</figcaption>
  <font size="-2">Photo from <a href="https://sketchplanations.com/the-coastline-paradox">Sketchplanations.com</a>
  </font>
</pre>

And swerving around a relative path can add significant length. The side-by-side swerving is why marathon runners ALWAYS end up running more than 26.2 miles (which often surprises and annoys first-time runners).

<pre align="center">
  <img src=".\Images\MarathonDistances.png" alt="marathon distances" width="400" height="600">
  <figcaption>Marathon runners often run significantly further than a marathon distance due to side-side swerving.</figcaption>
  <!--font size="-2">Photo from <a href="https://sketchplanations.com/the-coastline-paradox">Sketchplanations.com</a>
  </font -->
</pre>

Another issue with station coordinates is that the distance from a point to a path may be undefined, have multiple definitions, and/or is dependent on the history of the trajectory up to the point. This is illustrated in the figure below where a simple path connects points A, B, and C. The distance of point 1 to this path is unclear, as its position is outside the orthogonal projection of segment AB or segment BC. A solution may be to use projections of both segments, as illustrated with green dotted lines, but this projection can yield poor results if the point is near one projection but very far from another - an issue that can be solved by taking the maximum of the two distances, minimum, average, or weighted average (this issue is discussed in the functions below). Similarly, if a point is within the intersection of the orthogonal doimains of two segments, as is point 2 in the figure. Again, it is unclear how far this point is away from the path - whether to use segment AB or BC to answer this question. Again, this can be solved by taking the maximum of the two distances, minimum, average, or weighted average.

<pre align="center">
  <img src=".\Images\StationCoordinatesAreNotUnique.png" alt="station coordinates are not unique" width="520" height="370">
  <figcaption>Station Coordinates are not unique</figcaption>
  <font size="-2">Diagram created in Powerpoint</a>
  </font>
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Stations

In the codes that follow, there is often a need to distinguish between a single station and many stations. Therefore, in this library of codes, we treat station as a 1 x 1 scalar, and stations (plural) as a N x 1 vector. This way, if a function needs multiple stations, it can query the arguments for a vector 'stations' rather than a scalar 'station'.So a stations type variable consists of N `station' variables stacked together.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Path

In the codes that follow, a path is a set of [X Y] points as a N x 2 vector or array. These X Y rows denote the locations we are trying to follow which defines the path. For example:

```Matlab
    paths{1} = [
        8.1797    8.6006
        8.4101   15.8892
        10.0230   25.5102
        10.4839   39.7959];
```

A path type must have at least 2 rows so that a single 'path segment' can be defined. One connection between vertices is called a 'path segment', and a path must have at least ONE segment.

<pre align="center">
  <img src=".\Images\DefinitionOfPathSegment.jpg" alt="path segment" width="900" height="300">
  <figcaption>Definition of a path segment.</figcaption>
  <!--font size="-2">Photo from <a href="https://sketchplanations.com/the-coastline-paradox">Sketchplanations.com</a>
  </font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Paths

Some operations require at least 2 segments, and we refer to these variable types as 'paths'. The use of plural is to denote a type of variable that has multiple segments. An example operation that requires multiple segments is when we need to calculate the angles between segments.

We rarely use multiple paths as we have another variable type, traversals, discussed next that more clearly contains paths and properties of paths that are often used.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Traversal

Sometimes we want more information about a path than just XY, so we define something called a 'traversal'. It includes subfields that are routinely used, particularly the Station and Yaw. The codes below convert from Path to Traversals often, so a function exists that will fill in a Traversal given a Path.

```MATLAB
function traversal = fcn_Path_convertPathToTraversalStructure(path,varargin)
% fcn_Path_convertPathToTraversalStructure
% Takes a Path type and creates it to a traversal structure
%
% FORMAT: 
%
%       traversal = fcn_Path_convertPathToTraversalStructure(path,(fig_num))
%
% INPUTS:
%
%      path: an N x 2 vector of [X Y] positions, with N>=2
%
% OUTPUTS:
%
%      traversal: a sttructure containing the following fields
%            - X: an N x 1 vector that is a duplicate of the input X
%            - Y: an N x 1 vector that is a duplicate of the input Y
%            - Z: an N x 1 vector that is a zero array the same length as
%            the input X
%            - Diff: a [N x 2] array that is the change in X and Y
%            (front-padded with [0 0])
%            - Station: the XY distance as an N x 1 vector, representing
%            the distance traveled up to the current point (starting with 0
%            at the first point)
%            - Yaw: the calculated yaw angle (radians) of each path segment
%            (note: there are N-1 segments if there are N points, thus
%            there are N-1 Yaw points)
%            (thus - it is front-padded with NaN)
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.

```

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### Definition of Traversals

The last major type commonly used is 'traversals' which are collections of traversal variables within an array, e.g. many paths that are operated upon in code together. An example would be:

```MATLAB
Many_traversals.traversal{1} = some_traversal;
Many_ traversals.traversal{2} = some_other_traversal;

```

So, the variable 'Many_ traversals' is a 'traversals' type and contains an array of different traversals together. This variable allows us to create functions that operate on many traversals at the same time (such as plotting, averaging, variance calculations, etc.)

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Basic Getting Started Functions

#### fcn_Path_checkInputsToFunctions

 The function fcn_Path_checkInputsToFunctions is a general-use function whose only purpose is to confirm that variables that are passed to it meet the requirements of each type. It is often called repeatedly at the top of other codes to check that the inputs are the correct type. As of Dec. 2021, this function was replaced by those in the "Debug" repository, which are more powerful and general. However, this conversion to the new function has not yet been completed.

```MATLAB
%% Show how input arguments are checked, fcn_Path_checkInputsToFunctions
% TO-DO - in future versions, use debug tools as a utility and remove the checkinputs capability out of Path class
path_test = [4 1; 2 1];
fcn_Path_checkInputsToFunctions(path_test, 'path');

```

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_fillSamplePaths

The function fcn_Path_fillSamplePaths fills test data and includes a click-to-draw. The codes often require testing and creating paths, and this function is for dummy 'starter' paths to create sample test cases. It fills in (currently) three different paths in an array, which is easily converted into a 'traversals' type.

Note: a click-to-draw functionality also exists by setting a flag within the function itself. This allows interactive creation of paths in this function! An example call:

```MATLAB
% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
```

<pre align="center">
  <img src=".\Images\fcn_Path_fillSamplePaths.jpg" alt="fcn_Path_fillSamplePaths picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillSamplePaths creates test data sets for exercising PathClass functions.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Why do we use this example shown above? It has some challenging features that often break codes:

1. Each traversal has a VERY different number of points and station distances
2. The traversals cross each other repeatedly
3. They have different start/stop XY locations
4. Each traversal loops back on itself, which breaks Euclidean distance-based queries

Note that this code has, within it, the function, fcn_Path_fillPathViaUserInputs. This allows the user to click in a figure to define test paths. This is a function for the user to click on the figure to generate XY path. Points are collected and plotted until the user double clicks.

<pre align="center">
  <img src=".\Images\fcn_Path_fillPathViaUserInputs.jpg" alt="fcn_Path_fillPathViaUserInputs picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillPathViaUserInputs allows users to click into a figure to define test paths.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_calcDiffAnglesBetweenPathSegments

The function fcn_Path_calcDiffAnglesBetweenPathSegments calculates the change in angles between path segments. Specifically, it alculates the change in angles between path segments. It uses dot products and cross products to find the relative change in a path's angle. This avoids the use of arctan which does not work when segments are pointing toward -180 degrees.

```MATLAB
% Pick first path as reference_traversal structure
paths_array = fcn_Path_fillSamplePaths;
paths_to_check = paths_array{1};
fig_num = 11111; 
diff_angles = fcn_Path_calcDiffAnglesBetweenPathSegments(paths_to_check,fig_num);

```

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***
fcn_Path_calcYawFromPathSegments - Calculates yaw angle for each path segment.

#### fcn_Path_calcYawFromPathSegments

The function fcn_Path_calcYawFromPathSegments calculates yaw angle for each path segment. It does this by summing yaw changes to find the yaw angle along a path. An example call:

```MATLAB
%% Show how to calculate the yaw angles along a path, fcn_Path_calcYawFromPathSegments
 
% Basic call with one path
fig_num = 22222;
yaw_angles = fcn_Path_calcYawFromPathSegments(path_to_check,fig_num); %#ok<*NASGU>

```

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_convertPathToTraversalStructure

The function fcn_Path_convertPathToTraversalStructure converts paths into traversals. The traversal structure type is created because the row length of each traversal is different, and so we cannot save corresponding groups of paths as arrays (not easily, at least). So, we use structure types to save traversals to allow different lengths. Specificall, a traversal is a cell array of structures. For example, to plot many paths against each other, many paths are converted into a traversals type:

```MATLAB
% script_test_fcn_Path_plotTraversalsXY.m
% Tests fcn_Path_plotTraversalsXY
close all
clc
 
 
% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
 
 
% Convert paths into traversals
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_traveral});
    data.traversal{i_traveral} = traversal;
end
 
 
%% Call the plot command to show how it works. 
figure(11);
fcn_Path_plotTraversalsXY(data);

```

<pre align="center">
  <img src=".\Images\convertPathToTraversalStructure.jpg" alt="convertPathToTraversalStructure picture" width="400" height="300">
  <figcaption>The function convertPathToTraversalStructure creates a traversal structure type given a path input.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_plotTraversalsXY

The function fcn_Path_plotTraversalsXY plots the XY value of a traversals type. An example usage of the function is shown in the function above.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_plotTraversalsYaw

The function fcn_Path_plotTraversalsYaw plots the yaw values of the traversals, as we often want to plot the yaw within traversal arrays.

```MATLAB
% script_test_fcn_Path_plotTraversalsYaw.m
% tests fcn_Path_plotTraversalsYaw
 
clc
close all
 
% Fill in some dummy data
paths = fcn_Path_fillSamplePaths;
 
 
% Convert paths into traversals
for i_traveral = 1:length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(paths{i_traveral});
    data.traversal{i_traveral} = traversal;
end
 
 
%% Call the plot command to show how it works. First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.
figure(11);
fcn_Path_plotTraversalsYaw(data);
 
 
%% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 12;
fcn_Path_plotTraversalsYaw(data,fig_num);

```

<pre align="center">
  <img src=".\Images\fcn_Path_plotTraversalsYaw.jpg" alt="fcn_Path_plotTraversalsYaw picture" width="400" height="300">
  <figcaption>The function fcn_Path_plotTraversalsYaw plots the yaw angle of a traversals type.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Basic Path Operations

#### fcn_Path_findProjectionHitOntoPath

The function fcn_Path_findProjectionHitOntoPath finds the intersection between a sensor vector and a path. This is the most basic path operation in this library: finding the intersection of a path segment rooted at point p and extending to p+r, with a sensor (or path segment) rooted at q and extending to p+s. Diagramatically, it is shown as:

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_1.jpg" alt="intersection of 2 segments picture" width="200" height="200">
  <figcaption>The function fcn_Path_findProjectionHitOntoPath conducts the most basic and important analysis in this library: finding the intersection of two segments with each other.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This function is quite flexible, and includes many optional inputs and useful outputs:

```Matlab
function [distance,location,path_segment] = fcn_Path_findProjectionHitOntoPath(path,sensor_vector_start,sensor_vector_end,varargin)   
% fcn_Path_findProjectionHitOntoPath calculates hits between a sensor
% projection and a path, returning the distance and location of the hit.
%
% FORMAT: 
%
%      [distance,location,path_segment] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(fig_num))  
%
% INPUTS:
%
%      path: an N x 2 vector containing the X,Y points of the path to be
%      checked for intersections
%
%      sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
%      sensor's start location
%
%      sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
%      sensor's end location
%
%      (OPTIONAL INPUTS)
%      flag_search_type: an integer specifying the type of search.
%
%            0: return distance and location of first intersection only if
%            the given sensor_vector overlaps the path (this is the
%            default)
%
%            1: return distane and location of first intersection if any
%            projection of the sensor vector, in any direction, hits the
%            path (in other words, if there is any intersection). Note that
%            distance returned will be negative if the nearest intersection
%            is in the opposite direction of the given sensor vector.
%
%            2: returns distances and locations as M x 1 and M x 2 vectors
%            respectively, where the M rows represent ALL the detected
%            intersections. In cases where the sensor vector completely
%            overlaps a path segment, the start and end of overlap are
%            given.
%
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      distance: a 1 x 1 scalar representing the distance to the closest
%      intersection of the sensor with the path
%
%      location: a 1 x 2 vector of the X,Y location of intersection point
%
%      path_segment: the segment number of the path that was hit (1 is the
%      first segment, 2 is the second, etc)
```

Details on the derivation of this code can be found at: <a href="https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect">Stakoverflow: How do you detect where two line segments intersect?</a>. That discussion explains, in detail, many of the cases as well as provides the C-code implementation of line intersection checks in a manner that avoids division, and thus is very (very) fast. The details are important and so are repeated here.

First, the 2-d cross product is needed. While MATLAB includes a cross function, it is built for 3D implementation and thus requires extra steps. It is faster to define our own function that generates the 2-dimensional vector cross product $v \times w$ to be $v_x\cdot w_y - v_y \cdot w_x$. Note that the function is vectorized assuming column vectors. Suppose the path segment runs from $p$ to $p + r$ and the sensor starts from $q$ and extends to $q +s$. Then any point on the first line is representable as $p+ t\cdot r$ (for a scalar parameter $t$) and any point on the second line as $q + u\cdot s$ (for a scalar parameter $u$).

```Matlab
%% Calculate cross products
function result = crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end
```

The solution method used is to consider fractions of the path and of the sensor vectors. Specifically, the two line segments can be represented as vectors, which intersect if we can find $t$ and $u$ such that:

$p + t\cdot r = q + u\cdot s$

Pictorally, this looks like:

<pre align="center">
  <img src=".\Images\TwoSegmentsIntersect.jpg" alt="intersection of 2 segments picture 2" width="200" height="200">
  <figcaption>The intersection of 2 segments creates a vector equivalence at point of intersection</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

To solve, cross both sides with $s$, getting

$\left( p + t\cdot r \right) \times s = \left( q + u\cdot s \right) \times s$

And since  $s \times s = 0$ , this means:

$t \cdot \left( r \times s \right) = \left( q  - p \right) \times s$

And therefore, solving for $t$:

$\Large t = \Large \frac{\left( q  - p \right) \times s}{r \times s}$

We can solve for $u$ similarly:

$\left( p + t\cdot r \right) \times r = \left( q + u\cdot s \right) \times r$

And since  $r \times r = 0$ , this means:

$u \cdot \left( s \times r \right) = \left( p  - q \right) \times r$

And therefore, solving for $u$:

$\Large u = \Large \frac{\left( p  - q \right) \times r}{s \times r}$

Note that the calculation of u can reuse the $r \times s$ calculation from the $t$ equation, noting that $s \times r = - r \times s$. The result:

$\Large u = \Large \frac{\left( q  - p \right) \times r}{r \times s}$

There are several cases to consider, and following the notation on Stackoverflow:

1. If $r \times s = 0$ and $(q - p) \times r = 0$, then the two lines are collinear (more on this case shortly).
2. If $r \times s = 0$ and $(q - p) \times r \neq  0$, then the two lines are parallel and non-intersecting.
3. If $r \times s \neq 0$ and $0 \leq t \leq 1$ and $0 \leq u \leq 1$, the two line segments meet at the point $p + t r = q + u s$.
4. Otherwise, the two line segments are not parallel but do not intersect.

As noted therein: this method is the 2-dimensional specialization of the 3D line intersection algorithm from the article "Intersection of two lines in three-space" by Ronald Goldman, published in Graphics Gems, page 304. In three dimensions, the usual case is that the lines are skew (neither parallel nor intersecting) in which case the method gives the points of closest approach of the two lines.

Looking at the parallel case in more detail: $r \times s = 0$, we can check if they are colinear or disjoint. This is done by checking the amount of overlap. Specifically, we express the endpoints of the second segment ($q$ and $q +s$) in terms of the equation of the first line segment $(p + tr)$:

$t_0 = (q -p) \cdot r / (r \cdot r)$

$t_1 = (q + s - p) \cdot r / (r \cdot r) = t_0 + s\cdot r/ (r \cdot r)$

If the interval between $t_0$ and $t_1$ intersects in the interval [0, 1], then the line segments are collinear and overlapping (Case 1); otherwise, they are collinear and disjoint (Case 2). Note, as mentioned in that article, that if $s$ and $r$ point in opposite directions, then $s \cdot r < 0$ and so the interval to be checked is $[t_1, t_0]$ rather than $[t_0, t_1]$.

In programming the algorithm, one must be careful of intersections at the end-points of the path, as shown below.

<pre align="center">
  <img src=".\Images\IntersectionsAtEndOfSegments.jpg" alt="Intersections misisng End Of Segments" width="400" height="300">
  <figcaption>Intersections at the end of segments can create conditions where no intersections are detected, if not handled carefully</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

If intersections occur at the end-points of a path segment, then results can be ambiguous because $t = 0$ or $t = 1$.  In these cases, can have weird situations if we omit this, e.g. the following:

$((0 < t).*(1 > t).*(0<u).*(1>u))$

The result of this condition is that NO intersection would be detected!

If we correct the t-range to allow equality, then the intersection is found.

$ ((0 \leq t).*(1 \geq t).*(0<u).*(1>u))$

<pre align="center">
  <img src=".\Images\IntersectionsAtEndOfSegmentsGood.jpg" alt="Intersections At End Of Segments" width="400" height="300">
  <figcaption>Fixing the t range check fixes errors of intersections on the path.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

However, this form with the $u$ range not including 0 or 1 can miss paths that barely hit the sensor at start or end. For example:

<pre align="center">
  <img src=".\Images\IntersectionsMissingOriginOrEndpoints.jpg" alt="Intersections At End Of Segments" width="800" height="300">
  <figcaption>Intersections at the end or start of the sensor can have no intersections detected, if not handled carefully</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This can be fixed by correcting the inequality:

$ ((0 \leq t).*(1 \geq t).*(0 \leq u).*(1 \geq u))$

<pre align="center">
  <img src=".\Images\IntersectionsHittingOriginOrEndpoints.jpg" alt="Intersections At End Of Segments" width="800" height="300">
  <figcaption>Intersections at the end or start of the sensor can have no intersections detected, if not handled carefully</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

With this correct formulation, we can present several simple examples (See: script_test_fcn_Path_findProjectionHitOntoPath.m):

Example 1: A simple intersection

```Matlab
%% Simple test 1 - a simple intersection
fprintf(1,'Simple intersection result: \n');
path = [0 10; 10 10];
sensor_vector = [2 1; 5 15];
fig_debugging = 2343;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(path,sensor_vector,fig_debugging);
print_results(distance,location);
```

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex1.jpg" alt="Example 1: A simple intersection" width="400" height="300">
  <figcaption>Example 1: A simple intersection</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Example 2: not intersecting

```Matlab
%% Simple test 2 - no intersections
fprintf(1,'No intersection result: \n');
path = [-4 10; 2 10];
sensor_vector = [0 0; 5 12];
fig_debugging = 2343;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(path,sensor_vector,fig_debugging);
print_results(distance,location);
```

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex2.jpg" alt="Example 2: no intersection" width="400" height="300">
  <figcaption>Example 2: No intersection</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Example 3: multiple line segments, only the first detection is recorded

```Matlab
%% Simple test 3 - multiple intersections
fprintf(1,'Multiple intersections result: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector = [0 0; 5 12];
fig_debugging = 2343;
[distance,location] = fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector,fig_debugging);
print_results(distance,location);

```

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex3.jpg" alt="Example 3: multiple line segments" width="400" height="300">
  <figcaption>Example 3: Multiple line segments</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

In some applications, the length of the sensor vector is irrelevant. One simply wants to know any location, in any direction, where the first intersection would occur. The flag_search_type variable can be set to 1 (instead of 0, default) to specify that any intersection will work. Here is an example with a sensor vector of length 2. If an unlimited search is used (e.g. with a flag of 1), a hit is detected with either positive or negative distance.

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex4.jpg" alt="Example 4: the search type" width="800" height="500">
  <figcaption>Changing the search type to 1 gives the first detections at any distance, positive or negative</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

In some applications, we want to know all the intersection points, not just the first one. If the search type is set to 2, then this can be found:

```Matlab
%% Advanced test 2 - multiple intersections
fprintf(1,'Single intersections reporting only first result: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 23487;
flag_search_type = 0;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
 
fprintf(1,'Multiple intersections reporting all results: \n');
path = [0 10; 10 10; 0 6; 10 6; 0 2];
sensor_vector_start = [0 0]; 
sensor_vector_end   = [5 12];
fig_debugging = 23488;
flag_search_type = 2;
[distance,location] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_results(distance,location);
```

Which generates the following console result:

``` c
Single intersections reporting only first result: 
Distance   Location X   Location Y 
2.600      1.000    2.400
Multiple intersections reporting all results: 
Distance   Location X   Location Y 
10.833     4.167    10.000
7.800      3.000    7.200
6.500      2.500    6.000
2.600      1.000    2.400
```

And the following figure

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex5.jpg" alt="Example 5: the search type" width="800" height="500">
  <figcaption>Changing the search type to 2 gives all detections at any distance, positive or negative</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The algorithm is quite robust even to complex overlapping conditions such as multiple completely overlapping segments mixed with partial intersections.

```Matlab
%% Advanced Multihit Overlapping  test 32 - super overlapping colinear 1 BACKWARDS
fprintf(1,'Super overlapping colinear BACKWARDS result: \n');
path = [0 10; 10 10; 12 8; 14 10; 15 10];  
sensor_vector_end = [-3 10]; 
sensor_vector_start   = [15 10];
fig_debugging = 2343;
flag_search_type = 2;
[distance,location,path_segments] = ...
    fcn_Path_findProjectionHitOntoPath(...
    path,sensor_vector_start,sensor_vector_end,...
    flag_search_type,fig_debugging);
print_more_results(distance,location,path_segments);
```

which produces the following console result:

``` c
GIVES:
Super overlapping colinear BACKWARDS result: 
Distance     Location X      Location Y      PathSegment 
5.000        10.000              10.000          1
5.000        10.000              10.000          2
1.000        14.000              10.000          3
0.000        15.000              10.000          4
15.000       0.000               10.000          1
1.000        14.000              10.000          4
```

And the following figure

<pre align="center">
  <img src=".\Images\fcn_Path_findProjectionHitOntoPath_Ex6.jpg" alt="Example 4: the search type" width="400" height="300">
  <figcaption>The algorithm is quite robust to very difficult situations</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_findIntersectionsBetweenTraversals

The function fcn_Path_findIntersectionsBetweenTraversals finds where two traversals intersect. It is a traversal-type implementation of the prior function: fcn_Path_findProjectionHitOntoPath.

<pre align="center">
  <img src=".\Images\fcn_Path_findIntersectionsBetweenTraversals.jpg" alt="fcn_Path_findIntersectionsBetweenTraversals picture" width="800" height="300">
  <figcaption>The function fcn_Path_findIntersectionsBetweenTraversals finds where two traversals intersect.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

See script_test_fcn_Path_findIntersectionsBetweenTraversals for this demo.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_snapPointOntoNearestPath

The function fcn_Path_snapPointOntoNearestPath snaps a point onto the nearest path segment. This is one of the basic path operations. An illustration of the steps is shown below.

<pre align="center">
  <img src=".\Images\fcn_Path_snapPointOntoNearestPath_steps.jpg" alt="fcn_Path_snapPointOntoNearestPath picture" width="800" height="300">
  <figcaption>The function fcn_Path_snapPointOntoNearestPath snaps a point onto the nearest path segment.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The behavior of the snap depends on regions, whether or not the query point lands in the segment before the closest snap point, or the segment after the closest snap point. The code is separated by each of these 4 regions.

<pre align="center">
  <img src=".\Images\fcn_Path_snapPointOntoNearestPath_regions.jpg" alt="fcn_Path_snapPointOntoNearestPath regions" width="400" height="500">
  <figcaption>The function fcn_Path_snapPointOntoNearestPath snaps a point onto the nearest path segment.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The regions are defined to ensure continuity of distance calculations across the borders of regions. However, the interior of the region has a boundary where the solution must choose, and therefore jump, between either the "before" segment or "after" segment.

<pre align="center">
  <img src=".\Images\fcn_Path_snapPointOntoNearestPath_3D_distance.jpg" alt="fcn_Path_snapPointOntoNearestPath regions" width="400" height="300">
  <figcaption>The snap boundaries are mathematically defined to avoid jump discontinuities in distance.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Additional examples follow. See script_test_fcn_Path_snapPointOntoNearestPath.m for a comprehensive test suite.

Example 1: a typical query

```Matlab
%% BASIC example 1.2 - works
point = [1.4 1.3]; % Define the query point
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path
fignum = 112; % Define the figure number

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);
```

Produces the following results

``` c
Figure: 112
   Closest point is: 1.14 0.78 
   Matched to the path segment given by indices 3 and 4, 
   S-coordinate is: 1.61, 
   percent_along_length is: 0.40
```

<pre align="center">
  <img src=".\Images\fcn_Path_snapPointOntoNearestPath_Ex1.jpg" alt="fcn_Path_snapPointOntoNearestPath picture" width="400" height="300">
  <figcaption>The function fcn_Path_snapPointOntoNearestPath snaps a point onto the nearest path segment.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Example 2: a query outside both segments

```Matlab
%% BASIC example 1.4 - works, but shows that it is on neither segement
point = [0.9 1.4]; 
pathXY = [0 0; 0.5 0.2; 0.9 0.9; 1.5 0.6; 3 0]; % Define an XY path

fignum = 114;

% Snap the point onto the path
[closest_path_point,s_coordinate,...
    first_path_point_index,second_path_point_index,...
    percent_along_length] = ...
    fcn_Path_snapPointOntoNearestPath(point, pathXY,fignum);

% Print results to the workspace
fprintf(1,'Figure: %d\n',fignum);
fprintf(1,'\t\t Closest point is: %.2f %.2f \n',...
    closest_path_point(1,1),closest_path_point(1,2));
fprintf(1,'\t\t Matched to the path segment given by indices %d and %d, \n',...
    first_path_point_index,second_path_point_index);
fprintf(1,'\t\t S-coordinate is: %.2f, \n',...
        s_coordinate);
fprintf(1,'\t\t percent_along_length is: %.2f\n',percent_along_length);
```

Produces the following results

``` c
Figure: 114
   Closest point is: 0.90 0.90 
   Matched to the path segment given by indices 3 and 3, 
   S-coordinate is: 1.34, 
   percent_along_length is: 0.00
```

And the following figure:

<pre align="center">
  <img src=".\Images\fcn_Path_snapPointOntoNearestPath_Ex2.jpg" alt="fcn_Path_snapPointOntoNearestPath picture 2" width="400" height="300">
  <figcaption>The function fcn_Path_snapPointOntoNearestPath will use the nearest point if it is outside the orthogonal projection of a line segment.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_snapPointOntoNearestTraversal

The function fcn_Path_snapPointOntoNearestTraversal snaps point to nearest traversal. It is basically a function that calls fcn_Path_snapPointOntoNearestPath with traversal information converted to path format.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Functions that Trim Paths

#### fcn_Path_findTraversalStationSegment

The function fcn_Path_findTraversalStationSegment crops traversal by given station interval. It was written because, when snapping points onto paths or traversals, if we include too much of a path, particularly one turned back toward itself, we can get incorrect snap points.

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_IncorrectSnapPoints.jpg" alt="fcn_Path_findTraversalStationSegment with IncorrectSnapPoints" width="400" height="250">
  <figcaption>Fig.4 - The function fcn_AlignCoords_fitRotationKabsch performs regression fitting to find the best-fit rotation and translation that matches one set of points to another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

However, if we limit the traversal area being studied to only a region of interest, we may snap to the correct point.

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_CorrectSnapPoints.jpg" alt="fcn_Path_findTraversalStationSegment with Correct Snap Points" width="400" height="250">
  <figcaption>Fig.4 - The function fcn_AlignCoords_fitRotationKabsch performs regression fitting to find the best-fit rotation and translation that matches one set of points to another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The function fcn_Path_findTraversalStationSegment takes in the s-coordinates of the start and end area of interest, and returns the portion of the path that includes these station coordinates. As well, it can set a flag if the query coordinates are outside the path either on the start side or end side of the given traversal. Here's the setup code:

```Matlab
% script_test_fcn_Path_findTraversalStationSegment.m
% This is a script to exercise the function: 
% fcn_Path_findTraversalStationSegment.m
% This function was written on 2020_11_16 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
 
%      [traversal_trimmed,flag_outside_start, flag_outside_end] = ...
%      fcn_Path_findTraversalStationSegment(...
%      long_traversal, s_coord_start,s_coord_end, 
%      (fig_num))
 
% Revision history:
%     2021_01_09
%     -- updated name and types to take traversal inputs
%     -- added input checking
%     -- added flag_do_plots    
 
close all;
clear data;
 
% Fill in sample paths (as a starter)
paths_array = fcn_Path_fillSamplePaths;
 
% Convert paths to traversal structures
for i_Path = 1:1  % length(paths)
    traversal = fcn_Path_convertPathToTraversalStructure(...
        paths_array{i_Path});
    data.traversal{i_Path} = traversal;
end
 
% Plot the results?
if 1==1
    fig_num = 12;
    fcn_Path_plotTraversalsYaw(data,fig_num);
 
    fig_num = 13;
    fcn_Path_plotTraversalsXY(data,fig_num);
end
```

And here's the function call:

```Matlab
%% BASIC example 1
s_coord_start = 10;
s_coord_end   = 100;
fignum = 111;
[traversal_segment1,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Normal query - Test case #1');
```

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_Ex1.jpg" alt="fcn_Path_findTraversalStationSegment Example 1" width="400" height="300">
  <figcaption>The function fcn_Path_findTraversalStationSegment returns the portion of the path containing the given start and end station coordinates.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The function is robust to degenerate queries. If a degenerate query is given where there is only a single station, then it automatically finds the segment of the traversal portion containing that station.

```Matlab
%% BASIC example 4 - degenerate within
s_coord_start = 70;
s_coord_end   = 70;
fignum = 444;
[traversal_segment4,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query with exact same point');

```

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_Ex4.jpg" alt="fcn_Path_findTraversalStationSegment Example 4" width="400" height="300">
  <figcaption>If given a single station, the function fcn_Path_findTraversalStationSegment returns the path segment containing those station coordinates.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

And if the ENTIRE traversal is within the station limits, it simply returns the entire traversal back:

```Matlab
%% BASIC example 5 - entire traversal within
s_coord_start = -5;
s_coord_end   = 7000;
fignum = 555;
[traversal_segment5,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where entire traversal within start and end of s-limits');

```

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_Ex5.jpg" alt="fcn_Path_findTraversalStationSegment Example 5" width="400" height="300">
  <figcaption>If given a station range beyond that of the entire traversal, then the entire traversal is returned.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

If the station end query is past the end of the traversal, but start is within, it returns only the cut-out portion:

```Matlab
%% BASIC example 6 - end of traversal out
s_coord_start = 200;
s_coord_end   = 7000;
fignum = 666;
[traversal_segment6,flag_outside_start, flag_outside_end] = ...
    fcn_Path_findTraversalStationSegment(...
    traversal, s_coord_start,s_coord_end, fignum);
fprintf(1,...
    'Figure: %d, flag_outside_start is: %d, flag_outside_end is: %d \n',...
    fignum, flag_outside_start,flag_outside_end);
title('Query where end point outside of s-limits of traversal');

```

<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalStationSegment_Ex6.jpg" alt="fcn_Path_findTraversalStationSegment Example 6" width="400" height="300">
  <figcaption>If the station end query is past the end of the traversal, but start is within, it returns only the cut-out portion:
.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Similarly, if the station start query is before the start of the traversal, but end is within, it returns only the cut-out portion; this is not shown but is the 7th basic example in the test script. If both the start and end of query are BEFORE the start, then it returns only the first segment; this is the 9th basic example. And similarly, if both the start and end of query are AFTER the end, then it returns only the last segment; this is the 10th basic example. This function will still give the correct result even if the start and end points are out of order, but will generate a warning as it auto-corrects this issue.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_removePinchPointInTraversal

The function fcn_Path_removePinchPointInTraversal eliminates self-crossings in traversals.  Pinch points are where a path crosses back onto itself, creating a loop. The path without a pinch point is the minimum traversal moving only along the original path traveling from start to end. This traversal will avoid any loop or self-crossing area. The resulting path is one without "pinch points".

<pre align="center">
  <img src=".\Images\fcn_Path_removePinchPointInTraversal_illustration.jpg" alt="fcn_Path_removePinchPointInTraversal illustration" width="400" height="300">
  <figcaption>The function fcn_Path_removePinchPointInTraversal eliminates self-crossings in traversals.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This function is particularly useful to correct offset paths, which twist back onto themselves when created by an offset projection process. This function is also very useful to clean up weird polytopes that are self-crossing when traversing around an edge.

<pre align="center">
  <img src=".\Images\fcn_Path_removePinchPointInTraversal_pathOffsetRepair.jpg" alt="fcn_Path_removePinchPointInTraversal illustration" width="800" height="300">
  <figcaption>The function fcn_Path_removePinchPointInTraversal is particularly useful to correct offset paths, which twist back onto themselves .</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Functions That Project Paths

#### fcn_Path_findOrthogonalTraversalVectorsAtStations

The function fcn_Path_findOrthogonalTraversalVectorsAtStations calculates orthogonal vectors to a traersal at given stations. A common operation is to seek vectors that are ortho-normal to a path at given station points. This is used in calculating nearby paths, in projecting paths sideways, and similar operations. This function does this calculation of orthogonal vectors.

```Matlab
function [unit_normal_vector_start, unit_normal_vector_end] = ...
    fcn_Path_findOrthogonalTraversalVectorsAtStations(station_queries,central_traversal, varargin)
 
% fcn_Path_findOrthogonalTraversalVectorsAtStations
% Given a central traversal and a set of stations along that traversal,
% finds the unit normal vector on the central traveral at each station
% point.
% 
% FORMAT: 
%
%      [unit_normal_vector_start, unit_normal_vector_end] = ...
%        fcn_Path_findOrthogonalTraversalVectorsAtStations(...
%        station_queries,central_traversal,...
%        (flag_rounding_type),(fig_num));
%
% INPUTS:
%
%      station_queries: an N x 1 vector containing the station on the
%      central traversal where the projections should take place
%
%      central_traversal: a traversal structure that specifies the path
%      where projections to other paths are taking place.
%
%      (OPTIONAL INPUTS)
%      flag_rounding_type: a flag to indicate which type of projection is
%      used, especially when stations are located at the end-points of
%      segments within the nearby_traversal. When stations are at the
%      end-points of segments, the normal vector is undefined as it depends
%      on whether to use the prior or subsequent segment, or some
%      combination of these.
%
%      Note that the very first point always uses projections from the
%      following segement, and the very last point always uses the prior.
%      Otherwise, the flag determines behaviors for endpoints of internal
%      segments. The options include:
%
%          flag_rounding_type = 1;  % This is the default, and indicates
%          that the orthogonal projection of an endpoint is created by the
%          PRIOR segment leading up to each station query point.
% 
%          flag_rounding_type = 2;  % This indicates that the orthogonal
%          projection of an endpoint is created by the FOLLOWING segment
%          after each station query point.
% 
%          flag_rounding_type = 3;  % This indicates that the orthogonal
%          projection, ONLY if the station query falls at the joining point
%          between two segments (e.g. is on the "joint"), then the
%          projection is created by averaging the vector projections
%          created from the PRIOR segment and FOLLOWING segment.
% 
%          flag_rounding_type = 4;  % This indicates that the orthogonal
%          projections along segments should be calculated at the midpoints
%          of each segment, and then for each station qeuary, the vector
%          projections are interpolated from their prior and subsequent
%          vectors.
%
%      fig_num: a figure number to plot results. Turns debugging on.
%
% OUTPUTS:
%
%      unit_normal_vector start: a Nx2 vector containing [X1 Y1]
%      coordinates as columns, where the [X1 Y1] represents the location of
%      the start point of the vector, on the path.
%
%      unit_normal_vector_end: a Nx2 vector containing the [X2 Y2] location
%      of the end point of the unit vector. On both outputs, there are N
%      rows, one row for each station.
%
% DEPENDENCIES:
%
%      fcn_Path_checkInputsToFunctions
%
% EXAMPLES:
%      
% See the script: script_test_fcn_Path_findOrthogonalTraversalVectorsAtStations
% for a full test suite.
%
```

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalTraversalVectorsAtStations_illustration.jpg" alt="fcn_Path_findOrthogonalTraversalVectorsAtStations picture" width="400" height="300">
  <figcaption>The function fcn_Path_findOrthogonalTraversalVectorsAtStations calculates "normal"" vectors to a traversal at given stations.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Here are some examples:

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalTraversalVectorsAtStations_Ex1.jpg" alt="fcn_Path_findOrthogonalTraversalVectorsAtStations example 1" width="800" height="300">
  <figcaption>Examples of the function fcn_Path_findOrthogonalTraversalVectorsAtStations.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Why is "normal" in quotes? Because the definition of normal is unclear! Particularly at an internal vertex, it can be very unclear what to use as a normal vector.

Here are 4 options:

1. use previous segment so projection vector to define orthogonality,
2. use following segment,
3. average both previous and following, using average only at vertex, and
4. average both everywhere (which is smooth, but almost no vectors are actually orthogonal!)

The code allows all four to be selected via a flag option, and the differences between each of these options is illustrated below:

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalTraversalVectorsAtStations_Ex2.jpg" alt="fcn_Path_findOrthogonalTraversalVectorsAtStations example 1" width="900" height="200">
  <figcaption>Examples of the function fcn_Path_findOrthogonalTraversalVectorsAtStations.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_findOrthogonalHitFromTraversalToTraversal

The function fcn_Path_findOrthogonalHitFromTraversalToTraversal finds which traversals are hit at ortho projections from one traversal to another. Basically, it finds the points on traversals nearby a given traversal that are closest to given station points, assuming orthogonal projections or variants set by a flag consistent with the flag for fcn_Path_findOrthogonalTraversalVectorsAtStations. An illustration of this process is shown below; This is from the averaging example in the script: script_test_fcn_Path_findOrthogonalHitFromTraversalToTraversal.m

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalHitFromTraversalToTraversal_Ex1.jpg" alt="fcn_Path_findOrthogonalHitFromTraversalToTraversal Example 1 picture" width="900" height="200">
  <figcaption>The function fcn_Path_findOrthogonalHitFromTraversalToTraversal finds which traversals are hit at ortho projections from one traversal to another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Note that the result will give a positive or negative distance value, depending if it is to the left (positive cross-product) or right (negative) direction relative to the traversal.

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalHitFromTraversalToTraversal_Ex2.jpg" alt="fcn_Path_findOrthogonalHitFromTraversalToTraversal Example 2 picture" width="800" height="300">
  <figcaption>The function fcn_Path_findOrthogonalHitFromTraversalToTraversal can report both positive and negative intersection distances, with the sign determined by the cross-product convention.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

One can limit the search radius for intersections, in this case to 1.5 meters. Note that the results are different depending on the projection type.

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalHitFromTraversalToTraversal_Ex3.jpg" alt="fcn_Path_findOrthogonalHitFromTraversalToTraversal Example 3 picture" width="900" height="300">
  <figcaption>The function fcn_Path_findOrthogonalHitFromTraversalToTraversal has a flag to limit the search distance away from the central traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The search direction can also be negative. Here are the results for a search radius of 1.5 meters, showing that negative results are captured.

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthogonalHitFromTraversalToTraversal_Ex5.jpg" alt="fcn_Path_findOrthogonalHitFromTraversalToTraversal Example 4 picture" width="900" height="200">
  <figcaption>The function fcn_Path_findOrthogonalHitFromTraversalToTraversal is particularly useful to find nearby traversal points.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_findOrthoScatterFromTraversalToTraversals

The function fcn_Path_findOrthoScatterFromTraversalToTraversals finds closest points on many traversals to a given central traversal. Like the prior function, it gives positive values for positive (in the cross-product) sensor directions, negative for negative sensing.

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthoScatterFromTraversalToTraversals_Ex1.jpg" alt="fcn_Path_findOrthoScatterFromTraversalToTraversals picture" width="800" height="300">
  <figcaption>The function fcn_Path_findOrthoScatterFromTraversalToTraversals finds closest points on many traversals to a given central traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This is useful for finding all the nearby traversal points, including defining a search distance that ignores nearby traversals

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthoScatterFromTraversalToTraversals_Ex2.jpg" alt="fcn_Path_findOrthoScatterFromTraversalToTraversals picture" width="400" height="300">
  <figcaption>The function fcn_Path_findOrthoScatterFromTraversalToTraversals is useful to find common station points for repeated measurements along a reference traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

With dense projections, one can analyze the histogram of projection distances:

<pre align="center">
  <img src=".\Images\fcn_Path_findOrthoScatterFromTraversalToTraversals_Ex3.jpg" alt="fcn_Path_findOrthoScatterFromTraversalToTraversals picture" width="400" height="300">
  <figcaption>Dense projections allow error analysis of paths along a traversal, to find systematic bias, skew, etc.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_fillOffsetTraversalsAboutTraversal

The function fcn_Path_fillOffsetTraversalsAboutTraversal fills in an array of traversals about a reference traversal at user-defined offset distances.

<pre align="center">
  <img src=".\Images\fcn_Path_fillOffsetTraversalsAboutTraversal_Ex1.jpg" alt="fcn_Path_fillOffsetTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillOffsetTraversalsAboutTraversal fills in an array of traversals about a reference traversal at user-defined offset distances..</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This is particularly useful for finding the intersection points of road boundaries

<pre align="center">
  <img src=".\Images\fcn_Path_fillOffsetTraversalsAboutTraversal_Ex2.jpg" alt="fcn_Path_fillOffsetTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>This is particularly useful for finding the intersection points of road boundaries.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_convertTraversalXYtoSy

The function fcn_Path_convertTraversalXYtoSy calculates the SY (e.g ST or "station") coordinates for a given traversal, given a reference traversal and another query traversal.

<pre align="center">
  <img src=".\Images\fcn_Path_convertTraversalXYtoSy_Ex1.png" alt="fcn_Path_convertTraversalXYtoSy picture" width="400" height="300">
  <figcaption>The function fcn_Path_convertTraversalXYtoSy calculates the SY (e.g ST or "station") coordinates for a given traversal, given a reference traversal and another query traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>
The results can be very sensitive to station decimation, so this is a user-defined variable.
<pre align="center">
  <img src=".\Images\fcn_Path_convertTraversalXYtoSy_Ex2.png" alt="fcn_Path_convertTraversalXYtoSy picture" width="400" height="300">
  <figcaption>The results can be very sensitive to station decimation, so this is a user-defined variable.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_fillRandomTraversalsAboutTraversal

The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.

<pre align="center">
  <img src=".\Images\fcn_Path_fillRandomTraversalsAboutTraversal_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

One can choose different smoothness factors or the number of trajectories produced.

<pre align="center">
  <img src=".\Images\fcn_Path_fillRandomTraversalsAboutTraversal_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Choosing different smoothness factors or the number of trajectories produced.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

One can change the standard deviations

<pre align="center">
  <img src=".\Images\fcn_Path_fillRandomTraversalsAboutTraversal_Ex3.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Changing the standard deviations.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

One can also change the number of points to use to generate each random traversal

<pre align="center">
  <img src=".\Images\fcn_Path_fillRandomTraversalsAboutTraversal_Ex4.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Changing the number of points to use to generate each random traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Statistical Analysis of Paths

Basic path statistics operations include

Plotting a bound around a traversal
Calculating the default variance of a traversal
Plotting variance bands about a traversal

These basic functions are described below.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_plotTraversalXYWithUpperLowerBands

In the statistics functions, one often needs to plot bounds around a central traversal. This is done within the function fcn_Path_plotTraversalXYWithUpperLowerBands

<pre align="center">
  <img src=".\Images\fcn_Path_plotTraversalXYWithUpperLowerBands.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Tthe function fcn_Path_plotTraversalXYWithUpperLowerBands plot bounds around a central traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_calcSingleTraversalStandardDeviation

The variance of a single traversal is a measure of how much it bends at each segment, as a distance offset. 

This is useful to estimate the error in decimation of a path, as this variance becomes smaller with finer decimations.

<pre align="center">
  <img src=".\Images\fcn_Path_calcSingleTraversalStandardDeviation.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>This is useful to estimate the error in decimation of a path, as this variance becomes smaller with finer decimations.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_plotTraversalXYWithVarianceBands

The orthogonal projections can also be used to plot the standard deviation about a traversal

<pre align="center">
  <img src=".\Images\fcn_Path_plotTraversalXYWithVarianceBands_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The orthogonal projections can also be used to plot the standard deviation about a traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Note that this function, fcn_Path_plotPathXYWithVarianceBands inherets the band color from plotting, so that multiple traversals can be put on the same figure.

<pre align="center">
  <img src=".\Images\fcn_Path_plotTraversalXYWithVarianceBands_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>This function, fcn_Path_plotPathXYWithVarianceBands inherets the band color from plotting, so that multiple traversals can be put on the same figure.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Path Averaging Methods

Basic path averaging operations include

Averaging paths based on station
Averaging paths based on proximity
Averaging paths based on orthogonal projection
Find the traversal with the most data

These basic functions are described below.

One of the more useful ways to create permanent paths is to average ones followed earlier. These naturally occurring paths are called “desire lines”.

<pre align="center">
  <img src=".\Images\path_averaging_methods_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Path Averaging Methods</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The averaging process seems obvious, but it depends on metrics of distance

There are many!

1. Averaging points with the same station for different paths
2. Averaging the closest points from other paths to given points
3. Averaging points found by orthogonal projection from a path

Each of the averaging methods depends on how to define closeness of path points to each other

<pre align="center">
  <img src=".\Images\path_averaging_methods_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>TEach of the averaging methods depends on how to define closeness of path points to each other</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This discrepancy in station distances is why lanes on track and field have staggered starts. The inside line would be shorter than the outside lane if everyone started at the same line.

<pre align="center">
  <img src=".\Images\path_averaging_methods_Ex3.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The inside line would be shorter than the outside lane if everyone started at the same line.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Measuring nearness of a point to a path: three methods

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_findAverageTraversalViaStationAlignment

The “same station” method is just to average paths by station. This is done easily by just interpolation of all paths to same stations, then averaging.

It fails, however, when paths have different path lengths. As well, it tends to average points further and further apart from each other when moving along the path. 

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaStationAlignment_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The “same station” method is just to average paths by station. This is done easily by just interpolation of all paths to same stations, then averaging.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The “closest point” method finds the closest point on a nearby path to the central path. It then finds the associated line segment on the nearby path, and then projects FROM the line segment back to the central path.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaStationAlignment_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>TThe “closest point” method finds the closest point on a nearby path to the central path. It then finds the associated line segment on the nearby path, and then projects FROM the line segment back to the central path.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### fcn_Path_findAverageTraversalViaClosestPoint

The function, fcn_Path_findAverageTraversalViaClosestPoint implements closest-point averaging

The problem with closest point averaging is that the directions of the contributions can be unclear, e.g. the nearby path may be ahead or behind the station points on the central path.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function, fcn_Path_findAverageTraversalViaClosestPoint implements closest-point averaging.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Generating average trajectories by averaging orthogonal projections

Orthogonal projection takes a “central” trajectory, and then projects orthogonally from that trajectory at given stations to find where it hits nearby trajectories.

```MATLAB
%% BASIC example 1 - parallel lines, query is in middle area
stations = 1; % Define the station
 
% Create a dummy central path and convert it to a traversal
central_path = [0 0; 2 0];  
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
 
% Define a "nearby" path and convert it to a traversal
nearby_path = [0 4; 2 4];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
 
% Calculate the closest point and distance on the nearby path
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal);
 
print_results(stations,closest_path_point,distances);
 
```
<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>TOrthogonal projection takes a “central” trajectory, and then projects orthogonally from that trajectory at given stations to find where it hits nearby trajectories.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The approach gives the intersection point for arbitrary segments nearby the central trajectory

```MATLAB
%% BASIC example 2 - angled line segment adjacent to endpoint query
stations = 1;
central_path = [0 0; 2 0];
central_traversal = fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
 
nearby_path = [0 4; 2 7];
nearby_traversal =  fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
 
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal);
 
print_results(stations,closest_path_point,distances);

```
<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex3.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The approach gives the intersection point for arbitrary segments nearby the central trajectory.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

We define a “miss” to even include “grazing” from one trajectory to another

```MATLAB
%% BASIC example 3 - angled line segment adjacent to endpoint query but near-miss
stations = 10;
central_path = [0 0; 10 0];
central_traversal = ...
    fcn_Path_convertXYtoTraversalStructure(central_path(:,1),central_path(:,2));
nearby_path = [0 4; 10 7];
nearby_traversal = ...
    fcn_Path_convertXYtoTraversalStructure(nearby_path(:,1),nearby_path(:,2));
 
[closest_path_point,distances] = ...
    fcn_Path_FindOrthogonalHitFromPathToPath(stations,central_traversal,nearby_traversal);
```
<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex4.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Define a “miss” to even include “grazing” from one trajectory to another</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

A challenge with query stations is that the orthogonal projection is unclear. At the start and end, we use the segment ahead and behind these points.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex5.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>A challenge with query stations is that the orthogonal projection is unclear. At the start and end, we use the segment ahead and behind these points.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

As noted before in the orthogonal projection, the normal vectors are unclear at vertex points, and the path intersections can change depending on which type of averaging is used.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex6.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>As noted before in the orthogonal projection, the normal vectors are unclear at vertex points, and the path intersections can change depending on which type of averaging is used.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The same challenge arises for queries to the involuted area

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex7.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The same challenge arises for queries to the involuted area</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

One result of this process is that some areas of a nearby path may receive very poor sampling. Generally, portions of nearby paths that are parallel are sampled well, but perpendicular segments are not sampled at all.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex8.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>Generally, portions of nearby paths that are parallel are sampled well, but perpendicular segments are not sampled at all.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

For orthogonal projection, the option of averaging at the vertex appears to give the best results (e.g. flags 3 and 4 in orthogonal projections)

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex9.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>For orthogonal projection, the option of averaging at the vertex appears to give the best results (e.g. flags 3 and 4 in orthogonal projections).</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Of the three methods, the orthogonal one is most robust and clear. It is demonstrated next.

Using cross-cutting for path averaging

The particular advantage of the orthogonal projection is that every query “cuts” adjacent paths as one would expect. 

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaClosestPoint_Ex10.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The particular advantage of the orthogonal projection is that every query “cuts” adjacent paths as one would expect.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_newTraversalByStationResampling

The function, fcn_Path_newTraversalByStationResampling implements a station resampling process

```MATLAB
%% Finding the traversal with the most data, fcn_Path_findTraversalWithMostData
% finds the traversal index with the most amount of data (determined as the
% most elements in the X array)
 
% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
 
% Convert paths into traversals
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral});
    data.traversal{i_traveral} = traversal;
end
 
index_of_longest = fcn_Path_findTraversalWithMostData(data);
fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',...
    length(data.traversal),...
    index_of_longest,...
    length(data.traversal{index_of_longest}.X));

```

<pre align="center">
  <img src=".\Images\fcn_Path_newTraversalByStationResampling.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>


<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

*** 

#### fcn_Path_findAverageTraversalViaOrthoProjection

The function, fcn_Path_findOrthoScatterFromTraversalToTraversals implements this cutting process at given stations

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This allows one to average adjacent paths along each cut to produce an average path. 

This is implemented in fcn_Path_findAverageTraversalViaOrthoProjection.m

The process of averaging has to avoid averaging paths that are nearby but not adjacent in the s-coordinate. Otherwise, the average may bias toward the alternate paths and may even “wander” at the intersections.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

To force a query to include only a short part of a path, one can use: fcn_Path_findPathSXYSegment.m as mentioned earlier.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex3.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

A comparison of the averaging methods shown here reveals that orthogonal projection cross-sections appear to work the best

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex4.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Note that the orthogonal averaging function has the additional outputs of distance and XY hit points of the nearby paths, which allow statistics, detailed point analysis, and advanced plotting.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex5.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

When doing the average, the starting traversal for the average is important as this “seeds” the search. The user can specify this reference_traversal as an input. If nothing is specified, the longest input traversal is used.

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex6.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Because orthogonal averaging uses iterations, the function allows the user to specify the number of iterations to try (default is 40, which is usually too much)

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex7.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

By turning on the debug flag within the function, one can see the iteration results (shown here for N=40)

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex8.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Finally, the averaging process allows weighted smoothing to be used between iterations. A higher weight (from 0 to 1) means less of the new value is used. This prevents numerical oscillation (default is 0.8)

<pre align="center">
  <img src=".\Images\fcn_Path_findAverageTraversalViaOrthoProjection_Ex9.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

*** 

#### fcn_Path_findTraversalWithMostData

The function, fcn_Path_findTraversalWithMostData finds the traversal with the most data (e.g. most elements within an array)

This functionality is used to define the “default” reference traversal when doing orthogonal projection averaging.`

```MATLAB
%% Finding the traversal with the most data, fcn_Path_findTraversalWithMostData
% finds the traversal index with the most amount of data (determined as the
% most elements in the X array)
 
% Fill in some dummy data
paths_array = fcn_Path_fillSamplePaths;
 
% Convert paths into traversals
for i_traveral = 1:length(paths_array)
    traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_traveral});
    data.traversal{i_traveral} = traversal;
end
 
index_of_longest = fcn_Path_findTraversalWithMostData(data);
fprintf(1,'The longest path of the %.0d paths was path %.0d with %.0d elements\n',...
    length(data.traversal),...
    index_of_longest,...
    length(data.traversal{index_of_longest}.X));

```
<pre align="center">
  <img src=".\Images\fcn_Path_findTraversalWithMostData.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_cleanPathFromForwardBackwardJogs

The function, fcn_Path_cleanPathFromForwardBackwardJogs removes forward/backward oscillations caused by averaging data spatially

The main approach therein is to examine the change in relative yaw angle between segments, removing those that are too large.

```MATLAB
%% Example 1: Basic call 

fig_num = 1;
path_with_jogs = [0 0; 1 1; 2 2.2; 3.3 3; 2.5 2.7; 3.5 3.6; 5 5];
clean_path = fcn_Path_cleanPathFromForwardBackwardJogs...
(path_with_jogs,fig_num);
 
plot(path_with_jogs(:,1), path_with_jogs(:,2),'.-','Linewidth',2,'Markersize',25);  
plot(clean_path(:,1), clean_path(:,2),'.-','Linewidth',2,'Markersize',25);
title('Original path with jogs and cleaned path')
xlabel('X [m]')
ylabel('Y [m]')

```

<pre align="center">
  <img src=".\Images\fcn_Path_cleanPathFromForwardBackwardJogs.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_findClosestPointsToTraversal

The function, fcn_Path_findClosestPointsToTraversal implements the method described in the following slides to find the closest points on a traversal to a traversal

<pre align="center">
  <img src=".\Images\fcn_Path_findClosestPointsToTraversal_Ex1.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This method has the advantage of generating expected results even when the queries “graze” nearby paths

<pre align="center">
  <img src=".\Images\fcn_Path_findClosestPointsToTraversal_Ex2.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

But this method gives unexpected results when the segment’s closest location is not actually on the segment.

<pre align="center">
  <img src=".\Images\fcn_Path_findClosestPointsToTraversal_Ex3.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Generally, the “closest point” method gives expected results when applied to adjacent paths

<pre align="center">
  <img src=".\Images\fcn_Path_findClosestPointsToTraversal_Ex4.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### 3D and Elevated Paths

What is an elevated path?

An elevated path is a set of [X Y Z] points as a N x 3 vector or array. These X Y Z rows denote the locations we are trying to follow, in 3D, which defines the elevated path. For example:

```MATLAB
    paths{1} = [
        8.1797    8.6006  0.1
        8.4101   15.8892  0.2
        10.0230   25.5102  -0.1
        10.4839   39.7959  0.3];
```

An elevated path type must have at least 2 rows so that at least one “path segment” can be defined.



<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

#### fcn_Path_addElevationToPath

fcn_Path_addElevationToPath adds elevation to a "path" by finding two nearest neighbors in an "elevated path" and interpolating their elevations

```MATLAB
%% BASIC example 1.3 - works
point = [0.5 0.2; 1.4 1.3]; % Define the query point as an XY
reference_elevated_path = [0 0 0.1; 0.25 0.2 0.2; 0.9 0.9 0.3; 1.1 1.1 0.4; 2.3 2.7 0.5]; % Define an XYZ path
fignum = 113; % Define the figure number
 
% Snap the point onto the path
elevated_path = fcn_Path_addElevationToPath(point, reference_elevated_path, fignum);
```
<pre align="center">
  <img src=".\Images\fcn_Path_addElevationToPath.png" alt="fcn_Path_fillRandomTraversalsAboutTraversal picture" width="400" height="300">
  <figcaption>The function fcn_Path_fillRandomTraversalsAboutTraversal generates random traversals about a given traversal.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

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

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

### Examples

1. Run the main script to set up the workspace and demonstrate main outputs, including the figures included here:

   ```sh
   script_demo_AlignCoordinates
   ```

    This exercises the main function of this code: fcn_AlignCoords_fit2DCoordinates

2. After running the main script to define the included directories for utility functions, one can then navigate to the Functions directory and run any of the functions or scripts there as well. All functions for this library are found in the Functions sub-folder, and each has an associated test script. Run any of the various test scripts, such as:

   ```sh
   script_test_fcn_aligncoords_breakDataIntoLapIndices
   ```

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

## Major release versions

This code is still in development (alpha testing)

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary](https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary)

<a href="#pathplanning_pathtools_pathclasslibrary">Back to top</a>

***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->