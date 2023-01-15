<!--
The following template is based on:
Best-README-Template
Search for this, and you will find!
>
<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h2 align="center"> PathPlanning_PathTools_PathClassLibrary
  </h2>

<img src=".\Images\RaceTrack.jpg" alt="main laps picture" width="960" height="540">

  <p align="center">
    This is the Path class library repo for MATLAB and listing of all functions within this class, whose purpose is to abstract very common "path" operations. What's a "path"? Why do we need them? How do we use them? Read more below for some exciting tools!
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents">View Demo</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Report Bug</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Request Feature</a>
  </p>
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
	    <li><a href="#functions">Functions</li>
	    </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
	    <ul>
	    <li><a href="#examples">Examples</li>
	    <li><a href="#definition-of-endpoints">Definition of Endpoints</li>
	    </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

This "path" class library is one of the first and most commonly used libraries developed within the Intelligent Vehicles and Systems Group, as it contains the set of functions for processing directionality operations related to the "paths" that vehicles follow. This includes calculating intersections, offsets, projections, etc. For example: do you want to know the closest point a vehicle travels to an point in space (e.g. path to point)? Use this library. Or: what is the nearest road centerline to the current vehicle location (e.g. point to path)? Use this library. Or, figure out when a vehicle collides with another vehicle (intersection)? Use this library! Or even, what is the average trajectory of many different vehicles going down the road? You guessed it - use this library.

Most of the tools used herein are abstracted from video gaming algorithms, and thus are highly optimized for speed in an algorithm context, e.g. they are parallelizable, support vectorization of inputs, minimize internal looping, and are highly robust to user input errors and fault cases. Where possible, we try to use matrix forms. All codes are in MATLAB to support generalization into other languages (Python, C, C++, etc.) and some students have started this conversion as these codes are FAR faster. However, MATLAB still seems to be the best "scratchpad" to test algorithms against each other (with Python being great also) - and if an algorithm is poorly chosen, a fast programming language will rarely overcome this flaw.

### Before you start
To learn more about "paths" before starting, please watch the excellent YouTube video, "The Continuity of Splines", by Freya Holmer (her videos are amazing, among the best out there). Her introduction to splines can be found at this link: https://www.youtube.com/watch?v=jvPPXbo87ds&ab_channel=FreyaHolm%C3%A9r. 

Chapter 1 of her video introduces the blending "t" notation of interpolation (used throughout the functions in "Paths"), the "lerp" operation (e.g. linear interpolation), Bezier curve derivations as a function of lerp in quadratic and cubic forms via DeCasteljau's algorithm. It also introduces variants such as the Bernstein form using basis polynomials (which follows Pascal's triangle, see https://math.stackexchange.com/questions/447718/algorithm-for-bezier-curve-with-eleven-control-points) and is based on convex combinations. It also shows how to pre-calculate polynomial coefficients for Bezier curves, and the characteristic matrix form of the Bezier curve (which, by the way, is invertable!). It also introduces the problem of local versus global control of a path and why global control is problematic.

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

Most of the issues discussed in this video sneak into the code in this library (as well as MANY more). So 

























* Inputs: 
    * either a "traversals" type, as explained in the Path library, or a path of XY points in N x 2 format
    * the start, end, and optional excursions can be entered as either a line segment or a point and radius.  
* Outputs
    * Separate arrays of XY points, or of indices for the lap, with one array for each lap
    * The function also can return the points that were not used for laps, e.g. the points before the first start and after the last end

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1.  Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work.

2. Clone the repo
   ```sh
   git clone https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps
   ```
3. Unzip the zip files (DebugTools and PathClassLibrary) into a Utilities folder (.\Utilities), in locations .\Utilities\DebugTools and .\Utilities\PathClassLibrary
4. Confirm it works! Run script_demo_Laps from the root directory root location. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.


<!-- STRUCTURE OF THE REPO -->
### Directories
The following are the top level directories within the repository:
<ul>
	<li>/Documents folder: Descriptions of the functionality and usage of the various MATLAB functions and scripts in the repository.</li>
	<li>/Functions folder: The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
	<li>/Utilities folder: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>
</ul>

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools

* [PathPlanning_PathTools_PathClassLibrary](https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary) - the PathClassLibrary contains tools used to find intersections of the data with particular line segments, which is used to find start/end/excursion locations in the functions. The repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary

    Each should be installed in a folder called "Utilities" under the root folder, namely ./Utilities/DebugTools/ , ./Utilities/PathClassLibrary/ . If you wish to put these codes in different directories, the main call stack in script_demo_Laps can be easily modified with strings specifying the different location, but the user will have to make these edits directly. 
    
    For ease of getting started, the zip files of the directories used - without the .git repo information, to keep them small - are included in this repo.

<!-- FUNCTION DEFINITIONS -->
### Functions
**Basic Support Functions**
<ul>
	<li>
    fcn_Laps_plotLapsXY.m : This function plots the laps. For example, the function was used to make the plot below of the last Sample laps.
    <br>
    <img src=".\Images\fcn_Laps_plotLapsXY.png" alt="fcn_Laps_plotLapsXY picture" width="400" height="300">
    </li>
	<li>
    fcn_Laps_fillSampleLaps.m : This function allows users to create dummy data to test lap functions. The test laps are in general difficult situations, including scenarios where laps loop back onto themself and/or with separate looping structures. These challenges show that the library can work on varying and complicated data sets. NOTE: within this function, commented out typically, there is code to allow users to draw their own lap test cases.
    <br>
    <img src=".\Images\fcn_Laps_fillSampleLaps.png" alt="fcn_Laps_fillSampleLaps picture" width="400" height="300">
    </li>
    <li>
    fcn_Laps_plotZoneDefinition.m : Plots any type of zone, allowing user-defined colors. For example, the figure below shows a radial zone for the start, and a line segment for the end. For the line segment, an arrow is given that indicates which direction the segment must be crossed in order for the condition to be counted. 
    <br>
    <img src=".\Images\fcn_Laps_plotZoneDefinition.png" alt="fcn_Laps_plotZoneDefinition picture" width="400" height="300">
    </li>
    <li>
    fcn_Laps_plotPointZoneDefinition.m : Plots a point zone, allowing user-defined colors. This function is mostly used to support fcn_Laps_plotZoneDefinition.m 
    </li>
    <li>
    fcn_Laps_plotSegmentZoneDefinition.m : Plots a segment zone, allowing user-defined colors. This function is mostly used to support fcn_Laps_plotZoneDefinition.m 
    </li>
    
    

</ul>

**Core Functions**
<ul>
	<li>
    fcn_Laps_breakDataIntoLaps.m : This is the core function for this repo that breaks data into laps. Note: the example shown below uses radial zone definitions, and the results illustrate how a lap, when it is within a start zone, starts at the FIRST point within a start zone. Similarly, each lap ends at the LAST point before exiting the end zone definition. The input data is a traversal type for this particular function.
    <br>
    <img src=".\Images\fcn_Laps_breakDataIntoLaps.png" alt="fcn_Laps_breakDataIntoLaps picture" width="400" height="300">
    </li>	
	<li>
    fcn_Laps_checkZoneType.m : This function supports fcn_Laps_breakDataIntoLaps by checking if the zone definition inputs are either a point or line segment zone specification.
    </li>
	<li>
    fcn_Laps_breakDataIntoLapIndices.m : This is a more advanced version of fcn_Laps_breakDataIntoLaps, where the outputs are the indices that apply to each lap. The input type is also easier to use, a "path" type which is just an array of [X Y]. The example here shows the use of a segment type zone for the start zone, a point-radius type zone for the end zone. The results of this function are the row indices of the data. The plot below illustrates that the function returns 3 laps in this example, and as well returns the pre-lap and post-lap data. One can observe that it is common that the prelap data for one lap (Lap 2) consists of the post-lap data for the prior lap (Lap 1). 
    <br>
    <img src=".\Images\fcn_Laps_breakDataIntoLapIndices.png" alt="fcn_Laps_breakDataIntoLapIndices picture" width="600" height="300">
    </li>	
	<li>
    fcn_Laps_findSegmentZoneStartStop.m : A supporting function that finds the portions of a path that meet a segment zone criteria, returning the starting/ending indices for every crossing of a segment zone. The crossing must cross in the correct direction, and a segment is considered crossed if either the start or end of segment lie on the segment line. This is illustrated in the challenging example shown below, where the input path (thin blue) starts at the top, and then zig-zags repeatedly over a segment definition (green). For each time the blue line crosses the line segment, that portion of the path is saved as a separate possible crossing and thus, for this example, there are 5 possible crossings.
    <br>
    <img src=".\Images\fcn_Laps_findSegmentZoneStartStop.png" alt="fcn_Laps_findSegmentZoneStartStop picture" width="400" height="300">
    </li>	
	<li>
    fcn_Laps_findPointZoneStartStopAndMinimum.m : A supporting function that finds the portions of a path that meet a point zone criteria, returning the starting/ending indices for every crossing of a point zone. Note that a minimum number of points must be within the zone for it to be considered activated, which is useful for real-world data (such as GPS recordings) where noise may randomly push one point of a path randomly into a zone, and then jump out. This number of points threshold can be user-defined. In the example below, the threshold is 4 points and one can see that, for a path that crosses over the zone three times, that two of the crossings are found to meet the 4-point criteria.
    <br>
    <img src=".\Images\fcn_Laps_findPointZoneStartStopAndMinimum.png" alt="fcn_Laps_findPointZoneStartStopAndMinimum picture" width="400" height="300">
    </li>	
</ul>
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

<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

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
For more examples, please refer to the [Documentation](https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents)

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


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


## Major release versions
This code is still in development (alpha testing)


<!-- CONTACT -->
## Contact
Sean Brennan - sbrennan@psu.edu

Project Link: [hhttps://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps](https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[contributors-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[forks-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/network/members
[stars-shield]: https://img.shields.io/github/stars/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[stars-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/stargazers
[issues-shield]: https://img.shields.io/github/issues/ivsg-psu/reFeatureExtraction_Association_PointToPointAssociationpo.svg?style=for-the-badge
[issues-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues
[license-shield]: https://img.shields.io/github/license/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[license-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/blob/master/LICENSE.txt








