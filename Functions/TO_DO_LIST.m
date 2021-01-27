TO_DO List

script_main_laneChange.m - clean this up by taking Guangwei's functions and adding them into these functions too

LATER:
add more test cases to script_test_fcn_Path_cleanPathFromForwardBackwardJogs
figure out a search distance to find loops (see ortho averaging code)
Weighted averaging of paths (so we can add paths collected later)
a function that finds where intersections and lange changes occur within many paths by looking at change in yaw relative to mean?
a function to find where paths diverage?
pull out the spatial smoothing of offsets within function: fcn_Path_fillRandomTraversalsAboutTraversal into another stand-alone function: fcn_Path_filterOffsetsViaSpatialLowPassFilter(offsets,samplingDistance, smoothingDistance)
clean up fcn_Path_findAverageTraversalViaStationAlignment to match the inputs/header/style of ortho and closest functions
can we define Nth yaw angle same as N-1 th yaw angle? <-- maybe: not clear how to define yaw angle when at end of a traversal! (e.g., for N segments, there is really only N-1 "real" yaw angles)