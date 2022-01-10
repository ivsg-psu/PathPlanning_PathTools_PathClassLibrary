TO_DO List

script_main_laneChange.m - clean this up by taking Guangwei's functions and adding them into these functions too

LATER:
add more test cases to script_test_fcn_Path_cleanPathFromForwardBackwardJogs
Weighted averaging of paths (so we can add paths collected later)
a function to find where paths diverage, based on user-defined metric (e.g. a lane width)
a function that finds where intersections and lange changes occur within many paths by looking at change in yaw relative to mean?
pull out the spatial smoothing of offsets within function: fcn_Path_fillRandomTraversalsAboutTraversal into another stand-alone function: fcn_Path_filterOffsetsViaSpatialLowPassFilter(offsets,samplingDistance, smoothingDistance)
clean up fcn_Path_findAverageTraversalViaStationAlignment to match the inputs/header/style of ortho and closest functions. This should be a very simple function now that station resampling is done
