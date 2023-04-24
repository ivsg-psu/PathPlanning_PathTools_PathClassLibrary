<!--
The following template is based on:
Best-README-Template
Search for this, and you will find!
>
<!-- PROJECT LOGO -->
<br />
<p align="center">

<h2 align="center"> Errata_Tutorials_DebugTools  </h2>

<img src=".\Images\Debug_Main_Image.jpg" alt="main debug splash picture" width="960" height="540">
<br />
<font size="-2">Photo by Sigmund on Unsplash</font>

<br />
  <p align="left">

# DebugTools
This repo provides common tools used for debugging MATLAB codes within IVSG, and includes input checking, print to workspace, parsing user inputs, and similar functions.
    <br />
    <!-- a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents">View Demo</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Report Bug</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Request Feature</a>
    -->
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ul>
    <li>
      <a href="#about-the-library">About the Library</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#starting-examples">Starting Examples</a></li>
      </ul>
    </li>
    <li><a href="#repo-structure">Repo Structure</a>
	    <ul>
	    <li><a href="#directories">Directories</li>
	    <li><a href="#dependencies">Dependencies</li>
	    <li><a href="#function-structure">Function Structure</li>
	    </ul>
    </li>
    <li><a href="#functions-for-workspace-management">Functions for Workspace Management</a>
	    <ul>
      <li><a href="#automatically-installing-matlab-packages">Automatically installing MATLAB packages</li>      
	    <li><a href="#adding-subdirectories-to-the-path">Adding subdirectories to the pathTop-Level Directories</li>
	    </ul>
    </li>
    <li><a href="#functions-for-input-checking">Functions for Input Checking</a>
	    <ul>
	    <li><a href="#checking-inputs-to-functions">Checking inputs to functions</li>
        <li><a href="#checking-if-strings-partially-match">Checking if strings partially match</li>    
        <li><a href="#extracting-a-numeric-value-embedded-in-a-string">Extracting a numeric value embedded in a string </li>
        <li><a href="#converting-mixed-input-lists-into-comma-separated-string">Converting mixed input lists into comma sparated string</li>            
	    </ul>
    </li>
    <li><a href="#functions-for-output-formatting">Functions for Output Formatting</a>
	    <ul>
	    <li><a href="#converting-numbers-to-human-friendly-strings">Converting numbers to human friendly strings</li>
        <li><a href="#appending-arbitrary-values-to-a-string">Appending arbitrary values to a string</li>    
        <li><a href="#printing-results-to-fixed-length-strings">Printing results to fixed length strings</li>  
        <li><a href="#printing-matrices-to-fixed-length-columns">Printing matrices to fixed length columns</li>           
	    </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
	    <ul>
	    <li><a href="#examples">Examples</li>
	    <li><a href="#definition-of-endpoints">Definition of Endpoints</li>
	    </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ul>
</details>

<a href="#debugtools">Back to top</a>

<!-- ABOUT THE LIBRARY -->
## About the Library

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

Often the codes used within IVSG require very common, repeated tools to be used for debugging or general error catching. To avoid repitition of code, most of these functions are embedded here.

The general areas of functionality include:

* Workspace management
    * Path checking and creation
* Input checking: 
    * Testing whether inputs meet specifications required by code, such as whether or not they are integers, positive, contain 2 or 3 rows, only have 2 columns, etc.  
* Output formatting
    * Printing results to screen with a specified width.

<a href="#debugtools">Back to top</a>

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1.  Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work.

2. Clone the repo
```sh
git clone https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
```
3. Confirm it works! Run script_demo_DebugTools.m from the root directory root location. If the code works, the script should run without errors. This script produces numerous example cases, many of them shown in this README file.

<a href="#debugtools">Back to top</a>

### Starting Examples

1. Run the main script to set up the workspace and demonstrate main outputs, including the figures included in this README:

   ```sh
   script_demo_DebugTools
   ```

2. After running the main script to define the included directories for utility functions, one can then navigate to the Functions directory and run any of the functions or scripts there as well. All functions for this library are found in the Functions sub-folder, and each has an associated test script. Run any of the various test scripts, such as:

   ```sh
   script_test_fcn_DebugTools_extractNumberFromStringCell
   ```
For more examples, please refer to the [Documentation](https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents)

<a href="#debugtools">Back to top</a>

<!-- REPO STRUCTURE-->

## Repo Structure

### Directories
The following are the top level directories within the repository:
<ul>
	<li>/Data folder: Contains example data used for demonstrating the debugging tools.</li>
	<li>/Documents folder: Presentations of best practices for debugging - many of which were made by prior students.</li>
    	<li>/Example Code Snippets folder: Example codes demonstrating some very powerful capabilities in MATLAB including how to set up a class, how to set environmental variables, how to perform unit testing, and how to run multiple test scripts with logging.</li>
	<li>/Functions folder: The majority of the code functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
    	<li>/Images folder: Location of images used in this and or other README files.</li>
	<li>/Utilities folder: (empty) Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>
</ul>

<a href="#debugtools">Back to top</a>

### Dependencies

* (none)

<!--
* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools

* [PathPlanning_PathTools_PathClassLibrary](https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary) - the PathClassLibrary contains tools used to find intersections of the data with particular line segments, which is used to find start/end/excursion locations in the functions. The repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary

    Each should be installed in a folder called "Utilities" under the root folder, namely ./Utilities/DebugTools/ , ./Utilities/PathClassLibrary/ . If you wish to put these codes in different directories, the main call stack in script_demo_Laps can be easily modified with strings specifying the different location, but the user will have to make these edits directly. 
    
    For ease of getting started, the zip files of the directories used - without the .git repo information, to keep them small - are included in this repo.
-->

<a href="#debugtools">Back to top</a>

<!-- FUNCTION STRUCTURE -->
### Function Structure
Each of the functions described as follows has a consistent structure: each has an associated test script, using the convention

```sh
script_test_fcn_fcnname
```
where fcnname is the function name as listed above.

Also, each of the functions includes a well-documented header that explains inputs and outputs. These are supported by MATLAB's help style so that one can type:

```sh
help fcn_fcnname
```
for any function to view function details.

<a href="#debugtools">Back to top</a>

<!-- FUNCTIONS FOR WORKSPACE MANAGEMENT -->
## Functions for Workspace Management
### Automatically installing MATLAB packages
fcn_DebugTools_installDependencies.m : installs code packages that are specified by a URL pointing to a zip file into a default local subfolder, "Utilities", under the root folder. It also adds either the package subfoder or any specified sub-subfolders to the MATLAB path.

If the Utilities folder does not exist, it is created.

If the specified code package folder and all subfolders already exist, the package is not installed. Otherwise, the folders are created as needed, and the package is installed.

If one does not wish to put these codes in different directories, the function can be easily modified with strings specifying the desired install location.

For path creation, if the "DebugTools" package is being installed, the code installs the package, then shifts temporarily into the package to complete the path definitions for MATLAB. If the DebugTools is not already installed, an error is thrown as these tools are needed for the path creation.

Finally, the code sets a global flag to indicate that the folders are initialized so that, in this session, if the code is called again the folders will not be installed. This global flag can be overwritten by an optional flag input.


```Matlab
%% Basic test case
% NOTE: this installs under the current directory!
% Define the name of subfolder to be created in "Utilities" subfolder
dependency_name = 'DebugTools_v2023_01_18';

% Define sub-subfolders that are in the code package that also need to be
% added to the MATLAB path after install. Leave empty ({}) to only add
% the subfolder path without any sub-subfolder path additions.
dependency_subfolders = {'Functions','Data'};

% Define a universal resource locator (URL) pointing to the zip file to
% install. For example, here is the zip file location to the Debugtools
% package on GitHub:
dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_25.zip?raw=true';

% Call the function to do the install
fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
```

### Adding subdirectories to the path

fcn_DebugTools_addSubdirectoriesToPath.m : This function adds given subdirectories to the root path, and causes an error to be thrown if the directory is not found. It is typically used within a flag-checking if-statement, such as below

```Matlab
%% Demonstrate how to add subdirectories
if ~exist('flag_DebugTools_Folders_Initialized','var')
    fcn_DebugTools_addSubdirectoriesToPath(pwd,{'Functions','Data'});

    % set a flag so we do not have to do this again
    flag_DebugTools_Folders_Initialized = 1;
end
```
Note, the first time code is run, the location of this function is not going to be known. So typically, the initialization codes must be hard-coded first. The first part of the demo code for the DebugTools library does this, adn the code below is from the Laps class library to demonstrate how it is done for codes that are using DebugTools as a utility, which is more typical.

```Matlab
%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%%
% 
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
% * PathClassLibrary - the repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
% 
% Each should be installed in a folder called "Utilities" under the root
% folder, namely ./Utilities/DebugTools/ , ./Utilities/PathClassLibrary/ .
% If you wish to put these codes in different directories, the function
% below can be easily modified with strings specifying the different
% location.
% 
% For ease of transfer, zip files of the directories used - without the
% .git repo information, to keep them small - are included in this repo.
% 
% The following code checks to see if the folders flag has been
% initialized, and if not, it calls the DebugTools function that loads the
% path variables. It then loads the PathClassLibrary functions as well.
% Note that the PathClass Library also has sub-utilities that are included.
if ~exist('flag_Laps_Folders_Initialized','var')
    
    % add necessary directories for function creation utility 
    %(special case because folders not added yet)
    debug_utility_folder = fullfile(pwd, 'Utilities', 'DebugTools');
    debug_utility_function_folder = fullfile(pwd, 'Utilities', 'DebugTools','Functions');
    debug_utility_folder_inclusion_script = fullfile(pwd, 'Utilities', 'DebugTools','Functions','fcn_DebugTools_addSubdirectoriesToPath.m');
    if(exist(debug_utility_folder_inclusion_script,'file'))
        current_location = pwd;
        cd(debug_utility_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(debug_utility_folder,{'Functions','Data'});
        cd(current_location);
    else % Throw an error?
        error('The necessary utilities are not found. Please add them (see README.md) and run again.');
    end
    
    % Now can add the Path Class Library automatically
    utility_folder_PathClassLibrary = fullfile(pwd, 'Utilities', 'PathClassLibrary');
    fcn_DebugTools_addSubdirectoriesToPath(utility_folder_PathClassLibrary,{'Functions','Utilities'});
    
    % utility_folder_GetUserInputPath = fullfile(pwd, 'Utilities', 'GetUserInputPath');
    % fcn_DebugTools_addSubdirectoriesToPath(utility_folder_GetUserInputPath,{'Functions','Utilities'});

    % Now can add all the other utilities automatically
    folder_LapsClassLibrary = fullfile(pwd);
    fcn_DebugTools_addSubdirectoriesToPath(folder_LapsClassLibrary,{'Functions'});

    % set a flag so we do not have to do this again
    flag_Laps_Folders_Initialized = 1;
end
```

<!--img src=".\Images\fcn_Laps_plotLapsXY.png" alt="fcn_Laps_plotLapsXY picture" width="400" height="300"
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
  -->  
<a href="#debugtools">Back to top</a>    

## Functions for Input Checking

### Checking inputs to functions
<ul>
	<li> The function: fcn_DebugTools_checkInputsToFunctions checks to see if an input meets specified requirements. It is a VERY powerful tool that is commonly used at the top of most codes in the IVSG toolset. For example, the following code checks to see if the input has 2 columns, and 5 or less rows. If it does, it gives no error:

```Matlab
%% Check that an input has 2 columns, 
% Maximum length is 5 or less
Twocolumn_of_integers_test = [4 1; 3 9; 2 7];
fcn_DebugTools_checkInputsToFunctions(Twocolumn_of_integers_test, '2column_of_integers',[5 4]);
```

This function is quite flexible. To see all options, run the script file. An output option is also available that lists all possible inputs. To see these, run the following code:

```Matlab
options = fcn_DebugTools_checkInputsToFunctions;
fprintf(1,'Here are a listing of all active input checking options: \n');
for ith_option = 1:length(options)
    fprintf('\t"%s"\n',options(ith_option).name)
    fprintf('\t\t%s\n',options(ith_option).description)    
end
```

and it returns:

```
Here are a listing of all active input checking options: 
	"Mcolumn_of..."
		checks that the input type is K x M of ...
	"NorMcolumn..."
		checks that the input type is of minimum K x M or maximum K x N of ...
	"positive_..."
		checks that the input type is positive ...
	"_of_integers..."
		checks that the input type is of_integers...
	"_of_mixed..."
		checks that the input type is numeric but can include NaN..
	"_of_chars..."
		checks that the input type is a char type (uses ischar)
	"_of_strings..."
		checks that the input type is a string type (uses isstring)
	"DoesFileExist..."
		checks that the input type is an existing file
	"DoesDirectoryExist..."
		checks that the input type is an existing directory
	"polytopes..."
		checks that the input type is a polytope type, e.g. a structure with fields: vertices, xv, yv, distances, mean, area, max_radius.
	"mixedset..."
		checks that the input type is a structure matching a given template..
	"1column_of_numbers"
		checks that the input type is a structure with fields: name, settings, AABB.
	"positive_1column_of_numbers"
		checks that the input type is N x 1 and is a strictly positive number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N.
	"2column_of_numbers"
		checks that the input type is N x 2 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.
	"4column_of_numbers"
		checks that the input type is N x 4 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.
	"2or3column_of_numbers"
		checks that the input type is N x 2 or N x 3 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.
	"2column_of_integers"
		checks that the input type is N x 2 and is an integer. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.
	"polytopes"
		a 1-by-n seven field structure of polytopes within the boundaries, where n <= number of polytopes with fields:  vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is the number of the individual polytope vertices  xv: a 1-by-m vector of vertice x-coordinates  yv: a 1-by-m vector of vertice y-coordinates  distances: a 1-by-m vector of perimeter distances from one point to the next point, distances(i) = distance from vertices(i) to vertices(i+1) mean: centroid xy coordinate of the polytope area: area of the polytope max_radius: the largest distance from the centroid to any vertex
	"station"
		Path library type: checks that the station type is N x 1 and is a number.
	"stations"
		Path library type: checks that the station type is N x 1 and is a number, with N >= 2.
	"path"
		Path library type: checks that the path type is N x 2 with N>=2
	"path2or3D"
		Path library type: checks that the path type is N x 2 or N x 3, with N>=2
	"elevated_path"
		Path library type: checks that the elevated path type is N x 3 with N>=2
	"paths"
		Path library type: checks that the path type is N x 2 with N>=3
	"traversal"
		Path library type: checks if a structure with X, Y, and Station, and that each has an N x 1 vector within all of same length. Further, the Station field must be strictly increasing.
	"traversals"
		Path library type: checks if a structure containing a subfield that is a cell array of traveral{i}, e.g. "data.traversal{3}", with each traversal also meeting traversal requirements.
	"likestructure"
		Takes a structure input as the 3rd argument to serve as a template. Ensures that the input has the same structure fields.
```
</li>	
</ul>

<a href="#debugtools">Back to top</a>

### Checking if strings partially match
<ul>
<li> The function: fcn_DebugTools_doStringsMatch checks to see a string entry, usually from a prompt to a human user, matches or partially matches a given set of inputs, for example if "y" matches "Yes". This function can also handle cases where the user must select from a set of choices (a through d, for example) and confirms that one and only one of those choices was selected:

```Matlab
%% simple string comparisons, student answer is part of correct answer so returns true, ignoring case
student_answer = 'A';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result);

%% simple string comparisons, student answer is part of correct answer so true, checking to produce false result if student repeats (FALSE)
student_answer = 'aa';
correct_answers = 'abc';
result = fcn_DebugTools_doStringsMatch(student_answer,correct_answers);
assert(result==false);
```
</li>
</ul>

<a href="#debugtools">Back to top</a>

### Extracting a numeric value embedded in a string
<ul>
<li> The function fcn_DebugTools_extractNumberFromStringCell takes a char type in a cell, and finds the first number within that is numeric, and then returns the string for this number. It is robust in that weird entries also work, such as '-0004.2'. This function is particularly useful to parse human-input numbers.
	
```Matlab
%% Decimal number, negative, in cell array with leading zeros and text
result = fcn_DebugTools_extractNumberFromStringCell({'My number is -0000.4'});
assert(isequal(result,{'-0.4'}));
```
	
</li>
</ul>

<a href="#debugtools">Back to top</a>

### Extracting a numeric value embedded in a string
<ul>
<li> The function fcn_DebugTools_parseStringIntoCells parses a string containing comma-separated elements, parsing out the elements into cells.

```Matlab
%% Demonstrate fcn_DebugTools_parseStringIntoCells
% Choose a very Complex input
inputString = 'This,isatest,of';
result = fcn_DebugTools_parseStringIntoCells(inputString);
assert(isequal(result,[{'This'},{'isatest'},{'of'}]));
```

</li>
</ul>

<a href="#debugtools">Back to top</a>

### Converting mixed input lists into comma separated string
<ul>
<li> The function fcn_DebugTools_parseStringIntoCells parses a string containing comma-separated elements, parsing out the elements into cells.

```Matlab
%% Demonstrate fcn_DebugTools_convertVariableToCellString
% Multiple mixed character, numeric in cell array ending in string with commas
result = fcn_DebugTools_convertVariableToCellString([{'D'},{2},'abc , 123']);
assert(isequal(result,{'D, 2, abc , 123'}));
```

</li>
</ul>

<a href="#debugtools">Back to top</a>

## Functions for Output Formatting

### Converting numbers to human friendly strings
<ul>
<li>
The function: fcn_DebugTools_number2string.m prints a "pretty" version of a string, e.g avoiding weirdly odd numbers of decimal places or strangely formatted printing.

```Matlab
%% Basic case - example
stringNumber = fcn_DebugTools_number2string(2.333333333); % Empty result
assert(isequal(stringNumber,'2.33'));
fprintf(1,'%s\n',stringNumber);
```

Produces the following result:

```
2.33
```

</li>
</ul>

<a href="#debugtools">Back to top</a>

### Appending arbitrary values to a string
<ul>
<li>
The function: fcn_DebugTools_addStringToEnd.m appends a number, cell string, or string to the end of a string, producing a string.

```Matlab
input_string = 'test';
value_to_add = 2;
output_string = fcn_DebugTools_addStringToEnd(input_string,value_to_add);
assert(isequal(output_string,'test 2'));
```

Produces the following result:

```
test 2
```

</li>
</ul>

<a href="#debugtools">Back to top</a>

### Printing results to fixed length strings
<ul>
<li>
    The function: fcn_DebugTools_debugPrintStringToNCharacters converts strings into fixed-length forms, so that they print cleanly. For example, the following 2 basic examples:

```Matlab
% BASIC example 1 - string is too long
test_string = 'This is a really, really, really long string but we only want the first 10 characters';
fixed_length_string = fcn_DebugTools_debugPrintStringToNCharacters(test_string,10);
fprintf(1,'The string: %s\nwas converted to: "%s"\n',test_string,fixed_length_string);

% BASIC example 2 - string is too short
test_string = 'Tiny string but should be 40 chars';
fixed_length_string = fcn_DebugTools_debugPrintStringToNCharacters(test_string,40);
fprintf(1,'The string: %s\nwas converted to: "%s"\n',test_string,fixed_length_string);
```
Produces the following output:
```
The string: This is a really, really, really long string but we only want the first 10 characters
was converted to: "This is a "
The string: Tiny string but should be 40 chars
was converted to: "Tiny string but should be 40 chars      "
```
</li>
</ul>

<a href="#debugtools">Back to top</a>

### Printing matrices to fixed length columns
<ul>
<li>
    The function: fcn_DebugTools_debugPrintTableToNCharacters, given a matrix of data, prints the data in user-specified width to the workspace.

```Matlab
%% Fill in test data
Npoints = 10;
point_IDs = (1:Npoints)';
intersection_points = rand(Npoints,2);
s_coordinates_in_traversal_1 = rand(Npoints,1);
s_coordinates_in_traversal_2 = 1000*rand(Npoints,1);
table_data = [point_IDs, intersection_points, s_coordinates_in_traversal_1, s_coordinates_in_traversal_2];

%% Basic test case - constant column widths
header_strings = [{'Data ID'}, {'Location X'},{'Location Y'},{'s-coord 1'},{'s-coord 2'}];
formatter_strings = [{'%.0d'},{'%.12f'},{'%.12f'},{'%.12f'},{'%.12f'}];
N_chars = 15; % All columns have same number of characters
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings,N_chars);


%% Advanced test case - variable column widths
header_strings = [{'Data ID'}, {'Location X'},{'Location Y'},{'s-coord 1'},{'s-coord 2'}]; % Headers for each column
formatter_strings = [{'%.0d'},{'%.12f'},{'%.12f'},{'%.12f'},{'%.12f'}]; % How should each column be printed?
N_chars = [4, 15, 15, 5, 5]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings,N_chars);
```
Produces the following output:
```
Data ID         Location X      Location Y      s-coord 1       s-coord 2       
1               0.083497561789  0.466512012381  0.872812993043  231.09536311113 
2               0.279828920000  0.498104487123  0.938002006920  527.43396277663 
3               0.447007309335  0.487430654415  0.139689278548  724.99196135723 
4               0.587571263995  0.229468774748  0.393900144851  607.41579114493 
5               0.877634138415  0.085552232409  0.980562829715  588.36644579525 
6               0.469100520229  0.067383313609  0.644794025985  433.43484003305 
7               0.437418475702  0.888390934805  0.896409779454  244.17289350702 
8               0.746184939975  0.233167685670  0.482230405437  428.96035377800 
9               0.467910465808  0.861595759984  0.014093075189  10.177455521777 
10              0.860827351058  0.711735093008  0.622880344435  608.82144904436 


Data Location X      Location Y      s-coo s-coo 
1    0.083497561789  0.466512012381  0.872 231.0 
2    0.279828920000  0.498104487123  0.938 527.4 
3    0.447007309335  0.487430654415  0.139 724.9 
4    0.587571263995  0.229468774748  0.393 607.4 
5    0.877634138415  0.085552232409  0.980 588.3 
6    0.469100520229  0.067383313609  0.644 433.4 
7    0.437418475702  0.888390934805  0.896 244.1 
8    0.746184939975  0.233167685670  0.482 428.9 
9    0.467910465808  0.861595759984  0.014 10.17 
10   0.860827351058  0.711735093008  0.622 608.8 
```
</li>
</ul>
<a href="#debugtools">Back to top</a>


<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->



<a href="#debugtools">Back to top</a>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#debugtools">Back to top</a>

## Major release versions

This code is still in development (alpha testing). 
The 2023-01-25 release includes the following addition:
* Updated README
* Added fcn_DebugTools_installDependencies to support automated URL-referenced install of code packages

The 2023-01-18 release includes the following addition:
* Updated README
* Adding directory and file queries for input checking
* String output functions (fixed-length printing)
* Input functions (string parsing)
* Table formatted fixed-column-width output


<a href="#debugtools">Back to top</a>

<!-- CONTACT -->
## Contact
Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools)

<a href="#debugtools">Back to top</a>

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








