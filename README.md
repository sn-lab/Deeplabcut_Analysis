# deeplabcut_analysis


A MATLAB app to analyze .h5 tracking files from videos processed through either the "YMaze" or "BehaviorBox" DeepLabCut trained networks (which can be downloaded [here](https://cornell.box.com/s/wq6kjvpou3k9713ovjpwnh1ghhqiz0ww)). More info on DeepLabCut is available [here](https://deeplabcut.github.io/DeepLabCut/docs/intro.html).

This app simply functions to load a list of .h5 files to analyze, checks their compatibility with the selected analysis type, and processes each file through the appropriate analysis script according to the analysis parameters specified within the app. For each type of analysis, a specific script is used which loads the appropriate tracking labels, cleans up the dataset, and processes/plots the results in the desired way. These are the analysis types currently included in the app:

1. Y-maze test for spontaneous alternation
2. Open-field test for anxiety-like behaviors (i.e. staying close to walls)
3. Object-recognition test to quantify mouse interactions to detected objects (for each .h5 file, the original .avi video file must be in the same path and have the same base filename -- see examples)

(\)

## getting started

### If you already have Matlab installed:
Simply download this code (e.g. with Github Desktop) and add the folder to MATLAB's paths (from the MATLAB "HOME" tab: click on "Set Path", "Add with Subfolders", and select the folder where you downloaded this code). Then, open **deeplabcut_analysis.mlapp** to use the app.

### If you don't have MATLAB installed and/or don't have a MATLAB license:
First, download and install the free MATLAB runtime [here](https://www.mathworks.com/products/compiler/matlab-runtime.html) (**must use version 9.10, other versions aren't accepted**). Then, open **deeplabcut_analysis.exe** to use the app.


