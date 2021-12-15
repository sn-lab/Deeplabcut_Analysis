# deeplabcut_analysis


A MATLAB app to analyze .h5 tracking files from videos processed through either the "YMaze" or "BehaviorBox" DeepLabCut trained networks (which can be downloaded [here](https://cornell.box.com/s/wq6kjvpou3k9713ovjpwnh1ghhqiz0ww)). More info on DeepLabCut is available [here](https://deeplabcut.github.io/DeepLabCut/docs/intro.html).

This app simply functions to load a list of .h5 files to analyze, checks their compatibility with the selected analysis type, and processes each file through the appropriate analysis script according to the analysis parameters specified within the app. For each type of analysis, a specific script is used which loads the appropriate tracking labels, cleans up the dataset, and processes/plots the results in the desired way. These are the analysis types currently included in the app:

1. Y-maze test for spontaneous alternation
2. Open-field test for anxiety-like behaviors (i.e. staying close to walls)
3. Object-recognition test to quantify mouse interactions to detected objects (for each .h5 file, the original .avi video file must be in the same path and have the same base filename -- see examples)

&nbsp;

## downloading the app

### If you already have Matlab installed:
Simply download this code (e.g. with Github Desktop) and add the folder to MATLAB's paths (from the MATLAB "HOME" tab: click on "Set Path", "Add with Subfolders", and select the folder where you downloaded this code). Then, open **deeplabcut_analysis.mlapp** to use the app.

### If you don't have MATLAB installed and/or don't have a MATLAB license:
First, download and install the free MATLAB runtime [here](https://www.mathworks.com/products/compiler/matlab-runtime.html) (**must use version 9.10, other versions aren't accepted**). Then, open **deeplabcut_analysis.exe** to use the app.

&nbsp;

## analyzing h5 files
Before using this app, videos must first be processed through DeepLabCut using one of the networks listed at the top. The .h5 file that is output from DeepLabCut contains coordinates of all the tracked labels from the video, and that's what this app primarily uses to perform the analyses. When you start up the app and choose the type of analysis you want to perform from the dropdown box, click the "Add files" button and navigate to/select the .h5 file(s) you want to analyze. Click the "Analyze" button to process the selected files through the analysis script. **Note: for the "object interaction" analysis type, the original .avi video file needs to be in the same folder as the .h5 file you select, and have the same base filename (e.g. "video1.avi" and "video1DeepCut...h5" must be in the same folder)**

&nbsp;

## what if I want to use this app to analyze data from my own DeepLabCut networks instead of the ones listed at the top?
It is possible to train new networks and still use this app, but there are specific requirements in order to make a new network that is still compatible with this app as is.

### Requirements for the ymaze analysis type
In the config.yaml file of a new DLC network, these are the names and labels that must be used verbatim:
>Task: YMaze
>bodyparts:
>- rightarm
>- leftarm
>- middlearm
>- rightear
>- leftear
>- nose
>- tailbase

### Requirements for the open field and object interaction analysis types
In the config.yaml file of a new DLC network, this is the task name and labels that must be used verbatim:
>Task: BehaviorBox
>bodyparts:
>- toprightcorner
>- topleftcorner
>- botleftcorner
>- botrightcorner
>- rightear
>- leftear
>- nose
>- tailbase

## what if my DeepLabCut network doesn't conform to the requirements listed above?
For networks that are incompatible with the current analysis scripts of the app (e.g. if different bodyparts have to be labelled than the ones listed above), a new analysis script would have to be written and incorporated in various places in the app. Create a [new issue](https://github.com/sn-lab/Deeplabcut_Analysis/issues) to request a new analysis type, and we'll see what we can do.
