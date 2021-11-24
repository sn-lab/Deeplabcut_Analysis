# Deeplabcut_Analysis



A MATLAB app to analyze .h5 tracking files from DeepLabCut. Can be run as a standalone app on computers which don't have a MATLAB license; if MATLAB isn't installed, download and install the MATLAB runtime here: [Install and Configure MATLAB Runtime - MATLAB & Simulink (mathworks.com)](https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html).

The app simply functions to load a list of files to analyze, checks their compatibility with the selected analysis type, and processes each file through the appropriate analysis script. For each type of analysis, a specific script is needed which loads the appropriate tracking labels, cleans up the data, and processes/plots the results in the desired way. These are the analysis types currently included in the app:

1. Y-maze test for spontaneous alternation
2. Open-field test for anxiety-like behaviors (i.e. walking  next to walls)
3. Novel object-recognition test to quantify mouse interactions to detected objects (for each h5 file, the original avi file must be in the same path and have the same base filename)

- To be added: mouse wheel walking/running gait analysis

