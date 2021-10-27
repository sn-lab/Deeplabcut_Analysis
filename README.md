# Deeplabcut_Analysis



A MATLAB app (and scripts) to analyze .csv tracking files from deeplabcut. 

The app simply functions to load a list of files to analyze, checks their compatibility with the selected analysis type, and processes each file through the appropriate analysis script.

For each type of analysis, a specific script is needed which loads the appropriate tracking labels, cleans up the data, and processes/plots the results in the desired way. These are the types of analysis currently included:

1. Y-maze test for spontaneous alternation
2. Open-field test for anxiety-like behaviors (i.e. walking  next to walls)
3. Novel object-recognition test

- Not yet added: mouse wheel walking/running gait analysis
