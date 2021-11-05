@echo on
echo Starting DeepLabCut...
echo -------------------------------------------------------
echo DO NOT CLOSE THIS WINDOW (closing will exit DeepLabCut)
echo -------------------------------------------------------
call C:\ProgramData\Anaconda3\Scripts\activate.bat
call conda activate C:\Users\Schaf\anaconda3\envs\DEEPLABCUT
python C:\Users\Schaf\Documents\GitHub\Deeplabcut_Analysis\misc\start_dlc.py
