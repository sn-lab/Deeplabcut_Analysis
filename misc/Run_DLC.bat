@echo on
%windir%\System32\cmd.exe C:\ProgramData\Anaconda3\Scripts\activate.bat C:\ProgramData\Anaconda3 cmd /k 
cd DeepLabCut
conda activate DEEPLABCUT
python
import deeplabcut
deeplabcut.launch_dlc()