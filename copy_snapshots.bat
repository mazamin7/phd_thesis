@echo off
echo Copying FDTD snapshots...
xcopy "code\fdtd\scripts\snapshots\experiment_5\*" "latex\Images\fdtd_snapshots\experiment_5\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_6\*" "latex\Images\fdtd_snapshots\experiment_6\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\swr_parameters\*" "latex\Images\fdtd_snapshots\swr_parameters\" /Y /I /E

echo Copying PSTD snapshots...
xcopy "code\pstd\scripts\snapshots\experiment_5\*" "latex\Images\pstd_snapshots\experiment_5\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_6\*" "latex\Images\pstd_snapshots\experiment_6\" /Y /I /E

echo Renaming FDTD final_nodal_error.png to final_nodal_error_monolithic.png...
if exist "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error.png" move /Y "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error.png" "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error_monolithic.png"
if exist "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error.png" move /Y "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error.png" "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error_monolithic.png"

echo Snapshots successfully copied!
pause
