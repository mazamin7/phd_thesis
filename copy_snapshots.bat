@echo off
echo Copying FDTD snapshots...
xcopy "code\fdtd\scripts\snapshots\experiment_5\*" "latex\Images\fdtd_snapshots\experiment_5\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_6\*" "latex\Images\fdtd_snapshots\experiment_6\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_1_swr\*" "latex\Images\fdtd_snapshots\experiment_1_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_2_swr\*" "latex\Images\fdtd_snapshots\experiment_2_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_3_swr\*" "latex\Images\fdtd_snapshots\experiment_3_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_4_swr\*" "latex\Images\fdtd_snapshots\experiment_4_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_5_swr\*" "latex\Images\fdtd_snapshots\experiment_5_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\experiment_6_swr\*" "latex\Images\fdtd_snapshots\experiment_6_swr\" /Y /I /E
xcopy "code\fdtd\scripts\snapshots\swr_parameters\*" "latex\Images\fdtd_snapshots\swr_parameters\" /Y /I /E

echo Copying PSTD snapshots...
xcopy "code\pstd\scripts\snapshots\experiment_1\*" "latex\Images\pstd_snapshots\experiment_1\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_2\*" "latex\Images\pstd_snapshots\experiment_2\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_3\*" "latex\Images\pstd_snapshots\experiment_3\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_4\*" "latex\Images\pstd_snapshots\experiment_4\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_5\*" "latex\Images\pstd_snapshots\experiment_5\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_6\*" "latex\Images\pstd_snapshots\experiment_6\" /Y /I /E

xcopy "code\pstd\scripts\snapshots\experiment_1_order\*" "latex\Images\pstd_snapshots\experiment_1_order\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_2_order\*" "latex\Images\pstd_snapshots\experiment_2_order\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_3_order\*" "latex\Images\pstd_snapshots\experiment_3_order\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_4_order\*" "latex\Images\pstd_snapshots\experiment_4_order\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_5_order\*" "latex\Images\pstd_snapshots\experiment_5_order\" /Y /I /E
xcopy "code\pstd\scripts\snapshots\experiment_6_order\*" "latex\Images\pstd_snapshots\experiment_6_order\" /Y /I /E

echo Renaming FDTD final_nodal_error.png to final_nodal_error_monolithic.png...
if exist "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error.png" move /Y "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error.png" "latex\Images\fdtd_snapshots\experiment_5\final_nodal_error_monolithic.png"
if exist "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error.png" move /Y "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error.png" "latex\Images\fdtd_snapshots\experiment_6\final_nodal_error_monolithic.png"

echo Snapshots successfully copied!

