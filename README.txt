You need the following toolboxes to run all the files
- Aerospace Toolbox
- Navigation Toolbox
- Mapping Toolbox

Right click Sensors, Helpers and Full Loop Subsystems in MATLAB Files Tree --> Add to Path --> Selected Folder(s) and Subfolders
Should be able to run everything.


Changes to still make
- LiDAR characterisation must be done
- Consider all the update rates of each model of the sensors and the actual update steps in EKF
- Varying sigma for the low-g accelerometer to taper from 1 to 4
- Confirm all the sigma values and the LPF cutoff values for all the sensors in EKF
- Model actuators with step functions on servos and motor
- Controller model verification and clean up