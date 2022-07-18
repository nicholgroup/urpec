# urpec
urpec is a matlab-based proximity-effect-correction function for electron beam lithography.

urpec generates a proximity-effect-corrected pattern file for electron beam lithography. The corrected file is created by deconvolving a point spread function from an input .dxf or .mat pattern file. The output file has different colors, each of which receive a different dose. This function assumes that one unit in the input pattern file is one micron. urpec will recognize all closed polylines in all layers. 

urpec can be used with NPGS and can output a .dc2 file suitable for NPGS. Urpec can easily be adapted for use with other electron-beam lithography platforms, such as Elionix. Please contact us for more details about this.

To run urpec, download all files to a folder. Set your matlab path to this folder and run urpec_v4, optionally with arguments. See the documentation inside urpec_v4 for calling urpec_v4 with arguments. You can also run it from a script. This is useful for batch processing. See the script called run_urpec for more details on this and how to create an NPGS run file.

Current PSFs available include the following
- GaAs, 30 kV, 200 nm PMMA
- GaAs, 30 kV, 350 nm PMMA
- GaAs, 30 kV, 400 nm PMMA
- GaAs, 30 kV, 510 nm PMMA
- Si, 30 kV, 200 nm PMMA
- Si with 90 nm SiO2, 30 kV, 200 nm PMMA
- Si with 90 nm SiO2, 50 kV, 200 nm PMMA
- LiNbO3, 30 kV, 200 nm PMMA
- LiNbO3, 50 kV, 200 nm PMMA

If you need a different PSF, you should download Casino, run a simulation, save the data, and run the function casinoPSF2.

Casino can be found here:
http://www.gel.usherbrooke.ca/casino/What.html

An example Casino simulation is included in the "Examples" directory. Open this file in Casino, modify your sample and microscope as necessary, and run the simulation. It should already have the necessary settings enabled, but you need to change the number of electrons. Go to Setting>Set Up Microscope, and change the number of electrons to >2000. Also go to Settings>Run Time Options, and change the number of displayed trajectories to the same number. To change the sample go to Settings>Modify Sample, and change the sample as needed. Note that the top layer of your sample should be called PMMA or MMA in order for the analysis function to work. When the simulation has finished, run casinoPSF2, which extracts point spread function data.

See the script run_urpec for examples on how to run urpec_v4 and how to make run files.

The required Matlab tooboxes include:
- Image processing toolbox
- Statistics and machine learning toolbox
- Curve fitting toolbox
- Mapping toolbox

dxflib is (c) Grzegorz Kwiatek. 
It is distributed here for the sake of completeness. It can be found here:
https://www.mathworks.com/matlabcentral/fileexchange/33884-dxflib

dxf2coord_20 is based on a script of the same name written by Lukas Wischounig. 
It is distributed here for completeness. It can be found here:   
https://www.mathworks.com/matlabcentral/fileexchange/28791-dxf2coord-2-0

fitwrap is (c) Hendrik Bluhm

urpec is written by Adina Ripin, Elliot Connors, and John Nichol.







