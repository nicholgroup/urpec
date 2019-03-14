# urpec
urpec is a matlab-based proximity-effect-correction function for electron beam lithography.

urpec generates a proximity-effect-corrected .dxf file for electron beam lithography. The corrected file is created by deconvolving a point spread function from an input .dxf pattern file. The output file has different layers, each of which should receive a different dose. This function assumes that one unit in the input pattern file is one micron. urpec will recognize all closed polylines in all layers.

Right now this is intended for use with NPGS.

To run urpec, download all files to a folder. Set you matlab path to this folder and run urpec, optionally with arguments. See the documentation inside urpec for calling urpec with arguments. You can also run it from a script. This is useful for batch processing. See the script called run_urpec.

Current PSFs available are 
- GaAs, 30 kV, 200nm PMMA
- GaAs, 30 kV, 400 nm PMMA
- GaAs, 30 kV, 510 nm PMMA
- Si, 30 kV, 200 nm PMMA

If you need a different PSF, you should download Casino, run a simulation, save the data, and run the function casinoPSF.

Casino can be found here:
http://www.gel.usherbrooke.ca/casino/What.html

An example Casino simulation is included in the "Examples" directory. Open this file in Casino, modify your sample and microscope as necessary, and run the simulation. It should already have the necessary settings enabled, but you need to change the number of electrons. Go to Setting>Set Up Microscope, and change the number of electrons to 2000. Also go to Settings>Run Time Options, and change the number of displayed trajectories fo 2000. To change the sample go to Settings>Modify Sample, and change the sample as needed. Note that the top layer of your sample should be called PMMA in order for the analysis function to work. When the simulation has finished, save the data. Then open the file in excel, and resave as a .xlsx file. The resulting .xlsx file is compatible with the function casinoPSF, which extracts point spread function data.

dxflib is (c) Grzegorz Kwiatek. 
It is distributed here for the sake of completeness. It can be found here:
https://www.mathworks.com/matlabcentral/fileexchange/33884-dxflib

dxf2coord_20 is based on a script of the same name written by Lukas Wischounig. 
It is distributed here for completeness. It can be found here:   
https://www.mathworks.com/matlabcentral/fileexchange/28791-dxf2coord-2-0

fitwrap is (c) Hendrik Bluhm

urpec is written by Adina Ripin, Elliot Connors, and John Nichol.





