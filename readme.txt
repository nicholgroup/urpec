2019/03/07

urpec is a matlab-based proximity-effect-correction function for electron beam lithography.

urpec generates a proximity-effect-corrected .dxf file for electron beam lithography. The corrected file is created by deconvolving a point spread function from an input .dxf pattern file. The output file has different layers, each of which should receive a different dose. This function assumes that one unit in the input pattern file is one micron. urpec will recognize all closed polylines in all layers.

Right now this is intended for use with NPGS.

To run urpec, set your matlab path to Nichol Group\matlab\urpec and run urpec, optionally with arguments. See the documentation inside urpec for calling urpec with arguments. 

Current PSFs available are 
- GaAs, 30 kV, 200nm PMMA
- GaAs, 30 kV, 400 nm PMMA
- GaAs, 30 kV, 510 nm PMMA
- Si, 30 kV, 200 nm PMMA

If you need a different PSF, you should download Casino, run a simulation, save the data, and run the function casinoPSF.

Casino can be found here:
http://www.gel.usherbrooke.ca/casino/What.html

An example Casino simulation is included in the "Example stuff" directory. Open this file in Casino, modify your sample and microscope as necessary, and run the simulation. It should already have the necessary settings enabled. The top layer of your sample should be called PMMA in order for the analysis function to work. Set up the microscope to your desired voltage and make sure to simulate 2000 electrons. Make sure that your display is also set to display 2000 electrons. When the simulation has finished, save the data. Then open the file in excel, and resave as a .xlsx file. The resulting .xlsx file is compatible with the function casinoPSF, which extracts point spread function data. 

dxflib is (c) Grzegorz Kwiatek. 
It is distributed here for the sake of completeness. It can be found here:
https://www.mathworks.com/matlabcentral/fileexchange/33884-dxflib

dxf2coord_20 is written by Lukas Wischounig. 
It is distributed here for completeness. It can be found here:   
https://www.mathworks.com/matlabcentral/fileexchange/28791-dxf2coord-2-0

fitwrap is (c) Hendrik Bluhm

TODO:
1. Add the capability to generate an NPGS run file.
2. Add support for other file formats, like GDSII.
3. Make sure all layers are populated, even if they have nothing in them. Will this work with npgs?
4. The layer boundaries don't touch. The are separated by half the pixel size. Not a huge deal, but it would be nice to fix.
5. Be smart about the grid size
6. Be smart about the subfield size

templates test
