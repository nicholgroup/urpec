# urpec
Matlab-based proximity effect correction software for electron beam lithography

urpec generates a proximity-effect-corrected .dxf file for electron beam lithography. The corrected file is created by deconvolving a point spread function from an input .dxf pattern file. The output file has different layers, each of which should receive a different dose. This function assumes that one unit in the input pattern file is one micron. urpec will recognize all closed polylines in all layers.

Right now this is intended for use with NPGS.

To run urpec, download all files into a folder and set your matlab path to that folder. Once this is done, run the function urpec, optionally with arguments.

Current PSFs available are for use with Si and GaAs substrates at 30 kV. If you need a different PSF, you need to download Casino, run a simulation, save the data, and run the function casinoPSF.

Casino can be found here:
http://www.gel.usherbrooke.ca/casino/What.html

TODO:
Add the capability to generate an NPGS run file.
Add support for other file formats?




