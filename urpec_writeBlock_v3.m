function [ s ] = writeBlock( config )
%[ s ] = writeBlock( config )
%   Generates text for an NPGS run file for pattern writing.
%
% Version history
% v2: doesn't exist
% v3: belongs with urpec_v3. Assumes that all patterns are in the same
% layer but with different colors.
curDir=pwd;

cd(config(1).dir)

switch config.aperture
    case 30
        current = '438.0';
        aperture = '1';
        
    case 7.5
        current = '17.0';
        aperture = '2';
        
    case 10
        current = '38.0';
        aperture = '3';
        
    case 40
        current = '600';
        aperture = '4';
        
    case 60
        current = '1450';
        aperture = '5';
        
    case 120
        current = '6000';
        aperture = '6';       
end

%use the current specified in the input configuration
if ~isempty(config.current)
    current=num2str(config.current);
end

% select dose file
% display('Please choose layer doses...');
% [baseName, folder] = uigetfile('C:\NPGS\Projects\*.txt');
% file_doses = fullfile(folder, baseName);

fprintf('Loading doses...\n');
doses=load(config.doseFile);

% convert dose percentages to actual doses using config.dtc
doses = doses*str2num(config.dtc);

%find layers in use
fprintf('Loading cad file from %s \n',config.cadFile);
cad_t = fileread(config.cadFile);

%Jet-like colors
ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0]};

colorstrings = {};
for i=1:length(ctab)
    ctabmat = ctab{i};
    colorstrings{i}={[sprintf('%03d',ctabmat(1)) ' ' sprintf('%03d',ctabmat(2)) ' ' sprintf('%03d',ctabmat(3))]};
    %colorstrings{i} = {[num2str(ctabmat(1)) ' ' num2str(ctabmat(2)) ' ' num2str(ctabmat(3))]};
    %colorstrings{i} = strrep(colorstrings{i},'0','000');
end

%generate pattern writing text

sdose = [];
scol = {};
for i=1:length(doses)
        sdose(end+1) = doses(i);
        scol{end+1} = colorstrings{i};
end
tot_str_s = '';
nextnum = 2; %layer numbering starts at 2 and goes up with patterns created with urpec
% basic format for each pair of lines is 
% level #, layer #, mode, move, max mag, mag, spacing, spacing, aperture,
% current, color, dwell time, dose, ?
A=str2num(config.spacing{1})*str2num(config.spacing{2})*1e-20*1e4;

strline1 = ['lev_1 ' num2str(nextnum) ' w    0,0    1500    ' config.mag '    ' config.spacing{1} '    ' config.spacing{2} '    ' aperture '     ' current];
tot_str_s = strline1;
for i=1:length(ctab)
    dt=A*sdose(i)*1e-6/(str2num(current)*1e-12)*1e6; %autocompute dwell time.
    cnum=sprintf('%03d',i);
    %strline2 = ['col -' cnum ' ' char(scol{i}) ' 10.5239 ' num2str(sdose(i)) ' 0'];
    strline2 = ['col -' cnum ' ' char(scol{i}) ' ' num2str(dt) ' ' num2str(sdose(i)) ' 0'];

    tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
end


s=tot_str_s;

cd(curDir);


end

