function [ s ] = urpec_alignBlock( config )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%TODO: start new set somewhere

switch config.sm_aperture
    case 30
        sm_aperture_current = '445.0';
        sm_aperture = '1';
        
    case 7.5
        sm_aperture_current = '17.0';
        sm_aperture = '2';
        
    case 10
        sm_aperture_current = '38.0';
        sm_aperture = '3';
        
    case 40
        sm_aperture_current = '600';
        sm_aperture = '4';
        
    case 60
        sm_aperture_current = '1400';
        sm_aperture = '5';
        
    case 120
        sm_aperture_current = '6000';
        sm_aperture = '6';
        
end

switch config.lg_aperture
    case 30
        lg_aperture_current = '445.0';
        lg_aperture = '1';
    case 7.5
        lg_aperture_current = '17.0';
        lg_aperture = '2';
    case 10
        lg_aperture_current = '38.0';
        lg_aperture = '3';
    case 40
        lg_aperture_current = '600';
        lg_aperture = '4';
    case 60
        lg_aperture_current = '1400';
        lg_aperture = '5';
    case 120
        lg_aperture_current = '6000';
        lg_aperture = '6';
end


doses=0.*ones(1,100);

% convert dose percentages to actual doses using config.dtc

% choose .dc2 files
display('Please choose alignment file...');
[cad,dir] = uigetfile('C:\NPGS\Projects\*.dc2');
fullCad = fullfile(dir,cad);
cad=cad(1:end-4);

%find layers in use
cad_t = fileread(fullCad);

scad = strrep(cad_t,' ',''); % remove spaces
ind = regexp(scad,'DoNotUse'); %find index... layers are right after
scad = scad(ind:ind+1000); %shorten
scad(isspace(scad))='x'; %replace blanks with 'x'
indstart = 9;% starting layers index
indend = min(regexp(scad,'x212'))-1; %ending layers index
slstr = scad(indstart:indend); %string of layers used separated with x's
ind_list = regexp(slstr,'xx'); %list of x locations
ind_list = ind_list(2:end); % skip layer 0
slayers = {};
for i = 1:length(ind_list)
    if i<length(ind_list)
        slayers{length(slayers)+1} = slstr(ind_list(i)+2:ind_list(i+1)-1);
    else
        slayers{length(slayers)+1} = slstr(ind_list(end)+2:end-1);
    end
end
%EJC: 4/10/19 added this because was getting issue where more
%things are accidentally getting recognized as layers... this is
%not perfect either though... maybe need to find more robust way to
%pull out layers used
good = 1;
for i = 1:length(slayers)
    if good
        if length(slayers{i})>2
            slayers = slayers(1:i-1);
            good = 0;
        end
    end
end

%colortab from urpec
ctab={[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1] [1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1] [1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1]  };
colorstrings = {};
for i=1:length(ctab)
    ctabmat = 255.*ctab{i};
    colorstrings{i} = {[num2str(ctabmat(1)) ' ' num2str(ctabmat(2)) ' ' num2str(ctabmat(3))]};
    colorstrings{i} = strrep(colorstrings{i},'0','000');
end

%generate pattern writing text
%small
slogic={}; % logical cell array if layer name exists
for i=1:length(doses)
    if any(strcmp(slayers,num2str(i)))
        a = 1;
    else
        a = 0;
    end
    slogic{i} = a;
end
sdose = [];
scol = {};
for i=1:length(slogic)
    if slogic{i}
        sdose(end+1) = doses(i);
        scol{end+1} = colorstrings{i};
    end
end
tot_str_s = '';
nextnum = 2; %layer numbering starts at 2 and goes up with patterns created with urpec
for i=1:4
    strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' w    0,0    825    ' config.write_mag_sm '    528.717    528.717    ' sm_aperture '     ' sm_aperture_current];
    strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
    if i==1
        sm_str = sprintf('lev_%s %s w    0,0    29106    %s    42.2974    42.2974    %s     %s\ncol -001 %s 10.5239 %s 0',...
            slayers{i}, num2str(nextnum), config.write_mag_sm, sm_aperture, sm_aperture_current, char(scol{i}), num2str(sdose(i)));
    else
        %sm_str = sprintf('%s\nlev_%s %s w    0,0    29106    %s    42.2974    42.2974    %s     %s\ncol -001 %s 10.5239 %s 0',...
        %sm_str, slayers{i}, num2str(nextnum), config.write_mag_sm, sm_aperture, sm_aperture_current, char(scol{i}), num2str(sdose(i)));
    end
    nextnum = nextnum + 1;
    if i==1
        tot_str_s = strline1;
    else
        %tot_str_s = [tot_str_s newline strline1];
        tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline1);
    end
    %tot_str_s = [tot_str_s newline strline2];
    tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
end

for i=5:8
    if i==5
        strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' p    0,0    827    ' config.write_mag_sm '    264.359    264.359    ' sm_aperture '     ' sm_aperture_current];
        strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
    else
        strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' w    0,0    827    ' config.write_mag_sm '    264.359    264.359    ' sm_aperture '     ' sm_aperture_current];
        strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
    end
    if i==1
        sm_str = sprintf('lev_%s %s w    0,0    29106    %s    42.2974    42.2974    %s     %s\ncol -001 %s 10.5239 %s 0',...
            slayers{i}, num2str(nextnum), config.write_mag_sm, sm_aperture, sm_aperture_current, char(scol{i}), num2str(sdose(i)));
    else
        %sm_str = sprintf('%s\nlev_%s %s w    0,0    29106    %s    42.2974    42.2974    %s     %s\ncol -001 %s 10.5239 %s 0',...
        %sm_str, slayers{i}, num2str(nextnum), config.write_mag_sm, sm_aperture, sm_aperture_current, char(scol{i}), num2str(sdose(i)));
    end
    nextnum = nextnum + 1;
    if i==1
        tot_str_s = strline1;
    else
        %tot_str_s = [tot_str_s newline strline1];
        tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline1);
    end
    %tot_str_s = [tot_str_s newline strline2];
    tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
end


s=tot_str_s;


end

