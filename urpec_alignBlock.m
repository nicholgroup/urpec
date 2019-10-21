function [ s ] = urpec_alignBlock( config )
%[ s ] = urpec_alignBlock( config )
%   Generates text for an NPGS run file for alignment.
%   See run_urpec for an example on how to set this up.

%This function assumes that each alignment set consists of four marks, and
%that the layers are aranged in order. For example, layers 1-4 correspond
%to the coarse alignment, layers 5-8 correspond to the next alignment, etc.
%This function only handles 12 layers now, but it could easily be extended.

curDir=pwd;
cd(config(1).dir);

switch config.aperture
    case 30
        current = '445.0';
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
        current = '1400';
        aperture = '5';
        
    case 120
        current = '6000';
        aperture = '6';
        
end

doses=0.*ones(1,100);


cad_t = fileread(config.cadFile);

scad = strrep(cad_t,' ',''); % remove spaces
ind = regexp(scad,'DoNotUse'); %find index... layers are right after
scad = scad(ind:min([ind+1000,length(scad)])); %shorten
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
    strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' w    0,0    825    ' config.mag '    ' config.spacing{1} '    ' config.spacing{2} '    ' aperture '     ' current];
    strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
    nextnum = nextnum + 1;
    if i==1
        tot_str_s = strline1;
    else
        tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline1);
    end
    tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
end

try
    for i=5:8
        if i==5
            strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' p    0,0    827    ' config.mag '    ' config.spacing{3} '    ' config.spacing{4} '    ' aperture '     ' current];
            strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
        else
            strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' w    0,0    827    ' config.mag '    ' config.spacing{3} '    ' config.spacing{4} '    ' aperture '     ' current];
            strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
        end
        nextnum = nextnum + 1;
        if i==1
            tot_str_s = strline1;
        else
            tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline1);
        end
        tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
    end
catch
end


try
    for i=9:12
        if i==9
            strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' p    0,0    827    ' config.mag '    ' config.spacing{5} '    ' config.spacing{6} '    ' aperture '     ' current];
            strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
        else
            strline1 = ['lev_' slayers{i} ' ' num2str(nextnum) ' w    0,0    827    ' config.mag '    ' config.spacing{5} '    ' config.spacing{6} '    ' aperture '     ' current];
            strline2 = ['col -001 ' char(scol{i}) ' 2 ' num2str(sdose(i)) ' 1'];
        end
        nextnum = nextnum + 1;
        if i==1
            tot_str_s = strline1;
        else
            tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline1);
        end
        tot_str_s=sprintf('%s\r\n%s',tot_str_s,strline2);
    end
catch
end

s=tot_str_s;

cd(curDir);
end

