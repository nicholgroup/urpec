function [] = urpec_makeRunFile(entities)
% [] = urpec_makeRunFile(entities)
%   Makes a run file for NPGS.
%   See the run_urpec script for how to set up the entities argument.
%
%   This function assumes you have converted all pattern files to .dc2
%   format and have saved them for NPGS. Note that you will also need to
%   resave all .dc2 files for NPGS once you have made the run file with
%   this function. Resaving all pattern files will cause NPGS to reload all
%   design files and make sure the run file is ok.


lineChar='\r\n'; %Needed for all new lines

head=struct();

head.header.str=['yx01,0100nnyy0200,1xy',lineChar,'%d',lineChar];
head.header.val=length(entities);

head.magComment.str=['Comment T',lineChar,...
    '1',lineChar,...
    '<MAG>',lineChar,...
    '2',lineChar,...
    'Set mag to %dx on SEM',lineChar,...
    'Are you in NPGS mode?',lineChar]
head.magComment.val=120;

head.waitComment.str=['Comment T',lineChar,...
    '1',lineChar,...
    '<WAIT>',lineChar,...
    '1',lineChar,...
    'Wait to let the stage settle.',lineChar];

head.align.str=['Comment T',lineChar,...
    '1',lineChar,...
    '<WAIT>',lineChar,...
    '1',lineChar,...
    'Wait to let the stage settle.',lineChar];

head.align.str=['%s    %s',lineChar,...
    '1',lineChar,...
    '0,0',lineChar];
head.align.file='align';
head.align.mode='AL'; %A1 for auto1, AL for manual

head.write.str=['%s',lineChar,...
    '1',lineChar,...
    '0,0',lineChar];
head.write.file='saw_gaas30kv200_1';

head.move.str=['MoveOnly M',lineChar,...
    '%d,%d',lineChar];
head.move.val1=400;
head.move.val2=400;

rf=[];
writeBlock=[];

for i=1:length(entities)
    rf=[rf sprintf(head.(entities(i).type).str,entities(i).val{:})];
    
    if strmatch(entities(i).type,'write')
        writeBlock=[writeBlock urpec_writeBlock(entities(i)),lineChar,'*',lineChar];
    end
    
    if strmatch(entities(i).type,'align')
        writeBlock=[writeBlock urpec_alignBlock(entities(i)),lineChar,'*',lineChar];      
    end
end

runFile=[rf writeBlock];


curDir=pwd;
cd(entities(1).dir);
f = input('Please enter run file filename without a file extension (example: DD_L2_SiGe).' ,'s');
RunFile_Name = strcat(f,'.RF6');
fid = fopen(RunFile_Name,'w');
fprintf(fid,runFile);
fclose(fid);

cd(curDir);


end

