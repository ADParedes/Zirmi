%% Loading Keyhole modeled Tracks from Specific Experiment  
%-- Updated 2017-11-16 --ADP 
%-- Updates (e.g. ABO20) & Specific Position (e.g. P2)
% Version 1.3 11-20-17
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% Warning I

%% Identify COMPUTER  (added 2/12/2018)
boo_pc                     =   ispc;
disp                      ('Please Note the Following Warning : ')
warning                   ('Zirmi A needs to be ran on experiment directory before proceeding')
switch      boo_pc
    case {1}
        disp('Zirmi is Compatible with your PC') 
    otherwise
        disp('Zirmi is not yet compatible with your MAC...YET.')
end;
pause                     (2)
%% Warnings II
disp                      ('Please Note the Following Warning : ')
warning                   ('Zirmi A & B0 needs to be ran on experiment directory before proceeding')
pause                     (2)
clc

switch exist('PARAMETERS') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 1: PASS! PARAMETERS exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('CHECK 1:  FAIL --> No Registered  PARAMETER variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run Zirmi B0 before proceeding or directly load relevant .mat file')
        disp                      ('NOTE: You can load a .mat file by directly placing it in  Command Window')
        display                   ('Module 1 Part 2 DISCONTINUED') 
        return
end;
switch exist('POI') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 2: PASS! POI (other Parameters of Interest) exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('CHECK 1:  FAIL --> No Registered  POI variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run Zirmi B0 before proceeding or directly load relevant .mat file')
        disp                      ('NOTE: You can load a .mat file by directly placing it in  Command Window')
        display                   ('Module 1 Part 2 DISCONTINUED') 
        return
end;

%% Computer System Identification
!hostname
hostname                    =   char( getHostName( java.net.InetAddress.getLocalHost ) );
ip                          =   char( getHostAddress( java.net.InetAddress.getLocalHost ) );
user                        =   getenv('UserName');    
if strcmp(hostname,'ADP-Chicago'); %this is the desktop
    disp('Did NOT select Labtop Directory')
    boo_desktop             =   1; %manual save logic of computer type | "1" is Yes for this is desktop
elseif strcmp(hostname,'ADPAREDES') %this is the labtop
    disp('Selected Labtop Directory')
    boo_desktop             =   0;
else
    boo_desktop              =   2;
end;

ADP.adp1            = hostname;
ADP.adp2            = user; 
ADP.date            = date;
ADP.boo1            = 'C';                 % 'C' lets me know this is main confocal directory
ADP.boo2            = boo_desktop;         % labtop or not labtop 
ADP.boo3            = NaN;                 % checkconfocal not done until S1
ADP.boo4            = boo_pc;              % check to see if mac or PC

disp('END: Computer Identified')

%% Ensure Correct Directory
cd(POI.Parameter10c )
%-Selecting Specific Fish (Position e.g P0,P1,P2...) 
%-"P# Directory has Processed Data via PhagoSight
SelectedFolder_str  = uigetdir;
[pa5, name5, ext5]  = fileparts(SelectedFolder_str); %break up file name
%-Move to Directory of Interest
cd(pa5); cd(name5);
len_dir             = length(dir);
len_dir             = len_dir-2;
cd(pa5);
pause()
Parameter11c     = name5;
Parameter15      = 2;%%% % input_v=input('What Frame?')    <----------LOAD FROM THIS FRAME FOR EXAMPLE VISUALIZTION TIFS
Parameter14      = 'y';% input_OorNoO=input('Does it need an O or not  y/n','s');
%% Automatically Detecting .mat files of Interest
%-Paramter14 distinguishes between mat_ORe files and mat_Re files 
%-Paramter14 may vary on PC vs mac.  NOte: Here Written with PC.
if Parameter14=='n'
    mat     ='_mat_';
    Or      = strcat(Parameter11c,mat,'Or/T00001.mat');%'r/T00001.mat');
    Re      = strcat(Parameter11c,mat,'Re/T00001.mat');
    La      = strcat(Parameter11c,mat,'La/T00001.mat');
    Ha      = strcat(Parameter11c,mat,'Ha/handles.mat');%'Ha/handles_ANDRE.mat'
    ha2     = strcat(Parameter11c,mat,'Ha');
elseif Parameter14=='y'
    mat='_mat_O';
    if Parameter15>9
        Or=strcat(Parameter11c,mat,'r/T000',num2str(Parameter15),'.mat');
        Re=strcat(Parameter11c,mat,'Re/T000',num2str(Parameter15),'.mat');
        La=strcat(Parameter11c,mat,'La/T000',num2str(Parameter15),'.mat');
    else
        Or=strcat(Parameter11c,mat,'r/T0000',num2str(Parameter15),'.mat');
        Re=strcat(Parameter11c,mat,'Re/T0000',num2str(Parameter15),'.mat');
        La=strcat(Parameter11c,mat,'La/T0000',num2str(Parameter15),'.mat');            
    end;
    Ha      = strcat(Parameter11c,mat,'Ha/handles.mat');
    ha2     = strcat(Parameter11c,mat,'Ha');
else
    display('neither? Fine then i say YES')
    Parameter14='y';
    return
end;
%% loading .mat files of Interest from PhagoSight
load(Or)
load(Re)
load(La)
%% loading 'handles.mat'
str_FolderNewHandles   = strcat(Parameter11c,mat,'Ha-new');
Folder_newHandles      = exist(str_FolderNewHandles,'dir');
if Folder_newHandles ~= 0
    Ha_new=strcat(Parameter11c,mat,'Ha-new/handles.mat');
    load(Ha_new)
else
    load(Ha)
end;

%% Update Parameters
%- PhagoSight Specific Parameters
ch_Ph2                          = handles.ChannelDistribution(4)-1;
PhagoSight.ch_Ph2               = ch_Ph2;
ch_GFP                          = round(median(handles.ChannelDistribution(1):handles.ChannelDistribution(3)));
PhagoSight.ch_GFP               = ch_GFP;
%- move to example directory of current position in T0001
cd(POI.Parameter10c); cd(Parameter11c); cd('T0002');
%- Reloaded Parameteres
PhagoSight.ParameterZ           = length(dir('*GFP*'));
POI.Parameter11c                = name5;
%-Display Important Parameters
chr1                        = 'Your Registered Parameters have been Updated for Experiment';
chr2                        = POI.Parameter10a;
chr3                        = 'Experiment Positon:';
chr4                        = POI.Parameter11c;
chr5                        = 'SamplingFrequency(frames/minute)';
chr6                        = num2str(POI.Parameter4);
chr7                        = 'Number of Z-Positions';
chr8                        = num2str(PhagoSight.ParameterZ);
message2                    = sprintf([chr1 '\n' chr2 '\n' chr3 '\n' chr4 '\n' chr5 '\n' chr6...
                                    '\n' chr7 '\n' chr8]);
msgbox      (message2);% --> only works for futures -->   set(gcf,'Name','DataIn','Position',[s(3)/1.5 1 s(3)/3  (s(4))/2])
disp        ('END: Updating Parameters')
%% --------NEW SCRIPT----------------------------------------------------------------------------------
display('START2 EXPEDITE:  Performing Wound Region')% http://www.phagosight.org/NF/trackingManual4.php
cd(POI.Parameter10c)
handles.distMaps    = 0;
pause(2)
valsPh2             = POI.Parameter13 ;
Parameter11d        = str2num(POI.Parameter11c(2:end)); 
switch ADP.boo2
    case {0,1}
        POI.Parameter11d    = Parameter11d+1; %ADP uses P0-P11 (so i need to +1)
    otherwise        
end;

%-Need 'valsPh2' ~= valsPH2
if exist('valsPh2') %#ok<EXIST>
    OriginalImg = valsPh2{POI.Parameter11d}{Parameter15}; % 
    checkConfocal=1;
    valsPH2=valsPh2;
else
    OriginalImg = valsPH2{POI.Parameter11d}{Parameter15}; 
    checkConfocal=0;
end;


%%  *WoundGAP *Notochord Spatial Domains %update 2017/11/16
%- Limitation Spatial Domains are not dynamically determined.  The are
%- marked and assumed to be relatively accurate throughout imaging sessions.
%--First Frame USE - I
I_phag          = dataR(:,:,ch_Ph2);
I               = imadjust(OriginalImg);
PH2_8bit        = im2uint8(I);
PH2_re          = imresize(PH2_8bit, [256 256]); 
close all
figure(1)
if handles.rows>300
%     imagesc(I_phag)
    I = OriginalImg;  %Using Frame 2
else
    imshow(PH2_re)
    I = PH2_re;
end;
%---Last Frame USE I2
OriginalImg2    = valsPH2{POI.Parameter11d}{(handles.numFrames-2)}; %Using second to last Frame
optimalzstack   = 4;
I2              = imadjust(OriginalImg2);
PH2_8bit        = im2uint8(I2);
PH2_re          = imresize(PH2_8bit, [256 256]); 
if handles.rows>300
    I2 = OriginalImg2;
    disp('here')
else
    I2=PH2_re;
    disp('there')
end;
%%  USER OUTLINE TOOL
%- *Wound Gap Selection
figure(1);imagesc(I2);
disp('USER OUTLINE TOOL: Please outline WOUND GAP in (double click to finish) ')
wR1             = roipoly();
close all
figure(2);imagesc(I2);
%- *Tip Notochord Selection
disp('USER OUTLINE TOOL: Please outline tip of NOTOCHORD (double click to finish) ')
wR2             = roipoly();
woundRegion     = wR1|wR2;
close all
color_bw_wR     = imoverlay(I2,bwperim(woundRegion),[ 1 0 0]);
figure(3);imagesc(color_bw_wR);set(gcf,'Name','DataIn','Position',[POI.Parameter12(3)/1.5 1 POI.Parameter12(3)/3  (POI.Parameter12(4))/2])
%% Once this wound region has been selected it is passed as a variable
%- to the following functions
handles                 = effectiveDistance(handles,woundRegion);
%- To calculate the maps
handles                 = effectiveTracks(handles,woundRegion);
%- Viualised with the command "mesh"
figure(4);set(gcf,'Name','DataIn','Position',[POI.Parameter12(3)/1.5 1 POI.Parameter12(3)/3  (POI.Parameter12(4))/2])
mesh(handles.distMaps.oriDistMap)
displayActivationPoint(handles,dataR(:,:,optimalzstack));
dataR_WR                = dataR(:,:,ch_GFP).*(1-imdilate(zerocross(woundRegion),ones(3)));
%- Save in Structure Array
PhagoSight.wR1          = wR1;
PhagoSight.wR2          = wR2;
PhagoSight.woundRegion  = woundRegion;
disp('END:Normalizing Spatial Parameters & Wound Region')
%% Reinforce 
close all
ADP.boo3                = checkConfocal;
clearvars -except POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2
%% Save
switch ADP.boo2
    case {0,1}
        disp('Hello ADP - directory selected automatically')
        dir_metadata='C:\Users\AndreDaniel\OneDrive - University of Illinois at Chicago\PhD Work\PhD matlab\Zirmi_metaData';   
    case {2}
        disp('Hello USER: Select a directory you want to save metadata and Export Outcome Measures Database via Excel too')
        dir_metadata = uigetdir();
    otherwise 
        warning('Zirmi Cannot Recognize Your Computer')
end;
switch exist('POI') %#ok<EXIST>
    case{1}
        display(strcat('Registered PARAMETERS from : ',...
                        POI.Parameter10d));
    otherwise
        display('No Registered POI or PARAMETER variables')
end;

ADP.dir_metadat=dir_metadata;
cd(dir_metadata);
disp('Saving...')
POI.Parameter10d        = strcat(POI.Parameter10a,'_',POI.Parameter11c);
save(POI.Parameter10d)
fullFilename    =   fullfile(dir_metadata,POI.Parameter10d);
message         =   sprintf('Your  Parameters have been loaded from:\n%s',fullFilename);
fprintf             ('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',fullFilename,fullFilename);
msgbox              (message);

disp('FINISHED: Zirmi B1-B2_SelectingReformattedPreProcessedData')
 %% NAN 