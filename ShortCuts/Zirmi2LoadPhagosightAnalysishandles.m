%Zirmi2LoadPhagosightAnalysishandles
%% Loading Keyhole modeled Tracks from Specific Experiment  
%-- Updated 2017-11-16 --ADP 
%-- Updates (e.g. ABO20) & Specific Position (e.g. P2)
disp('START1: Only RUN after S1 script & SA script')
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
%%-Computer System Identification
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
cd(pa5)

Parameter11c     = name5;
Parameter13      = 2;%%% % input_v=input('What Frame?')    <----------LOAD FROM THIS FRAME FOR EXAMPLE VISUALIZTION TIFS
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
    if Parameter13>9
        Or=strcat(Parameter11c,mat,'r/T000',num2str(Parameter13),'.mat');
        Re=strcat(Parameter11c,mat,'Re/T000',num2str(Parameter13),'.mat');
        La=strcat(Parameter11c,mat,'La/T000',num2str(Parameter13),'.mat');
    else
        Or=strcat(Parameter11c,mat,'r/T0000',num2str(Parameter13),'.mat');
        Re=strcat(Parameter11c,mat,'Re/T0000',num2str(Parameter13),'.mat');
        La=strcat(Parameter11c,mat,'La/T0000',num2str(Parameter13),'.mat');            
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
ch_Ph2                          = handles.ChannelDistribution(4);
PhagoSight.ch_Ph2               = ch_Ph2;
ch_GFP                          = round(median(handles.ChannelDistribution(1):handles.ChannelDistribution(3)));
PhagoSight.ch_GFP               = ch_GFP;
%- move to example directory of current position in T0001
cd(POI.Parameter10c); cd(Parameter11c); cd('T0002');
%- Reloaded Parameteres
PhagoSight.ParameterZ                  = length(dir('*GFP*'));
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
message2                = sprintf([chr1 '\n' chr2 '\n' chr3 '\n' chr4 '\n' chr5 '\n' chr6...
                                    '\n' chr7 '\n' chr8]);
msgbox(message2);% --> only works for futures -->   set(gcf,'Name','DataIn','Position',[s(3)/1.5 1 s(3)/3  (s(4))/2])
clearvars -except POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2
disp('END: Updating Parameters')
disp('NOT SAVING')
POI.Parameter10d            = strcat(POI.Parameter10a,'_',POI.Parameter11c);
disp('FINISHED')