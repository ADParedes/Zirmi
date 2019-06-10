%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update 2017-11-15 - ADP
%Computer hostname changed - previously LLIULM --> ADPAREDES 
%Computer hostname reflects labtop in which CODE was created
% Version 1.3 11-20-17
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% Close all/Clear Variables/Warning
clear; clc; close all;
boo_pc                     =   ispc;

switch      boo_pc
    case {1}
        disp('Zirmi is Compatible with your PC')
    otherwise
        disp('Zirmi is not yet compatible with your MAC...YET.')
        disp('Look for Future updates')
        display                   ('Module 1 DISCONTINUED')
        return
        
end;

%% Warning
disp                    ('Please Note the Following Warnings : ')
warning                 ('Zirmi A (Module 1 Part 1) needs to be ran so as to load (.mat) file To Register Parameters in a Standardized fashion')
pause                     (3)
clc

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
%% Directory Identifications          
switch boo_desktop
    case {1,0}
        if exist('input_16bit')
        else
                input_16bit        ='E';
        %     input_16bit         =input('C for 16bit else 'n' for not? E for external','s');
                if input_16bit=='n' %8-bit (reduced images from Matlab R2013b academic version)
                    disp('Directory Supplementary 1: MatlabR2013b based codes for where you keep data')
                    str_Presumptious            ='C:\Users\AndreDaniel\Documents\SkyDrive\PHD\Prelim\VisualBasic\Z-Stack\Zebrafish';
                elseif input_16bit=='E'  %External Drive
                    disp('Directory Supplementary 2: External Drive E:\ for where you keep data')
                    str_ProcessedData_str       ='F:\AndreDParedes Personal Data\Paredes Data';
                    str_Presumptious            = str_ProcessedData_str;
                    str_Experiment_str          ='E:\'; %Experiment Folder- Starting Point to find Experiment to Analyze
                    str_Functions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\LLIULM\MATLAB\Image Processing\Fun');
                    str_RawData_str             ='F:\AndreDParedes Personal Data\Paredes Data\Processed_MrBig(.mat&.tif))';
                    %---For the Zirmi Function Selection
                    if strcmp(hostname,'ADPAREDES'); %this is the labtop, Formally known as 'LLIULM'
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\Zirmi\Script_E_Functions'); 
                    elseif strcmp(hostname,'ADP-Chicago'); %this is the desktop
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Dropbox\Zirmi\Script_E_Functions');
                    end;   
                    str_Presumptuous            = str_ProcessedData_str;
                    checkConfocal               = 1;
                elseif input_16bit  =='C' 
                    disp('Selected Directory 1: Main Files for where you keep data')
                    str_Functions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\LLIULM\MATLAB\Image Processing\Fun');
%                     str_ExtendedFunctions_str   = strcat('C:\Users\',user,'\OneDrive\m_Scripts\phagosight-master');  
                    str_ProcessedData_str       = strcat('C:\Users\',user,'\OneDrive\MrBigTest\Confocal Data\Processed Data');
                    %---For the Zirmi Function Selection
                    if strcmp(hostname,'ADPAREDES'); %this is the labtop, Formally known as 'LLIULM'
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\Zirmi\Script_E_Functions'); 
                    elseif strcmp(hostname,'ADP-Chicago'); %this is the desktop
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Dropbox\Zirmi\Script_E_Functions');
                    end;   
                    str_Presumptuous            = str_ProcessedData_str;
                    checkConfocal               = 1;
                else
                    disp('Nothing was Selected...genius')
                end;
        end;
    otherwise
        disp('Hello Zirmi User')
        disp('First we Need to Determine Directories')
        str_Presumptuous                = uigetdir('*.mat*', 'Please Select Directory Where you keep Experiment .mat files');
        str_Zfunctions_str              = uigetdir('*', 'Please Select Directory where Zirmi_E_Functions exists');
%         str_ExtendedFunctions_str   = uigetdir('*', 'Please Select Directory where PhagoSight Functions exist');
        
end;
disp        ('END: Directory Selections')
%% Find Folder with Registered input Parameters
cd(str_Presumptuous)
if input_16bit  =='C' 
else
    exp_folder      = uigetdir('*.mat*', 'Please select your uniquely Experiment .mat file');
    cd(exp_folder)
end;
startingFolder      =  str_Presumptuous;
defaultFileName     =  fullfile(startingFolder,'*.*');
% [baseFileName,folder] =uiputfile(defaultFileName,'Just try to specify in a different folder - I dare you') %This is to Overwrite save on a file of choosing.
% if baseFileName == 0
%     % User Clicked the cancel button.
%     return;
% end;
disp                   ('Note: Select the .mat file of the Experiment you wish to analyze')
xFile               =  uigetfile('*.mat*', 'Pick unique Experiment ID MATLAB file with Registration Parameters');
load(xFile);
temp_Xfolder        = xFile(1:5);
xFolder             = strcat(cd,'\',temp_Xfolder);
%-Determined ti_d risky way
posnum              = length(PositionName);
disp                ('Position#/Sampling Frequency Reference Sheet: 7/1; 8/1.08; 9;1.2; 10/1.4'); %ADP specific - not necessarily compatible for other users.
display             ('END: Experimental Registration Parameters')
%% Parameters of Interest
POI.Parameter1              = 2^16-1;   %MaxPixelIntensity - BPP - dependent on Bits 
POI.Parameter2              = 1.64;     %Lateralpixelresolution  - micronsperpixel micronsPerPixel =1.660; %microns per pixel (512x512) BUT image is in 256x256 so reduced by 2   3.102
POI.Parameter3              = 10;       %ZstepMicrons
POI.Parameter4              = ti_d;     %SamplingFrequency
POI.Parameter5              = t_plate;  %Minutes Post Injury of Image Start (MPI)
POI.ParameterA              = 0.70;     %Trackability i.e Minum Cell Track length Thresholdper SM a.k.a MajorityTracksPercent=.70; %<--HIGH% selection - ADP
POI.Parameter_gtA           = 0.9;      %Trackability per GT
POI.ParameterB              = 1;        %StaticLimit (1 provides best result)  i.e. 0.8519um. 
POI.ParameterC              = 65;       %Distance from Wound Margin for Standardizatin of Wound Region 
POI.ParameterS              = 100;      %Leukocyte Spatial Interval (150um)
POI.ParameterZ              = znum;     %Number of Z Positions
%- ADP Specific 
ADP.adp1            = hostname;
ADP.adp2            = user; 
ADP.date            = date;
ADP.boo1            = input_16bit;         % 'C' lets me know this is main confocal directory
ADP.boo2            = boo_desktop;         % labtop or not labtop 
ADP.boo3            = NaN;                 % checkconfocal not done until S1
ADP.boo4            = boo_pc;              % check to see if mac or PC
%-Other Parameters
POI.Parameter10a            = name;                 % str of unique Experiment name (batch imaging-set)
POI.Parameter10b            = xFile;                % metadata file by Experiment (batch imaging-set)
POI.Parameter10c            = xFolder;              % dir of Unique Experiment Processed data
POI.Parameter10d            = NaN;                  % Position String Combined
POI.Parameter10e            = str_Zfunctions_str;   % Zirmi_Script_E_Functions
POI.Parameter11a            = length(PositionName); % val number of positions of this experiment group
POI.Parameter11b            = PositionName;         % str of Position  numbers being analyzed
POI.Parameter11c            = NaN;                  % str of Position  numbers being analyzed done S1
POI.Parameter12             = get(0,'ScreenSize');  % s= [ 1 1 1920 1080]  --> Width-x-(1920) & Height-y-(1080)
POI.Parameter13             = valsPh2;              % This is complex array with all the BF images so I can easy access BF images
%% -GUI for Metadata Inputs. (Diaglog Box)
prompt                  = {'Imaging Input #1 : MaxPixelIntensity:'...
                            ,'Imaging Input #2: LateralPixelResolution (pixel per micron (um))'...:'...
                            ,'Imaging Input #3 : Z Pixel Resolution (um ; micrometer)'...:'...
                            ,'Imaging Input #4: SamplingFrequency(min):'...                            
                            ,'Imaging Input #5: MPI - Image Start - Minutes post injury'...
                            ,'Imaging Input #6: # Z stacks:'...
                            ,'Enter Parameter 9: Spatial Interval Distance S2,S3, and S4(um; micrometer):'...
                            ,'Enter Parameter 5 : Track inclusion detectable Threshold( 1.0 - 0.0)'...
                            ,'Enter Parameter 6: StaticLimit (um ; micrometer)'...'...
                            ,'Enter Parameter 3: Wound Margin Distance for Wound ROI (um ; micrometer)'...
                            };
dlg_title               = 'Input Parameters';
num_lines               = 1;
defaultans              = {num2str(POI.Parameter1),num2str(POI.Parameter2),num2str(POI.Parameter3),num2str(POI.Parameter4),num2str(POI.Parameter5)...
                            ,num2str(POI.ParameterZ),num2str(POI.ParameterS),num2str(POI.ParameterA),num2str(POI.ParameterB),num2str(POI.ParameterC)};
answer                  = inputdlg(prompt,dlg_title,num_lines,defaultans);
PARAMETERS.Parameter1              = str2num(answer{1}); % Max Pixel Intensity:  Bits per pixel (BPP)
PARAMETERS.Parameter2              = str2num(answer{2}); % LateralPixelResolution
PARAMETERS.Parameter3              = str2num(answer{3}); % Z pixel resolution
PARAMETERS.Parameter4              = str2num(answer{4}); % Sampling Frequency
PARAMETERS.Parameter5              = str2num(answer{5}); % MPI
PARAMETERS.ParameterA              = str2num(answer{8}); % Trackability
PARAMETERS.ParameterB              = str2num(answer{9}); % Static Ratio
PARAMETERS.ParameterC              = str2num(answer{10}); % Standardizatin of Wound Region
PARAMETERS.ParameterS              = str2num(answer{7}); % Spatial Intervals from wound
PARAMETERS.ParameterZ              = str2num(answer{6}); % # Z stacks

POI.Parameter1              = PARAMETERS.Parameter1; % MAX PIXEL
POI.Parameter2              = PARAMETERS.Parameter2; % Later Pixel Resolution
POI.Parameter3              = PARAMETERS.Parameter3; % Z pixel resolution
POI.Parameter4              = PARAMETERS.Parameter4; % sampling Frequency
POI.Parameter5              = PARAMETERS.Parameter5; % MPI
POI.ParameterA              = PARAMETERS.ParameterA; % "Trackability"
POI.ParameterB              = PARAMETERS.ParameterB; % Static limit
POI.ParameterC              = PARAMETERS.ParameterC; % Sandardization of wound region
POI.ParameterS              = PARAMETERS.ParameterS; % Space domain S2, S2, S3 distance
POI.ParameterZ              = PARAMETERS.ParameterZ; % Num Z stacks

disp(strcat('Registered Maximum Pixel Intensity:',num2str(PARAMETERS.Parameter1)))
disp(strcat('Registered LateralPixelResolution(um/pixel):',num2str(PARAMETERS.Parameter2)))
disp(strcat('Registered ZstepMicrons(um):',num2str(PARAMETERS.Parameter3)))
disp(strcat('Registered Sampling Frequency:',num2str(PARAMETERS.Parameter4)))
disp(strcat('Registered MPI (minutes) ):',num2str(PARAMETERS.Parameter5)))
disp(strcat('Registered Z stacks:',num2str(PARAMETERS.ParameterZ)))
disp(strcat('Registered wound region distance from wound perimeter:',num2str(PARAMETERS.ParameterB)))
disp(strcat('Registered Track Inclusion Criteria (detectable centriod positions):',num2str(PARAMETERS.ParameterA)))
disp(strcat('Registered StaticLimit:',num2str(PARAMETERS.ParameterC)))
disp(strcat('Registered Spatial Intervals:',num2str(PARAMETERS.ParameterS)))

disp('Saving...')
%% append save pameters to current registered metadata
switch exist('POI') %#ok<EXIST>
    case{1}
        display(strcat('Registered PARAMETERS from : ',...
                        POI.Parameter10a));
    otherwise
        display('No Registered POI or PARAMETER variables')
end;

save                (xFile,'ADP','PARAMETERS','POI','-append')
fullFilename    =   fullfile(str_Presumptuous,xFile);
message         =   sprintf('Your Registration Parameters have been loaded from:\n%s',fullFilename);
fprintf             ('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',fullFilename,fullFilename);
msgbox              (message);

clearvars -except POI PARAMETERS ADP 
display('FINISHED: Zirmi B0_RegisteredMetadata.m ')