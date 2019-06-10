%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zirmi Module B0 Version 1.4 
% Update: 2019-06-03  
% USER: Andre Daniel Paredes 
% email: andre.paredes@ymail.com
% Description:Script/shortcut creates loads Zirmi Module A metadata and confirms 
% with USER that the metadata is correct. Additionally identifies important
% directories to addpath for function accessibility, loading mat files, 
% and for saving mat files and excel files.
%% Close all/Clear Variables/Warning
clear ; clc; close all;
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
disp('END: Computer Identified')
%% Directory Identifications          
switch boo_desktop
    case {1,0}
        if exist('input_16bit')
        else
            input_16bit = 'C';
                if input_16bit  =='C' 
                    disp('Selected Directory 1: Main Files for where you keep data')
                    %str_Functions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\LLIULM\MATLAB\Image Processing\Fun');
%                     str_ExtendedFunctions_str   = strcat('C:\Users\',user,'\OneDrive\m_Scripts\phagosight-master');  
                    str_ProcessedData_str       = strcat('C:\Users\',user,'\OneDrive - University of Illinois at Chicago\PhD Work\PhD matlab\mpx');
                    %---For the Zirmi Function Selection
                    if strcmp(hostname,'ADPAREDES'); %this is the labtop, Formally known as 'LLIULM'
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\Zirmi\Script_E_Functions'); 
                            str_phagosight_str           = strcat('C:\Users\',user,'OneDrive - University of Illinois at Chicago\PhD Work\PhD matlab\m_Scripts\GoldenPhagosight\phagosight-master'); 

                            
                    elseif strcmp(hostname,'ADP-Chicago'); %this is the desktop
                            str_Zfunctions_str           = strcat('C:\Users\',user,'\Documents\Dropbox\Zirmi\Script_E_Functions');
                    end;   
                    str_Presumptuous            = str_ProcessedData_str;
                    checkConfocal               = 1;
                else
                    disp('Nothing was Selected...genius')
                end;               
        end;
    otherwise
        usb_dirs = {'D','E','F','G','H','J','K','L','M','N'};
        count_usb = 0;
        boo_usb   = 0;
        while boo_usb ==0
            count_usb = 1 + count_usb;
            str_usb = strcat(usb_dirs{count_usb});
            boo_usb = exist(strcat(str_usb,':\AndreDParedes'),'dir');
        end;
        switch str_usb
            case usb_dirs
                input_16bit ='E'; 
                disp(strcat('Directory Supplementary 2: ADP External Drive ',.... 
                                    str_usb,':\ for where you keep data'))
                str_Presumptuous            =strcat(str_usb,...
                    ':\AndreDParedes\Paredes Data\Processed_mpx(.mat&.tif)');
                str_Zfunctions_str          =strcat(str_usb,...
                    ':\AndreDParedes\Paredes Data\Zirmi Phagosight Scripts\Script_E_Functions');
                str_phagosight_str          =strcat(str_usb,...
                    ':\AndreDParedes\Paredes Data\Zirmi Phagosight Scripts\GoldenPhagosight');
            otherwise
                input_16bit ='NotADP';
                disp('Hello Zirmi User')
                disp('First we Need to Determine Directories')
                str_Presumptuous            = uigetdir('*.mat*', 'Please Select Directory Where you keep Experiment .mat files');
                str_Zfunctions_str          = uigetdir('*', 'Please Select Directory "Zirmi_E_Functions" ');
                str_phagosight_str          = uigetdir('*', 'Please Select Directory of Phagosight Functions');
        end;
        addpath(genpath(str_phagosight_str))
end;
disp        ('END: Directory Selections')
%% Find Folder with Registered input Parameters
cd(str_Presumptuous)
switch input_16bit
    case{'C','NotADP','E'}
    otherwise
        disp('ADP you are using directory outside C Drive') 
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
xFile               =  uigetfile('*.*', 'Pick unique Experiment ID MATLAB file (.mat) with Registration Parameters');
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
POI.Parameter_gtA           = 0.8;      %Trackability per GT
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
prompt                  = {'Enter Parameter1: MaxPixelIntensity:'...
                            ,'Enter Parameter2: LateralPixelResolution:'...
                            ,'Enter Parameter3: Z Pixel Resolution (um):'...
                            ,'Enter Parameter4: SamplingFrequency(min):'...
                            ,'Enter Parameter5: MPI - Image Start - Minutes post injury'...
                            ,'Enter ParameterA: Trackability Threshold(%)'...
                            ,'Enter ParameterB: StaticLimit (Pixel)'...
                            ,'Enter ParameterC: Wound Margin Distance for Wound ROI (Pixel)'...
                            ,'Enter ParameterS: Leukocyte S# from Following Distance Interval From wound(um):'...
                            ,'Enter ParameterZ: # Z stacks:'};
dlg_title               = 'Input Parameters';
num_lines               = 1;
defaultans              = {num2str(POI.Parameter1),num2str(POI.Parameter2),num2str(POI.Parameter3),num2str(POI.Parameter4),num2str(POI.Parameter5)...
                            ,num2str(POI.ParameterA),num2str(POI.ParameterB),num2str(POI.ParameterC),num2str(POI.ParameterS),num2str(POI.ParameterZ)};
answer                  = inputdlg(prompt,dlg_title,num_lines,defaultans);
PARAMETERS.Parameter1              = str2num(answer{1}); % Bits per pixel (BPP)
PARAMETERS.Parameter2              = str2num(answer{2}); % LateralPixelResolution
PARAMETERS.Parameter3              = str2num(answer{3}); % Z pixel resolution
PARAMETERS.Parameter4              = str2num(answer{4}); % Sampling Frequency
PARAMETERS.Parameter5              = str2num(answer{5}); % MPI
PARAMETERS.ParameterA              = str2num(answer{6}); % Trackability
PARAMETERS.ParameterB              = str2num(answer{7}); % Static Ratio
PARAMETERS.ParameterC              = str2num(answer{8}); % Static Ratio
PARAMETERS.ParameterS              = str2num(answer{9}); % Spatial Intervals from wound
PARAMETERS.ParameterZ              = str2num(answer{10}); % # Z stacks

POI.Parameter1              = PARAMETERS.Parameter1;
POI.Parameter2              = PARAMETERS.Parameter2;
POI.Parameter3              = PARAMETERS.Parameter3;
POI.Parameter4              = PARAMETERS.Parameter4;
POI.Parameter5              = PARAMETERS.Parameter5;
POI.ParameterA              = PARAMETERS.ParameterA;
POI.ParameterB              = PARAMETERS.ParameterB;
POI.ParameterC              = PARAMETERS.ParameterC;
POI.ParameterS              = PARAMETERS.Parameter1;
POI.ParameterZ              = PARAMETERS.Parameter1;

disp(strcat('Registered Bits per Pixel (BPP):',num2str(PARAMETERS.Parameter1)))
disp(strcat('Registered LateralPixelResolution(um/pixel):',num2str(PARAMETERS.Parameter2)))
disp(strcat('Registered ZstepMicrons(um):',num2str(PARAMETERS.Parameter3)))
disp(strcat('Registered Sampling Frequency:',num2str(PARAMETERS.Parameter4)))
disp(strcat('Registered MPI (minutes) ):',num2str(PARAMETERS.Parameter5)))
disp(strcat('Registered Trackability:',num2str(PARAMETERS.ParameterA)))
disp(strcat('Registered StaticLimit:',num2str(PARAMETERS.ParameterB)))
disp(strcat('Registered Spatial Intervals:',num2str(PARAMETERS.ParameterS)))
disp(strcat('Registered Z stacks:',num2str(PARAMETERS.ParameterZ)))
disp('NOT Saving...')
%% append save pameters to current registered metadata
% save(zfile,'PARAMETERS','POI','-append')
fullFilename    =   fullfile(str_Presumptuous,xFile);
message         =   sprintf('Your Registration Parameters have been loaded from:\n%s',fullFilename);
msgbox(message);
addpath(genpath(POI.Parameter10e),genpath(str_phagosight_str)); % ZIRMI E Functions directory
clearvars -except POI PARAMETERS ADP
display('FINISHED: S0_RegisteredMetadata.m ')