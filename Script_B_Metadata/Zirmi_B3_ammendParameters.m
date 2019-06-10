%% Default Parameters
POI.Parameter1              = 2^16-1;   %MaxPixelIntensity - BPP - dependent on Bits 
POI.Parameter2              = 1.64;     %Lateralpixelresolution  - micronsperpixel micronsPerPixel =1.660; %microns per pixel (512x512) BUT image is in 256x256 so reduced by 2   3.102
POI.Parameter3              = 10;       %ZstepMicrons
POI.Parameter4              = ti_d;     %SamplingFrequency
POI.Parameter5              = t_plate;  %Minutes Post Injury of Image Start (MPI)
POI.ParameterA              = 0.70;     %Trackability i.e Minum Cell Track length Thresholdper SM a.k.a MajorityTracksPercent=.70; %<--HIGH% selection - ADP
POI.Parameter_gtA           = 0.9;      %Trackability per GT
POI.ParameterB              = 1;        %StaticLimit (1 provides best result)  i.e. 0.8519um. 
POI.ParameterC              = 65;       %Distance from Wound Margin for Standardizatin of Wound Region 
POI.ParameterS              = 150;      %Leukocyte Spatial Interval (150um)
POI.ParameterZ              = znum;     %Number of Z Positions
%-Other Parameters
POI.Parameter10a            = name;                 % str of unique Experiment name (batch imaging-set)
POI.Parameter10b            = zfile;                % metadata file by Experiment (batch imaging-set)
POI.Parameter10c            = zfolder;              % dir of Unique Experiment Processed data
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
                            ,'Enter ParameterS: Leukocyte S# from Following Distance Interval From wound(um):'...
                            ,'Enter ParameterZ: # Z stacks:'};
dlg_title               = 'Input Parameters';
num_lines               = 1;
defaultans              = {num2str(POI.Parameter1),num2str(POI.Parameter2),num2str(POI.Parameter3),num2str(POI.Parameter4),num2str(POI.Parameter5)...
                            ,num2str(POI.ParameterA),num2str(POI.ParameterB),num2str(POI.ParameterS),num2str(POI.ParameterZ)};
answer                  = inputdlg(prompt,dlg_title,num_lines,defaultans);
PARAMETERS.Parameter1              = str2num(answer{1}); % Bits per pixel (BPP)
PARAMETERS.Parameter2              = str2num(answer{2}); % LateralPixelResolution
PARAMETERS.Parameter3              = str2num(answer{3}); % Z pixel resolution
PARAMETERS.Parameter4              = str2num(answer{4}); % Sampling Frequency
PARAMETERS.Parameter5              = str2num(answer{5}); % MPI
PARAMETERS.ParameterA              = str2num(answer{6}); % Trackability
PARAMETERS.ParameterB              = str2num(answer{7}); % Static Ratio
PARAMETERS.ParameterS              = str2num(answer{8}); % Spatial Intervals from wound
PARAMETERS.ParameterZ              = str2num(answer{9}); % # Z stacks
disp(strcat('Registered Bits per Pixel (BPP):',num2str(PARAMETERS.Parameter1)))
disp(strcat('Registered LateralPixelResolution(um/pixel):',num2str(PARAMETERS.Parameter2)))
disp(strcat('Registered ZstepMicrons(um):',num2str(PARAMETERS.Parameter3)))
disp(strcat('Registered Sampling Frequency:',num2str(PARAMETERS.Parameter4)))
disp(strcat('Registered MPI (minutes) ):',num2str(PARAMETERS.Parameter5)))
disp(strcat('Registered Trackability:',num2str(PARAMETERS.ParameterA)))
disp(strcat('Registered StaticLimit:',num2str(PARAMETERS.ParameterB)))
disp(strcat('Registered Spatial Intervals:',num2str(PARAMETERS.ParameterS)))
disp(strcat('Registered Z stacks:',num2str(PARAMETERS.ParameterZ)))