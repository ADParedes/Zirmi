%% Unique to ADP
clear 
dirScript_Z_Supplementary      ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_Z_Supplementary';
dirScript_E_Functions          ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_E_Functions'; 
%% Supplementary_ScriptS3 
%   Zirmi Source Information can be found at:
%   <https://github.com/ADParedes/Zirmi>
%   Written By: Andre Daniel Paredes | email: andre.paredes@ymail.com
%   Matlab Source Information can be found at:
%   <https://www.mathworks.com/help/images/image-analysis.html>
%   Description:  Manual tracing technique to image segment wound region and acquire
%   raw pixel intensity
%% Define the following Directories appropriately
cd          (dirScript_Z_Supplementary)
addpath     (dirScript_E_Functions)
%% Section A: Read Image Data and Preprocess
%-Step A1: Read Image data into workspace from preselected directory
fluorFish                   = imread('fluorescent_fish.tif');
%-Step A2: Format fluorescen_fish image to a standardized bit format to streamline processing
bytes                       = 2^16-1;
fluorFishImg                = im2uint16(fluorFish);
%% Section B: Manually Image Segment Wound Region and Display
disp        ('Outline the Wound Region')
close all; 
imagesc     (fluorFishImg); 
colormap    ('gray'); 
movegui     ('northeast');
title       ('ROS REGION');
RosRegion                   = roipoly();
imshow      (RosRegion)
%% Section C:  Compute Fluorescent Intensity in ROI
areaROI                     = sum(sum(RosRegion));
meanPixelIntensity          = sum(sum(uint16(RosRegion).*fluorFishImg))/areaROI;
display(strcat('Raw Pixel Intensity (abu):',num2str(meanPixelIntensity)))
display(strcat('Region area (pixel^2):',num2str(areaROI)));
