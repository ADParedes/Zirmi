%% Unique to ADP
clear 
dirScript_Z_Supplementary      ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_Z_Supplementary';
dirScript_E_Functions          ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_E_Functions'; 
%% Supplementary_ScriptS1 
%   Zirmi Source Information can be found at:
%   <https://github.com/ADParedes/Zirmi>
%   Written By: Andre Daniel Paredes | email: andre.paredes@ymail.com
%   MATLAB Source Information can be found at:
%   <https://www.mathworks.com/help/images/image-analysis.html>
%   Description:  Single Image Processing Workflow commonly performed through available MATLAB tool. Requires  manual modifications of thresholding and morphological function values per image.
%% Define the following Directories appropriately
cd          (dirScript_Z_Supplementary)
addpath     (dirScript_E_Functions)
%% Section A: Read Image Data and Preprocess
%-Step A1: Read Image data into workspace from preselected directory
fluorFish                   = imread('fluorescent_fish.tif');
%-Step A2: Format fluorescen_fish image to a standardized bit format to streamline processing
bytes                       = 2^16-1;
I                = im2uint16(fluorFish);
%-Step A3: adjust Image to allow for more rigid image segmentations 
background              = imopen(I,strel('disk',15));
I2                      = I - background;
I3                      = imadjust(I2);
figure                  ();imshow(I3);
%% Section B:  Threshold Fish Tail
% Step B1. Detect entire Fish Tail
[~, threshold]      = edge(I, 'sobel');
fudgeFactor         = .55;
BWs                 = edge(I,'sobel', threshold * fudgeFactor);
 figure(), imshow(BWs), title('binary gradient mask');
% Step B2. Morphological Operate Binary Image
se90                = strel('line', 3, 90);
se0                 = strel('line', 3, 0);
BWsdil              = imdilate(BWs, [se90 se0]);
BWdfill             = imfill(BWsdil,'holes');
BW2                 = bwareaopen(BWdfill, 1000);
seD                 = strel('diamond',1);
BWfinal             = imerode(BW2,seD);
BWfinal             = imerode(BWfinal,seD);
BWfinal             = bwareaopen(BWfinal,1000);
%-Step B3.  Show fish tail image segmentation
BWoutline           = bwperim(BWfinal);
thickBWperim        = bwmorph(BWoutline,'thicken',1);
Segout              = I; 
Segout(thickBWperim)= bytes;
imgFinalPerim       = imoverlay(I,thickBWperim,[1 1 0]);   
figure(); 
imshow(imgFinalPerim);
title('outlined original image');
%% Section C:  Extrapolate Binary Image Components for further Analysis
binaryImgComp       = bwconncomp(BWfinal,4);


