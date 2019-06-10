%% Unique to ADP
clear 
dirScript_Z_Supplementary      ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_Z_Supplementary';
dirScript_E_Functions          ='C:\Users\ilyaVM\Dropbox\Zirmi\Script_E_Functions'; 
%% Supplementary_ScriptS2 
%   Zirmi Source Information can be found at:
%   <https://github.com/ADParedes/Zirmi>
%   Written By: Andre Daniel Paredes | email: andre.paredes@ymail.com
%   Description: This uses “mindthegap” function incorporated into Zirmi as 
%   a means to reproduce a wound region, acquire raw pixel intensity values, and eliminate 
%   confounding fluorescence within the wound gap that is difficult and time consuming
%   to perform manually.
%% Define the following Directories appropriately
cd      (dirScript_Z_Supplementary)
addpath (dirScript_E_Functions)
%% Section A:  Read Image Data and Preprocess
%-Step A1: Read Image data into workspace from preselected directory
fluorFish                   = imread('fluorescent_fish.tif');
%-Step A2: Format fluorescen_fish image to a standardized bit format to streamline processing
bytes                       = 2^16-1;
fluorFishImg                = im2uint16(fluorFish);
%-Step A3: adjust Image to allow for more rigid image segmentations 
image_background            = imopen(fluorFishImg,strel('disk',15));
fluorFishImg_v2             = fluorFishImg - image_background;
fluorFishImg_v3             = imadjust(fluorFishImg_v2);
%% Section B:  Threshold Fish Tail
%-Step B1: Threshold Fish Fluoresence to ascertain binary image
[~, threshold_fluor]        = edge(fluorFishImg, 'sobel');
fudgeFactor                 = .55;
binaryFishImg               = edge(fluorFishImg,'sobel', threshold_fluor * fudgeFactor);
%-Step B2. Apply Morphological Operations to enhance image segmentation
structElement_90            = strel('line', 3, 90);
structElement_0             = strel('line', 3, 0);
binaryFishImg_dil           = imdilate(binaryFishImg, [structElement_90 structElement_0]); 
binaryFishImg_fill          = imfill(binaryFishImg_dil,'holes');
binaryFishImg_open          = bwareaopen(binaryFishImg_fill, 1000);
structElement_D             = strel('diamond',1); 
binaryFishImg_erode_v1      = imerode(binaryFishImg_open,structElement_D);
binaryFishImg_erode_v2      = imerode(binaryFishImg_erode_v1,structElement_D);
binaryFishImg_Final         = bwareaopen(binaryFishImg_erode_v2,1000);
%-Step B3.  Show fish tail image segmentation
fishImg_perim               = bwperim(binaryFishImg_Final);
fishImg_perim_v2            = bwmorph(fishImg_perim,'thicken',1);
Segout                      = fluorFishImg; 
Segout(fishImg_perim_v2)    = bytes;
disp_fish_image_seg         = imoverlay(fluorFishImg,fishImg_perim_v2,[1 1 0]);   
figure1                     = figure(); 
imshow(disp_fish_image_seg); 
title ('flourescent fish tail image segmentation');
%% Section C:  Standardize image segmentation of "V-Shaped" injury 
%-Step C1: Determine Image Size
[width_x,height_y]          = size(fluorFishImg);
%-Step C2: Trace Wound Gap
disp    ('Please trace a "V" along te wound perimeter')
warning ('Be sure to extend "V" line to open water')
warning ('DO NOT close the V or make any closed shape')
f2                          = figure(2);
imagesc (disp_fish_image_seg); 
title   ('Trace Wound Perimeter');
%-Step C3: Identify and check wound gap region in reference to image
[nodes_xPos, nodes_yPos]                      = getline(f2);
binary_woundgap_seg         = poly2mask(nodes_xPos,nodes_yPos,width_x,height_y);
exist_woundgap              = sum(sum(binary_woundgap_seg));
if exist_woundgap           >  0
    exist_woundgap              = 1;
    bw_woundgap_outline         = bwperim(binary_woundgap_seg);
    Segout_woundgap             = imoverlay(fluorFishImg,...
                                    binary_woundgap_seg,[1 0 0]);
    f3                          = figure(3);    
    imshow(Segout_woundgap);
    title   ('Trace Wound Perimeter');
else
    warning ('Wound Perimeter was not Properly traced')
    disp    ('Please rerun script')
    disp    ('Supplementary_ScriptS2 Terminated')
    return
end;
%% Section D:  Automate wound Region ROI image segmentation 
clc; close all
%-Step D1:  Predefine variables
%           User-determiend radial Distance from wound Gap
roi_rad_dist                = 15; % user based
%           Arrays  used to store computations
Arr_posCirCent              ={0};
Arr_xCirCent                ={0};
Arr_yCircCent               ={0};
Binary_circMASK                       ={0}; %  = combo;
Binary_perimROI                      ={0};%  = raw_combo;
point_dists                 =(0);
%-Step D2-D6:  Employ for loop to address each node             
lensNode                    = length(nodes_xPos);
countTraceNodes           = 0;
for countTraceNodes       = 1:(lensNode-1)
%-Step D2: Calculate Line equation between traced nodes
    temp_cirCenter                  = 0; % clear variable
    %       First or iniial node position
    temp_xPos_i                     = nodes_xPos(countTraceNodes); 
    temp_yPos_i                     = nodes_yPos(countTraceNodes);
    %       Second or final node position 
    if countTraceNodes == lensNode        
        temp_xPos_f                 = nodes_xPos(1);
        temp_yPos_f                 = nodes_yPos(1);
    else
        temp_xPos_f                 = nodes_xPos(countTraceNodes+1);
        temp_yPos_f                 = nodes_yPos(countTraceNodes+1);
    end;
    %       Display trace line between the two nodes 
    f2                              = figure(2+countTraceNodes);
    imagesc (fluorFishImg); 
    hold    on
    plot    (nodes_xPos,nodes_yPos,'LineWidth',2,'Color','r')
    plot    (temp_xPos_i,temp_yPos_i,'c*','LineWidth',2)
    plot    (temp_xPos_f,temp_yPos_f,'y*','LineWidth',2)
    hold    on
    %       Calculate Euclid Distance Between Nodes
    point_dists(countTraceNodes)    = pdist([temp_xPos_i,...
                                        temp_yPos_i;temp_xPos_f,...
                                        temp_yPos_f]);
    temp_point_dist                 = point_dists(countTraceNodes);
    %       Number of circle centers for every pixel between nodes     
    numCircles                      = round(temp_point_dist); 
    %       Calculate line equation in reference to image
    lineCoeffs                      = polyfit([temp_xPos_i, temp_xPos_f],...
                                        [temp_yPos_i, temp_yPos_f],1);
    aCoeff                          = lineCoeffs (1);
    bCoeff                          = lineCoeffs (2);
    min_xPos    = min(nodes_xPos(countTraceNodes:countTraceNodes+1));
    max_xPos    = max(nodes_xPos(countTraceNodes:countTraceNodes+1));
    min_yPos    = min(nodes_yPos(countTraceNodes:countTraceNodes+1));
    max_yPos    = max(nodes_yPos(countTraceNodes:countTraceNodes+1));
%-Step D3: Map circle centers for wound region ROI   
    %      Determine circle center positions respective to image
    cirCent_xPos                    = linspace(min_xPos,max_xPos,...
                                        numCircles); 
    cirCent_yPos                    = aCoeff*cirCent_xPos+bCoeff;  
    temp_cirCenter                  =[cirCent_xPos',cirCent_yPos']; 
    Arr_posCirCent{countTraceNodes} = temp_cirCenter;
    Arr_xCirCent{countTraceNodes}   = cirCent_xPos;
    Arr_yCircCent{countTraceNodes}  = cirCent_yPos;
    %     Define Variables for sub for loop  
    arr_combomask                   ={0};
    binary_indyCirImages            ={0};
    mean_maskPixInten               =(0);   
%-Step D4-D6: Sub for loop defining finite circles comprising ROI
    for countCirc                   = 1:numCircles
%-Step D4: Create circle boundaries at center position
        indy_xCirCent                   = cirCent_xPos(countCirc);
        indy_yCirCent                   = cirCent_yPos(countCirc);
        THETA                           = linspace(0, 2 * pi, 1000);
        RHO                             = ones(1, 1000) * roi_rad_dist;
        [xBounds,yBounds]               = pol2cart(THETA, RHO);
        pos_xBounds                     = xBounds + indy_xCirCent;
        pos_yBounds                     = yBounds + indy_yCirCent;     
        plot(pos_xBounds, pos_yBounds, '-',...
                'linewidth',2,'color','g');
%-Step D5: Create binary mask of circle ROI in fish tissue
        temp_xPos                       = indy_xCirCent  ;   
        temp_yPos                       = indy_yCirCent;     
        [xMesh,yMesh]                   = meshgrid(-(temp_xPos-1):(width_x-temp_xPos),...
                                            -(temp_yPos-1):(height_y-temp_yPos));                                        
        binary_circMask                 =((xMesh.^2+yMesh.^2)<=roi_rad_dist^2); 
        %          Corrected circle mask to be inside fish tissue
        binary_corrMask        = binary_circMask & binaryFishImg_Final; 
%-Step D6: Calculate mean mask pixel intensity of combined circle masks
        maskPixNum                      = sum(sum(binary_corrMask));        
        maskPixInten                    = sum(sum(uint16(binary_corrMask).*fluorFishImg))/maskPixNum;
        mean_maskPixInten(countCirc)    = maskPixInten;
        binary_indyCirImages{countCirc} = binary_corrMask;
        %  Combine Individual Circle Masks in fish issue to complete ROI
        if countCirc==1
            binary_circComb     = binary_circMask;
            binary_corrCombo    = binary_corrMask;
        else
            binary_circComb     = binary_circComb + binary_circMask;
            binary_circComb     = im2bw(binary_circComb);
            binary_corrCombo    = binary_corrCombo + binary_corrMask;
            binary_corrCombo    = im2bw(binary_corrCombo);
        end;
        arr_combomask{countCirc}= binary_corrCombo;
        hold on
    end;
    Binary_circMASK{countTraceNodes}   = binary_corrCombo;
    Binary_perimROI{countTraceNodes}   = binary_circComb;
    
    
end;
%% Section E:  Display and Check Final wound perimeter Image Segmentation
clc; close all;
%-Step E1:  Combine ROI masks
dispImg                   = fluorFishImg;
blankImg                  = binaryFishImg_Final;
for countROIs             = 1:countTraceNodes      
     if countROIs         == 1
         woundPerimMask         = Binary_circMASK{countROIs};
         circRawMask            = Binary_perimROI{countROIs};
     else
         woundPerimMask         = woundPerimMask + Binary_circMASK{countROIs};
         woundPerimMask         = im2bw(woundPerimMask);
         circRawMask            = circRawMask + Binary_perimROI{countROIs};
         circRawMask            = im2bw(circRawMask);
     end;
     blankImg             = imoverlay(blankImg,...
                                bwperim(Binary_circMASK{countROIs}),...
                                0.3.*rand(1,3));
     dispImg              = imoverlay(dispImg,...
                                bwperim(Binary_perimROI{countROIs}),...
                                0.6.*rand(1,3));
 end;
%-Step E2:  Display ROI masks
circImgPerim        = imoverlay(fluorFishImg,bwperim(circRawMask),[1 1 1]);
finalImgPerim       = imoverlay(fluorFishImg,bwperim(woundPerimMask),[1 1 0]); 
figure(); 
imshow(circImgPerim); 
title('Complete Circles Final Image'); 
hold on  
plot(nodes_xPos,nodes_yPos,'LineWidth',2,'Color','r')   
hold off
figure(); imshow(finalImgPerim); title('Wound Region Final Image'); hold on  
plot(nodes_xPos,nodes_yPos,'LineWidth',2,'Color','r')   
hold off
%% Section F: Display and Check Final wound ROI Image Segmentation  
%-Step F1: Remove Confounding fluorescent debris between wound gap
if exist_woundgap       >0
    exist_woundgap          = 1;
    binaryFinal             = ~woundPerimMask | binary_woundgap_seg;
    binaryFinal             = ~binaryFinal;
else
end; 
%-Step F2: Display final ROI Image Segmentation
FinalSegout             = imoverlay(binaryFishImg_Final,...
                            bwperim(binaryFinal),[1 0 1]);
figure  (); 
imshow  (FinalSegout); 
title   ('Fish Mask'); 
hold    on  
plot    (nodes_xPos,nodes_yPos,'LineWidth',2,'Color','r')   
hold    off
figure(); 
imshow  (binaryFinal); 
title   ('only regions of interest'); 
hold    on
%% Section G:  Compute Fluorescent Intensity in ROI
areaROI                 = sum(sum(binaryFinal));
meanPixelIntensity      = sum(sum(uint16(binaryFinal).*fluorFishImg))/areaROI;
