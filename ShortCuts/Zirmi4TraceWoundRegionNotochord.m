%% Zirmi B4 (SHORTCUT) - Outline Wound region and notochord
% Update 2019-06-10  Version 1.4
% created by: Andre Daniel Paredes, PhD
% email: andre.paredes@ymail.com
%% Update 2019-06-10
% When ch_Ph2 is selected incorrectly it can error.  
% This portion is customized only for ADP experiment setup.
switch ch_Ph2
    case{8}
    otherwise
        ch_Ph2 = 8;
        PhagoSight.ch_Ph2 = 8;
end;
%% Define variables
Parameter15      = 2;%%% % input_v=input('What Frame?')    <----------LOAD FROM THIS FRAME FOR EXAMPLE VISUALIZTION TIFS
display('START2 EXPEDITE:  Performing Wound Region')% http://www.phagosight.org/NF/trackingManual4.php
cd(POI.Parameter10c)
handles.distMaps    = 0;
pause(2)
valsPh2             = POI.Parameter13 ;
Parameter11d        = str2num(POI.Parameter11c(2:end)); 
switch ADP.boo2
    case {0,1,2}
        POI.Parameter11d    = Parameter11d+1; %ADP uses P0-P11 (so i need to +1)
    otherwise 
        warning('Potential Error')
        disp('Not calibrated for MACs, may not work')
        POI.Parameter11d    = Parameter11d;
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
I=imadjust(OriginalImg);
PH2_8bit = im2uint8(I);
PH2_re = imresize(PH2_8bit, [256 256]); 
close all
figure(1)
if handles.rows>300
    imagesc(dataR(:,:,ch_Ph2))
    I=OriginalImg;  %Using Frame 2
else
    imshow(PH2_re)
    I=PH2_re;
end;
%- *Wound Gap Selection
disp('USER OUTLINE TOOL: Please outline WOUND GAP in (double click to finish) ')
wR1  = roipoly();

OriginalImg2 = valsPH2{POI.Parameter11d}{(handles.numFrames-2)}; %Using second to last Frame
optimalzstack=4;
I2=imadjust(OriginalImg2);
PH2_8bit = im2uint8(I2);
PH2_re = imresize(PH2_8bit, [256 256]); 
close all
figure(2)
if handles.rows>300
    I2=OriginalImg2;
    disp('here')
else
    I2=PH2_re;
    disp('there')
end;
imagesc(I2);
%- *Tip Notochord Selection
disp('USER OUTLINE TOOL: Please outline tip of NOTOCHORD (double click to finish) ')
wR2  = roipoly();
woundRegion=wR1|wR2;
close all
color_bw_wR=imoverlay(I2,bwperim(woundRegion),[ 1 0 0]);
figure(3);imagesc(color_bw_wR);set(gcf,'Name','DataIn','Position',[POI.Parameter12(3)/1.5 1 POI.Parameter12(3)/3  (POI.Parameter12(4))/2])
%- Save in Structure Array
PhagoSight.wR1          = wR1;
PhagoSight.wR2          = wR2;
PhagoSight.woundRegion  = woundRegion;
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

disp('END:Normalizing Spatial Parameters & Wound Region')
%%
pause();
close all
ADP.boo3                        = checkConfocal;
clearvars -except woundRegion POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 hPlotNet
disp('FINISHED: B3_SelectingReformattedPreProcessedData - Wound Region')
