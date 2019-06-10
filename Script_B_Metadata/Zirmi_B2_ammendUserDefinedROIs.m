% http://www.phagosight.org/NF/trackingManual4.php
handles.distMaps=0;
pause(2)
%- Update 2017-11-17
%% Define Variables 
s          = POI.Parameter12;
ParameterZ = PhagoSight.ParameterZ;
ch_PH2=handles.ChannelDistribution(4);
ch_GFP=round(median(handles.ChannelDistribution(1):handles.ChannelDistribution(3)))
%%%%%%%%%%%%
    if isempty(ParameterZ)
        disp(name5)
        ParameterZ=input('What position is this (to avoid reference error) ? '); %Update 2017-03-23
        ParameterZ=ParameterZ+1;
    else
    end;
%%%%%%%%%%%
if exist('valsPh2')

    OriginalImg = valsPh2{ParameterZ}{2};
    checkConfocal=1;
    valsPH2=valsPh2;
else
    OriginalImg = valsPH2{ParameterZ}{2};
    checkConfocal=0;
end;
%ch_PH2=4;
I=imadjust(OriginalImg);
PH2_8bit = im2uint8(I);
PH2_re = imresize(PH2_8bit, [256 256]); 
close all
figure(1)
if handles.rows>300
%  imagesc(dataR(:,:,ch_GFP))
    imagesc(dataR(:,:,ch_PH2))
    I=OriginalImg;
else
    imshow(PH2_re)
    I=PH2_re;
end;

% woundRegion  = roipoly();
%%%%%%%%%%%%%update 2016/09/27
disp('This is for epicenter of wound gap -therefore- outline wound gap')
wR1  = roipoly();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

OriginalImg2 = valsPH2{ParameterZ}{(handles.numFrames-2)};
ch_PH2=4;
I2=imadjust(OriginalImg2);
PH2_8bit = im2uint8(I2);
PH2_re = imresize(PH2_8bit, [256 256]); 
close all
figure(2)
if handles.rows>300
%  imagesc(dataR(:,:,ch_GFP))
%     imagesc(dataR(:,:,ch_PH2))
    I2=OriginalImg2;
    disp('here')
else
%     imshow(PH2_re)
    I2=PH2_re;
    disp('there')
end;
imagesc(I2);
disp('This is for epicenter of notochord -therefore- outline tip of notochord')
wR2  = roipoly();

woundRegion=wR1|wR2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end
close all
color_bw_wR=imoverlay(I2,bwperim(woundRegion),[ 1 0 0]);
figure(3);imagesc(color_bw_wR);
%% Once this wound region has been selected it is passed as a variable
%% to the following functions
figure(4);
handles       = effectiveDistance(handles,woundRegion);
% % % % % % % % % % % % % % % % % % % % % % % disp (handles)

%% To calculate the maps
handles      = effectiveTracks(handles,woundRegion);
% % % % % % % % % % % % % % % % % % % % % % % disp (handles)
% % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % disp (handles.distMaps)

%% Viualised with the command "mesh"
 mesh(handles.distMaps.oriDistMap)

% displayActivationPoint(handles,dataR(:,:,2),woundRegion);
%  displayActivationPoint(handles,PH2_re,woundRegion);
% 
% displayActivationPoint(handles,wounRegion);
displayActivationPoint(handles,dataR(:,:,ch_PH2));

dataR_WR = dataR(:,:,ch_GFP).*(1-imdilate(zerocross(woundRegion),ones(3)));
figure(5);set(gcf,'Name','DataIn','Position',[s(3)/1.5 1 s(3)/3  (s(4))/2])
% imagesc(dataR_WR) 
imshow(I2)


