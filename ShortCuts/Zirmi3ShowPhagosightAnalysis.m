% Shortcut summary goes here
%- Update 2017-11-17
%% Define Variables 
ch_PH2          = PhagoSight.ch_Ph2;
ch_GFP          = PhagoSight.ch_GFP;
s               = POI.Parameter12; %s=get(0,'ScreenSize'); % s= [ 1 1 1920 1080]  --> Width-x-(1920) & Height-y-(1080)
micronsPerPixel = POI.Parameter2; %LateralPixelResolution %old % micronsPerPixel=1.550 %microns per pixel (512x512) BUT image is in 256x256 so reduced by 2   3.102 
framesPerSec    = POI.Parameter4; %ti_d SamplingFrequency
%%- Load Functions
str_oldfunctions='C:\Users\AndreDaniel\Documents\Dropbox\LLIULM\MATLAB\Image Processing\Fun'; %Functions Folder-Upload Functions Andre Created      
addpath(str_oldfunctions);%add path to *.mfile folder and *.tif foldersss;
%%- Figures
% figure(1);set(gcf,'Name','DataIn','Position',[s(3)/2 1 s(3)/2  (s(4))]); imagesc(dataIn(:,:,ch_GFP))
%     figure(2);set(gcf,'Name','DataR'); imagesc(dataR(:,:,ch_GFP))
    figure(3);set(gcf,'Name','dataL'); imshow(dataL(:,:,ch_GFP))


figure(4);set(gcf,'Name','1  = highlights longer (in distance) branches');
        plotTracks(handles,1,:,micronsPerPixel,framesPerSec);
%     figure(4);set(gcf,'Name','1  = highlights longer (in distance) branches');
%     plotTracks(handles,1,1,micronsPerPixel,framesPerSec);
figure(5);set(gcf,'Name','2  = highlights faster branches');%,'Position',[s(3)/2 1 s(3)/2  (s(4)-100)]);
        plotTracks(handles,2,:,micronsPerPixel,framesPerSec);        
%     figure(6);set(gcf,'Name',' 9  = plot ONLY those branches crossing the present Frame');
%         plotTracks(handles,9,:,micronsPerPixel,framesPerSec);
%     figure(7);set(gcf,'Name','10 = with labels (numbers) for the tracks');
%         plotTracks(handles,10,:,micronsPerPixel,framesPerSec); 
%     figure(8);set(gcf,'Name','4  = highlights shorter branches');
%         plotTracks(handles,4,:,micronsPerPixel,framesPerSec);  
%     figure(10);set(gcf,'Name','1  = Velocity','Position',[s(3)/2 1 s(3)/2  (s(4)-200)]);
%         plotTrackStats(handles,10,:,micronsPerPixel,framesPerSec);
% for i=1:19
%     figure();
%     plotTracks(handles,i,1,micronsPerPixel,framesPerSec); %1 is z stack
% end;
%% ******PhagoSight Refernce Sheet **************
%   [-14:-1]    Will be the same options as the positive ones but with X,Y,Z
%   [1:14]      Will plot the tracks with the following options for X,Y,T
%                       1  = highlights longer (in distance) branches
%                       2  = highlights faster branches
%                       3  = highlights longer (in number of frames) branches
%                       4  = highlights shorter branches
%                       5  = highlights slower branches
%                       6  = highlights smaller branches
%                       7  = discard branches with small total distance, i.e. 30% of upper half average distance
%                       8  = discard branches with less than 3 nodes
%                       9  = plot ONLY those branches crossing the present Frame
%                       10 = with labels (numbers) for the tracks
%                       11 = all tracks in green
%                       12 = all tracks in red
%                       13 = top in red, bottom in green
%                       14 = top in green, bottom in red
%                       15 = plot [x,volume,time] as type 1
%                       16 = plot [y,volume,time] as type 1
%                       17 = Merge into plot [x,volume,time] one colour per handle
%                       18 = Merge into plot [y,volume,time] one colour per handle
%                       19 = Merge several separate handles, p