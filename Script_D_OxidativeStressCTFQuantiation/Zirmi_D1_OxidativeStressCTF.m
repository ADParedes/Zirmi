%% Calculating ROS Generation at Wound Region 
%2016-10-07 -Updated 2017-09-20- For DICE Indexing
%STEP 1:  Loading matifiles
%STEP 2: Load Analysis (selecting Experiment and Position)
%ONCE THAT IS DONE THEN RUN THIS SCRIPT
%  Version 1.3 09-20-17
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% Contingency #0 :Computer System Identification
clc %DO NOT need to have ran PhagoSight
close all
clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 time Zirmi equationCTF
disp                   ('Please Note the Following Warnings : ')
warning                ('Zirmi A & B0 needs to ran first or script will discontinue') % don't need handles to do this
pause                  (2)
%  editingscript         = 1; %This is for OutputSegout (Automatic)
% editingscript         = 2; %This is for ReferenceSegout (Manual)
%% PC Check
boo_pc                     =   ispc;
switch      boo_pc
    case {1}
        disp('Zirmi is Compatible with your PC')
    otherwise
        disp('Zirmi is not yet compatible with your MAC...YET.')
        return
end;
ADP.boo4            = boo_pc;              % check to see if mac or PC
%% Parameters
PARAMETER_WOUND = 65; 
%% To Run independetly separately from Zirmi B1,& C
user                        =  getenv('UserName');
hostname                    =   char( getHostName( java.net.InetAddress.getLocalHost ) );
switch exist('POI') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('Registered PARAMETERS from : ',...
                                        POI.Parameter10a,'_',POI.Parameter11c)); %Will be empty if NaN
        dir_Xfolder                 = POI.Parameter10c;   % dir of Unique Experiment Processed data
        BPP                         = POI.Parameter1;     % MaxPixelIntensity - BPP - dependent on Bits 
        cd                          (dir_Xfolder);
        name                        = POI.Parameter10a;
        MPI_start                   = PARAMETERS.Parameter5;   % Minutes post injury (Start of Imaging) ; old: t_plate
        SamplingFrequency           = PARAMETERS.Parameter4;   % Sampling Frequency; old: ti_d
        if strcmp(hostname,'ADPAREDES'); %this is the labtop, Formally known as 'LLIULM'
            str_Zirmi_E_Functions           = strcat('C:\Users\',user,'\Documents\Dropbox\Zirmi\Script_E_Functions'); 
        elseif strcmp(hostname,'ADP-Chicago'); %this is the desktop
            str_Zirmi_E_Functions           = strcat('C:\Users\',user,'\Dropbox\Zirmi\Script_E_Functions');
        else
            str_Zirmi_E_Functions       = POI.Parameter10e;
        end; 
        
        boo_pc                      = ADP.boo4;
        str_user                    = ADP.adp2;
        boo_user                    = ADP.boo2;
        
        addpath                     (str_Zirmi_E_Functions); % Allows access to Zirmi_E_Functions
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('No Registered POI or PARAMETER variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run Zirmi B0 before proceeding')
        display                   ('SCRIPT RUN DISCONTINUED')
        return
end;
%--Select the folder again
cd                                  (dir_Xfolder);
disp                                (strcat('Analyzing |Batch:',POI.Parameter10a,'|Unique Fish:',POI.Parameter11c,'|')) 
dir_positionFolder                  = fullfile(dir_Xfolder,POI.Parameter11c);
[dir_Xfolder, str_positionName, ~]  = fileparts(dir_positionFolder); %break up file name % AS DONE IN 'B:Zebra3-LoadAnalysis'
disp                                (strcat('START CTF : ',name,'_',str_positionName))
disp                                ('END: Selected Fish For CTF Analysis') 
%% Determining where to Save
disp            (strcat('Hello Zirmi User :',str_user));
switch boo_user
    case {0,1} %User is ADP
        dir_saveExcel               = strcat('C:\Users\',user,'\OneDrive - University of Illinois at Chicago\PhD Work\PhD matlab\Zirmi_CTF\viaZirmi_xlsx');
        dir_saveMatfiles            = strcat('C:\Users\',user,'\OneDrive - University of Illinois at Chicago\PhD Work\PhD matlab\Zirmi_CTF\viaZirmi_mat');
    otherwise  %User is not ADP
        boo_saveCTF                 =  isfield(ADP,'boo5');
        switch boo_saveCTF
            case {1}
                dir_saveExcel                = ADP.boo5; %Specific to CTF
                dir_saveMatfiles             = ADP.boo6; %Specific to CTF
                
            otherwise               
                disp('Select Directories where files will be saved')
                % Added Update Version 1.4  
                [PATHSTR2,NAME2,EXT2] = fileparts(POI.Parameter10e);
                [PATHSTR,NAME,EXT] = fileparts(PATHSTR2);
                cd(PATHSTR)
                ADP.dir_metadat = strcat(PATHSTR,'\Zirmi Data');
                mkdir(ADP.dir_metadat) 
                
                dir_saveExcel                   = ADP.dir_metadat;
                dir_saveMatfiles                = ADP.dir_metadat;
                ADP.boo5                        = dir_saveExcel;
                ADP.boo6                        = dir_saveMatfiles;
                save(POI.Parameter10b,ADP,'-append')
        end;
end;
disp            ('END: Selected Directory to Save')
%% Calculate Time Parameters from Dir of Interest
%--Access the Position Folder (E.g. P1 or P2..) to determine Frames used
cd              (dir_Xfolder)
cd              (str_positionName)
clear            arr_T tfolder
%--Save the directory and access folder names to determine Frames
d_positionName                                      = dir('*T*'); % Removes anything that is not a T0000 style folder
[arr_T,arr_boo_GT,arr_GT,arr_LoopGT,num_frames ]    = function_calculateGT(d_positionName,SamplingFrequency,MPI_start);
disp            ('END: Time Arrays were Calculated unique to this Fish')
%% Part 1: Perform Quantitation of ROS Baseline (PreInjury) at Region of Interest 
% input_preinjury                             = input('Calculate PreInjry CTF? 1=yes ~=no');
input_preinjury                             = 0; %Assume User did not take a PreInjury
switch input_preinjury
    case {1}
        input_40c                                   ='n'; 
        [ num_PreDf,num_PreBg,num_PreRw,num_PreMF,Img0,Tr0]  = function_preInjuryCTF(dir_Xfolder,str_positionName,BPP,input_40c);
        Prex                                        =[num_PreDf,num_PreRw,num_PreBg];
    otherwise
        warning ('There is no PreInjury Images')
        num_PreDf = NaN;
        num_PreBg = NaN;
        num_PreRw = NaN;
        num_PreMF = NaN;
        Prex      = NaN;
        Img0      = NaN;
        Tr0       = NaN;
end;
%% Part 2: Define/Assign Variables
clc
cd              (dir_positionFolder)
%----------------Strings
dir_positionName= cd;
d_positionName  = dir('*T*');
len_d           = length(d_positionName);
%----------------Single values
arr_time        = 0; %Corresponding time to frame
arr_timepoints  = 0; % selected time points for editing
arr_T_p         = 0; % respective frames correspoinding to 'timepoints'
%----------------Cell Arrays
Arr_Bkgd        ={0};
Arr_Ifish       ={0};%mock overlay
Arr_Raw         ={0};
Arr_Difference  ={0};
I_fish          ={0};%RAW TexRed 
I_rr            ={0};%caudal fin ROI 
I_br            ={0};
I_brs           ={0};%all backgrounds
BWFINAL         ={0};%wound ROI of bwfinal done via mindthegap function 
BWFISH          ={0};
%---------------Matrices
MeanBgs         =(0);
IntDen          =(0);
areaROI         =(0);
CTCF            =(0);
i_z             = 1;
%% Selected Times We use for Manual Tracings
arr_time        = round((arr_T(1:num_frames)*SamplingFrequency+MPI_start)-SamplingFrequency);
default_tP      = [60,90,120]; %This is what is done specifically
switch user
    case{'AndreDaniel'}
        arr_timepoints  = default_tP;
        len_n           = length(arr_timepoints);
    otherwise
        %-Selecting By TimePoint (MPI)
        arr_timepoints  = default_tP;
        len_n           = length(arr_timepoints);
        %-Selecting By Frames
        everyframe      =(1:num_frames);
        every20         =(1:20:num_frames);
end;
% This Loop selected timepoints near arr_timepoints, but does not
% extrapolate data that is not there, therefore if we don't have a 120
% frame or above we will perform a 120 manual check.
tp              =(0);
for j=1:len_n
    val            = arr_timepoints(j); % value to find
    tmp            = abs(arr_time-val);
    [~, idx]       = min(tmp); % idex of closest value
    closest        = arr_time(idx); %closest value
    tp(j)          = idx;
end;
arr_T_p         = cat(2,(1),tp);
iteration       = 200;
everyframe      =(1:num_frames);
every20         =(1:20:num_frames);
everyiteration  =(1:iteration:num_frames);
golditeration   = cat(2,arr_T_p,everyiteration); % (1:5:frames); %This is the way to change everother
%---------------Announcements
frameuser       = golditeration;  %based on frame thus idx
timeuser        = arr_timepoints;
disp            (strcat('Processing...',POI.Parameter10d))
disp            ('The following minutes post injury time points will be checked : ')
disp            (timeuser);
pause           (2)
% disp            ('NOTE: To edit modify "arr_timepoints"  and change accordingly')
%% Part 2: LOOP Quantitate ROS in Region of Injury after Injury
%for i=(arr_LoopGT) %This was the do the 30 minute interations
for i=everyframe; % This is to do every single frame 
    %% Set Parameters
    cd(dir_positionName) %Go back to experiment folder
    dir_time         = d_positionName(i).name; %select Frame folder 
    cd                  (dir_time) %move to this directory
    str_TexRed      = strcat('T0*ex*'); 
    str_GFP         = strcat('T0*GFP*');
    str_BF          = strcat('T0*Ph2*');
    f_GFP           = dir(str_GFP); %select only GFP named files (tiffs assuming nothing else saved in this folder)
    f_TexRed        = dir(str_TexRed); %select only TexRed named files (TIFFS assuming nothing else saved in this folder)
    f_BF            = dir(str_BF);
    arr_Diff        = 0;
    arr_Raw         = 0;
    %% Select TexRed In Focus  
    if length(f_TexRed)>1
        %This section is obsolete because we combine all TexRed Images
        disp('Determining best z_stack');
        for ii=1:length(f_TexRed)
            figure(ii);imagesc(imread(f_TexRed(ii).name))
        end;
        i_z = input('Whats the best z-stack');
    else
        i_z = 1;
    end;
    everyother          = i;       
    seeMeTR             = f_TexRed(i_z).name;
    seeMeBF             = f_BF(i_z).name;
    I_TR                = imread(seeMeTR); 
    I_BF                = imread(seeMeBF);
    I_BF                = imadjust(I_BF);
    %% Qualitative
    %--Variables for If loop
    OriginalImgTR               = I_TR; 
    input_happy                 ='n';
    input_threshold             = 7000;
    if ~isempty(find(i==golditeration))%-User Defined Selected Iterations
        disp                (arr_T_p) %disp        (i)       
        %--Thresholding While Loop -------------------------------THRESHOLD
        while input_happy ~= 'y' %while it DOES NOT equal 'y'
            disp                    ('Current Threshold: Not Accepted')
            close                   all
            A                       = I_TR;
            A(A>=input_threshold)   = 1;
            A(A>1)                  = 0;
            figure();imagesc                 (A);movegui('northeast'); colormap default;
            save_Threshold          = input_threshold;
            input_threshold         = input(strcat('Set the Low Threshold |',...
                                        'Currently Set at:',...
                                        '[',num2str(input_threshold),']',...
                                        '~ to keep:  '));
            if isempty              (input_threshold)
                input_threshold         = save_Threshold;
            else
            end;
            close                   all
            A                       = I_TR;
            A(A>=input_threshold)   = 1;
            A(A>1)                  = 0;
            figure();imagesc                 (A);movegui('northeast'); colormap default;
            input_happy             = input('Are you happy With Threshold y or ~ to  ACCEPT: ','s');
        end;
        disp                    ('Current Threshold: ACCEPTED')
        %--REGION Selection -------------------------------------------ROIs
        if i==1            
            Img1            = I_BF;
            disp            ('Outline the Wound Region of Interest or Entire Caudal Fin; Used for Reference "Mock" Values')
            close all; imagesc(I_BF); colormap('gray'); movegui('northeast');title('ROS REGION');
            RosRegion       = roipoly();
        %--ROI Background----------------------------------------BACKGROUND
            disp            ('Select Background Region')
            boo_bg          = 'y';
            a               = 0;
            temp_Bgs        =(0);
            [x,y]           = size(A);   
            blankimage      = zeros(x,y);
            while (isempty(boo_bg)) || boo_bg=='y' 
                a       = a+1;
                if a==1
                    BgRegions           = blankimage;
                else
                    boo_bg=input('Another Background Region, Zirmi recommends 2-3?  ~ or y = yes, n = break only: ','s');
                    if boo_bg=='n'
                        break
                    elseif isempty(boo_bg)
                    else
                    end;
                end;
                disp                    ('Select Background Region')
                close all; imagesc(I_BF);colormap('gray'); movegui('northeast');title('BACKGROUND REGION');
                BackgroundRegion         = roipoly();
                OriginalImgTR            = I_TR; 
                BgRegions                = BgRegions | BackgroundRegion;

                BackgroundPixelNum       = sum(sum(BackgroundRegion));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/
                BackgroundPixelIntensity = sum(sum(uint16(BackgroundRegion).*OriginalImgTR))/BackgroundPixelNum;
                temp_Bgs(a)=BackgroundPixelIntensity;
                disp(temp_Bgs)
                close all;

            end;
            close all;
        else
            %-None of This is used for Zirmi CTF Quantitation
            disp(strcat('Keeping Previous Background Region on frame: ',num2str(i)))
            close all; 
            temp_outlinebg  = bwperim(I_brs{i-1}); %background ROI's
            temp_outline    = bwperim(I_rr{i-1});
            temp_segout     = imoverlay(I_TR,temp_outline,[1 0 0]);
            temp_segouts    = imoverlay(temp_segout,temp_outline,[1 1 0]);
            temp_segout_BF  = imoverlay(I_BF,temp_outline,[1 0 0]);
            close all;
            figure();imagesc(temp_segout);movegui('southeast');title('RED IMAGE');
            figure();imagesc(temp_segout_BF);colormap('gray'); movegui('northeast');title('BRIGHT FIELD IMAGE');
            rr_input=input(strcat('Do ROS again? y = yes, all else breaks,',...
                                    fprintf('\n'),...
                                    'NOTE: This is not included in ZIRMI CTF Quantitation'),'s');
            if rr_input=='y'
                RosRegion       = roipoly();
            else
            end;             
        end;      
    else
        %-These are automated processing.
        %-Thresholding
            A               = I_TR;
            A(A>=input_threshold) = 1;
            A(A>1)          = 0;
        %- Region Selection
            RosRegion       = I_rr{i-1};
            BackgroundRegion= I_br{i-1};
            BgRegions       = I_brs{i-1};               
            temp_Bgs        = MeanBgs(i-1);
                
    end;
    B                        = A; 
    switch i
        case {1}
            % Zirmi does not use MeanBgs or IntDen variables for CTF
            % calculations. Alternatively, this can ben done more
            % efficently
            MeanBgs(i)                = NaN;
            IntDen (i)                = NaN;
        otherwise
            RawPixelIntensity           = sum(sum(uint16(RosRegionMask).*OriginalImgTR))/CellPixelNum;
            BackgroundPixelNum          = sum(sum(BackgroundRegion));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/
            BackgroundPixelIntensity    = sum(sum(uint16(BackgroundRegion).*OriginalImgTR))/BackgroundPixelNum;               
            % Quantitated cc_rr
            MeanBgs(i)                  = mean(BackgroundPixelIntensity);        
            cc_rr                       = regionprops(RosRegion);
            IntDen(i)                   = cc_rr.Area * RawPixelIntensity; %Integrated Density
    end;
   %% Archived ROI spatial references
    I_fish{i}                = I_TR; %Raw Image
    I_rr{i}                  = RosRegion; %Caudal Fin ROI
    I_br{i}                  = BackgroundRegion; %Background ROI
    I_brs{i}                 = BgRegions; %obsolete
%     disp('End Spatial Referencing ROIs')
    %% Quantitation
    arr_Bkgd  =(0);
    arr_Raw   =(0);
    arr_Diff  =(0);
    for j=1:length(f_TexRed)
        seeMe_TR1           = f_TexRed(j).name;
        seeMeBF1            = f_BF(i_z).name;
        I_TR1               = imread(seeMe_TR1);  %read image to beMe
        I_BF1               = imread(seeMeBF1);
        I_BF1               = imadjust(I_BF1);
        %% Binary of Fish with Threshold
        se90                = strel('line', 3, 90); se0 = strel('line', 3, 0);            
        Bdil                = imdilate(B, [se90 se0]); %figure, imshow(Adil), title('dilated gradient mask');
        BWdfill             = imfill(Bdil, 'holes'); %figure, imshow(BWdfill); title('binary image with filled holes')             
        bw                  = bwareaopen(BWdfill, 5000);%Eliminating objects fewer than 5000 Pixels
        bw2                 = bwmorph(bw,'clean'); %clean not sure, can't hurt
        bw2                 = bwmorph(bw2,'thicken'); %thickens by just a little
        bw                  = bw2;
        bwO                 = bwperim(bw);
        bwO2                 = imdilate(bwO,[se90,se0]);
        bw_color_realsegout = imoverlay(I_TR1,bwO,[1 0 0]);%high is red
        %%%Converse is the Background
        se90bck             = strel('line', 20, 90);
        se0bck              = strel('line', 20, 0);
        bw_background       = imdilate(BWdfill, [se90bck se0bck]);
        bw_background       =~bw_background;
        bw_fish             = xor(bw,bwO); %This is the fish bw    
        %% Automating wound
        %*****These are the variables used for Zirmi CTF Quantification*****        
        if  ~isempty(find(i==golditeration))%i==1
            diameter                    = PARAMETER_WOUND; % approximately 36 um on this confocal setting 
            %Update: Changed I_TR to I_BF - note sometimes I_BF is not easy
            %to see.  Could have a quick algorthm implemented to determine
            %which is best
            [mmp_I bwfinal FinalSegout ]= mindthegap(I_BF1,bw_color_realsegout,bw_fish,diameter);%I_TR
            mfi(i)                      = mmp_I; %ROI mask pixel intensity (mask f intensity)
            BWFINAL{i}                  = bwfinal;
            BWFISH{i}                   = bw_fish;
%             disp ('Zirmi Check Point'); pause()
        else
            BWFINAL{i}                  = BWFINAL{i-1};   
            BWFISH{i}                   = BWFISH{i-1};
        end;
        pixelNum(i)                 = sum(sum(BWFINAL{i}));
        pixelIntensity              = sum(sum(uint16(BWFINAL{i}).*I_TR1))/pixelNum(i);
        mfi(i)                      = pixelIntensity;
        temp_outlinebg              = bwperim(BWFINAL{i});
        temp_segouts                = imoverlay(I_TR,temp_outlinebg,[1 1 0]);
        close all;
        figure();imagesc(temp_segouts);
        %% ROI spatial referencing        
        rR                          = RosRegion&bw; %caudal fin ROI
        rR_outline                  = bwperim(rR); %outline of caudal fin ROI
        rBG_outline                 = bwperim(BackgroundRegion);% %Manual ROI location of background
        segoutIntensity             = I_TR1;
        segoutIntensity(rR_outline) = BPP;
        I_red                       = imadjust(I_TR1, stretchlim(I_TR1), []);
        I_fishRR                    = imoverlay(I_red,rR_outline,[1 0 0]);
        I_fishROIs                  = imoverlay(I_fishRR,rBG_outline,[1 1 0]);
        imshow(I_fishROIs);  
        text( 5, 10, 'Wound Region','Color','red')
        text( 5,30,'Background Region','Color','yellow')
        %% CTCF Calculation
         cc_rr=regionprops(rR);
         temp_area= [cc_rr.Area];
         areaROI(i)= sum(temp_area);
         CTCF(i)   = IntDen(i) - (areaROI(i)*MeanBgs(i));%Corrected Total Cell Fluorescence for Background 
        %%
        %-PreSet
        OriginalImgTR               = I_TR1; % Raw Image
        RosRegionMask               = rR;  %Caudal Fin
        RawPixelIntensity           = 0;
        BackgroundPixelIntensity    = 0;
        CellPixelNum                = sum(sum(RosRegionMask));
        BackgroundPixelNum          = sum(sum(BackgroundRegion));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/
        %-Calculate Pixel Intensity Aveage
        RawPixelIntensity           = sum(sum(uint16(RosRegionMask).*OriginalImgTR))/CellPixelNum;
        BackgroundPixelIntensity    = sum(sum(uint16(BackgroundRegion).*OriginalImgTR))/BackgroundPixelNum;
        DifferenceIntensity         = [RawPixelIntensity-BackgroundPixelIntensity];
        %-Save in arrays
        arr_Bkgd(j)                 = BackgroundPixelIntensity;
        arr_Raw(j)                  = RawPixelIntensity;
        arr_Diff(j)                 = DifferenceIntensity;
        %(sum(sum(OriginalImgBF))-sum(sum(uint16(CellMask).*OriginalImgBF)))/BackgroundPixelNum;
    end;
    %-Finalize Array Savings
    Arr_Ifish{i}             = I_fishROIs; %overaly of segout ROI (1.Bkgd, 2. caudal fin) on Raw image
    Arr_Raw{i}               = arr_Raw; %Discrete Value 
    Arr_Bkgd{i}              = arr_Bkgd;%Discrete Value
    Arr_Difference{i}        = arr_Diff; %Discrete Value
    Bg(i,1)                   = arr_Bkgd(i_z); %Reformated arr version
    Df(i,1)                   = arr_Diff(i_z);%Reformated arr version
    Rw(i,1)                   = arr_Raw(i_z);%Reformated arr version
    disp                        (strcat('Completed: ',num2str((i/num_frames)*100),...
                                    '% [',num2str(i),'] out of [',num2str(num_frames),']'))
%    copyfile(combPh2img.name,combinePh2cd_s*tr)
      
end;
disp('END: ROS Detection Injury Loop')
%% Mock Values ------These values are not used in final calculation by Zirmi to calculate CTF
disp                (strcat('Background : ',num2str(Bg)))
disp                (strcat('Raw        : ',num2str(Rw)))
disp                (strcat('Difference : ',num2str(Df)))
%--------------------Adjust Values
adjBg               = Bg;
idx                 = find(Bg>(std(Bg)+mean(Bg)));
adjBg(idx)          = mean(Bg);
len_Bgnoise         = length(idx);
if len_Bgnoise      > 6
    disp    ('Length of Bg noise higher than SD is more than 5%')
elseif len_Bgnoise  > 12
    disp    ('Length of Bg noise higher than SD is more than 10%')
end;
%-------------------Adjusted for Noise
adjDF               =(0);
adjRw               = Rw;
adjDf               = adjRw-adjBg;
adjPixValues        =[adjDf,adjRw,adjBg];
dBaseline           =((adjDf-num_PreDf)/num_PreDf)*100; %percent change from baseline
disp                ('END: Mock  Values')
%% CTCF WOUND
%-- Search/Find CTF
[f,n,e]                 = fileparts(ADP.dir_metadat);
cd                      (f)
mkdir                   ('Zirmi_CTF')
cd                      ('Zirmi_CTF')
str_fishName            = strcat(name,'_',str_positionName);
str_fishFile            = strcat(str_fishName,'_Wound.mat');
str_dirseaerch          = strcat('*',str_fishFile,'*');
% str_dirseaerch          = strcat('*','AB020P0Wound','*');testusingdirabov
d_uniqueFish            = dir(strcat(str_dirseaerch)); 
boo_uniqueFish          = length(d_uniqueFish);         
wound_IntDen            =(pixelNum).* mfi;  % Integrated Density
arr_time                = round((arr_T(1:num_frames)*SamplingFrequency+MPI_start)-SamplingFrequency); %For Plotting
%% Loading & Determination save(-append,variables)
switch boo_uniqueFish
    case {0} % There is no such file
        input_appendsave                = 1; % YES save over
        boo_CTF                         = 1; % It will exist
    otherwise % There is a file - LOAD IT 
        load(str_fishFile,'-mat')
        clc
        boo_zeroPixelValues             = isfield(Zirmi.CTF,'zeroPixelValues'); % 1 means it exists
        switch boo_zeroPixelValues
            case {1}
                disp                    ('Baseline fish values are registered')
            otherwise
                disp                    ('Baseline fish values are NOT registered')
        end;
        boo_CTF                         = isfield(Zirmi.CTF,'CTF');             % 1 means it exists
        switch boo_CTF
            case {1}
                disp                    ('CTF fish values are registered')
            otherwise
                disp                    ('CTF fish values are NOT registered')
        end;
        input_appendsave                = input('Overwrite Zirmi CTF struct? 1 = yes, 0-9 = n :  ');
end;
%% Saving Zirmi CTF values
switch input_appendsave
    case {1} % Yes overwrite equation CTF      
%         save(str_fishFile,'mmpi0','wound_CTCF','Time');
        % Zirmi Values Used For Quantification
        Zirmi.CTF.CTF_background        = (wound_IntDen - ((pixelNum).*MeanBgs));
        Zirmi.CTF.areaROI               = pixelNum;         % Real Wound ROI Pixel Area (e.g. pix^2) 
        Zirmi.CTF.pixelValues           = mfi;              % Real  mask pixel intensity of MPI
        Zirmi.CTF.intDen                = wound_IntDen;     % Real Wound Integrrated Density 
        Zirmi.CTF.time                  = arr_time;         % MPI time points
        Zirmi.CTF.golditeration         = unique(golditeration);
        Zirmi.CTF.idx                   = arr_T_p;        
%         Zirmi.CTF.zeroPixelValues       = mmpi0;            % In this case uninjured Fish
        % Mock Values not used for Quantification
        Zirmi.Mock.CTF                  = CTCF;             % Mock CTF Using Manual ROI  - Corrected Total Cell Fluorescence for Background        CTCF(i)   = IntDen(i) - (areaROI(i)*MeanBgs(i));%Corrected Total Cell Fluorescence for Background 
        Zirmi.Mock.areaROI              = areaROI;          % MOCK area ROI
        Zirmi.Mock.intDen               = IntDen;           % MOCK 
        Zirmi.Mock.zeroPixelValues      = MeanBgs;          % In This case background
        % Plot values used for plotting and scaling 
        % Remember arr_T is the 1:frames vs arr_time
        Zirmi.Plot.I_preBF_1            = Img0;             % PreInjury BrightField Image   
        Zirmi.Plot.I_preRed_1           = Tr0;              % Preinjury RED IMAGE
        Zirmi.Plot.I_bf_1               = Img1;             % BF at time point 1
        Zirmi.Plot.I_red                = I_fish;
        Zirmi.Plot.bwFINAL              = BWFINAL;          % Mindthegap background subtraction
        Zirmi.Plot.bwFISH               = BWFISH;
        Zirmi.Plot.bwBKGD               = I_brs;            % BWbackground
        Zirmi.Plot.str_x                = name;
        Zirmi.Plot.str_p                = str_positionName;
        Zirmi.Plot.str_f                = str_fishFile; 
        Zirmi.Plot.axisX                = arr_time;
        Zirmi.Plot.frames               = arr_T;
        Zirmi.Plot.sp_frames            = golditeration;    % extra T_p  this is the index of interest
        Zirmi.Plot.idx_frames           = arr_T_p;          % use to be T_p ; Special Frames used to manually select ROI
        Zirmi.Plot.num_frames           = num_frames;       % use to be frames   
        save                            (str_fishFile,'ADP','PARAMETERS','POI','Zirmi') %overwrites - Does not append
        fullFilename                    =   fullfile(dir_saveMatfiles,str_fishFile);
        message                         = sprintf('New file has been saved:\n%s',fullFilename);
        fprintf                         ('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',fullFilename,fullFilename);
        msgbox                          (message);      
    otherwise % Do Not Overwrite, Everything Stays the Same
        disp                ('Zirmi values were not changed (.mat file NOT overwritten)')
end;
%% Plot (conditional)
% Note: Zirmi.CTF.CTF_background = (wound_IntDen - ((pixelNum).*MeanBgs));
switch boo_CTF
       case {1}
           disp('plot CTF');
           close all
           plot(Zirmi.CTF.CTF_background);
       otherwise
           disp ('CTF cannot be calculated without BaselineValues')
end;

disp    ('End: Conditional Plot')
%% Editing Script
lens  = num_frames;
if exist('editingscript')==0

elseif editingscript==1
    OutSeg=(0);
    for K=1:lens
        OutSeg(K)=sum(BWFINAL{K}(:));
    end;
    OutputSegout=BWFINAL;
    disp('END: OutSeg for DICE Index')
elseif editingscript==2
    RefSeg=(0);
    for K=1:lens
        RefSeg(K)=sum(BWFINAL{K}(:));
    end;
    ReferenceSegout=BWFINAL;
    disp('END: RefSeg for DICE Index') 
    %% Common Segout
    CommonSegout={0};
    common      =(0);
    for K=1:lens
        CommonSegout{K} = ReferenceSegout{K} & OutputSegout{K};
        common(K)       = sum(CommonSegout{K}(:));        
    end;
    %% DICE 
    DICE        = 2*(common./(RefSeg+OutSeg));
    DICE.DSC    = DICE';
    DICE.ref    = RefSeg';
    DICE.out    = OutSeg';
    DICE.common = common';
else
end;
%% Fish Region DICE
    FishSegout  ={0};
    outFish     =(0);         
    for K=1:lens
%         CommonSegout{K} = ReferenceSegout{K} & OutputSegout{K};
        outFish(K)       = sum(BWFISH{K}(:));        
    end;
    DICE.out.fish=outFish';
disp('END: TISSUE DICE')
%% Moving Average &&& Create Moving
%check shortcut %
% cd('C:\Users\AndreDaniel\OneDrive\PhD Data\ROS Movie')
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% loops = frames;
% clear F
% F(loops) = struct('cdata',[],'colormap',[]);
% for j = 1:loops
%     figure(j);set(gcf,'name',strcat(name,'-',name5));
%     imshow(Arr_Ifish{j})
%     
%     txt1=strcat(num2str((i+tplate)-(ti_d*(60/60))),' min');
%     text(100,100,txt1);
%     drawnow
%     F(j) = getframe(gcf);
%     close();
% end
%- Playback the movie two times (works but not useful.    
% fig = figure;
% movie(fig,F,2)
%% save as movie 
%- save as avi 
% mkdir(strcat(name,'-',name5))
%     try
%         movie2avi(F,strcat(name,'-',name5),'compression','none');%strcat(name,'-',name5,'VI/video_1.avi'),'compression','none')
%     catch
% %         %an error was detected whilst saving the movie, 
% %         [numRows,numCols,numChannels]= size(F(1).cdata);
% %         for counterFrame = 1:size(F,2)
% %             % iterate over frames and crop if necessary to the first one
% %             F(counterFrame).cdata(numRows+1:end,:,:)=[];
% %             F(counterFrame).cdata(:,numCols+1:end,:)=[];
% %             % to guarantee the same size, add a zero in the last pixel
% %             F(counterFrame).cdata(numRows,numCols,3)=0;
% %         end
% %         movie2avi(F,strcat(name,'-',name5,'VI/video_1.avi'),'compression','none')
%     end
%% save the movie as a GIF
%     [imGif,mapGif] = rgb2ind(F(1).cdata,256,'nodither');
%     numFrames = size(F,2);
% 
%     imGif(1,1,1,numFrames) = 0;
%     for k = 2:numFrames 
%       imGif(:,:,1,k) = rgb2ind(F(k).cdata,mapGif,'nodither');
%     end
%%
%     imwrite(imGif,mapGif,strcat(name,name5),...
%             'DelayTime',0,'LoopCount',inf) %g443800
%% Test indy threshold
% La=bwlabel(bw,4);
% [r,c]=find(La==2); %provides the row and column postions for value "2" 
% 
% % SubarrayIdx
% 
% [B,L] = bwboundaries(bw,'noholes');
% figure(100);set(gcf,'Name','Area of Objectes')
% imshow(label2rgb(L,'copper',[0 0 0]))% Display the label matrix and draw each boundary
% hold on
% %creates * in middle of object
% stats = regionprops(L,'Area','Centroid','PixelIdxList'); %Determines Area of Shapes
% centroids = cat(1, stats.Centroid); 
% plot(centroids(:,1), centroids(:,2), 'y*') 
% hold on
% 
% for k = 1:length(B)
%   boundary = B{k};
%   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', .5)
%     area = stats(k).Area;
%   area_string = sprintf('%2.2f',(area));
%   
%   text(boundary(1,2)-35,boundary(1,1)+13,area_string,'Color','r',...
%        'FontSize',8,'FontWeight','bold');
% end
% 
% %   boundary= B{k}; % obtain (X,Y) boundary coordinates corresponding to label 'k'
%% clear variables
%% clear variables
clearvars -except POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 Zirmi DICE
