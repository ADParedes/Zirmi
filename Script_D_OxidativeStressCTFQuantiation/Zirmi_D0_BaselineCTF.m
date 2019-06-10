%% Calculating ROS Generation at Uninjured Caudal Fin (BASELINE)
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
disp                                (strcat('Analyzing |Batch:',POI.Parameter10,'|Unique Fish:',POI.Parameter11c,'|')) 
dir_positionFolder                  = fullfile(dir_Xfolder,POI.Parameter11c);
[dir_Xfolder, str_positionName, ~]  = fileparts(dir_positionFolder); %break up file name % AS DONE IN 'B:Zebra3-LoadAnalysis'
disp                                (strcat('START CTF : ',name,'_',str_positionName))
disp                                ('END: Selected Fish For CTF Analysis') 
%% Determining where to Save
disp            (strcat('Hello ',str_user));
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
                dir_saveExcel                   = uigetdir('*.mat*', 'Please Select Directory Where you keep Experiment .mat files');
                dir_saveMatfiles                = uigetdir('*.mat*', 'Please Select Directory Where you keep Experiment .mat files');
                [dir_umbrella, ~, ~]            = fileparts(dir_Xfolder); %Identify Presumptuous 
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
input_preinjury                             = 0;
switch input_preinjury
    case {1}
        input_40c                                   ='y'; %I am assuming no 40 circles on preinjury
        [ num_PreDf,num_PreBg,num_PreRw,num_PreMF,Img0,Tr0]  = function_preInjuryCTF(dir_Xfolder,str_positionName,BPP,input_40c);
        Prex                                        =[num_PreDf,num_PreRw,num_PreBg];
    otherwise
        warning ('There is no PreInjury Images')
        num_PreDf = NaN; %use to be PreDf --need to change this
        num_PreBg = NaN;
        num_PreRw = NaN;
        num_PreMF = NaN;
        Prex      = NaN;
        Img0      = NaN;
        Tr0       = NaN;
end;

disp        ('END: Part 1: PreInjury ROS Calculation')
%% Quantitate ROS in Region of Injury after Injury
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
cd(dir_positionFolder)
%--Strings
dir_positionName        = cd;
d_positionName          = dir('*T*');
len_d                   = length(d_positionName);
%----
Arr_Bkgd        = {0};
Arr_Ifish       = {0};
Arr_Raw         = {0};
Arr_Difference  = {0};
z_i             = 1;
I_wr            ={0};
I_br            ={0};
I_brs           ={0};
MeanBgs         =(0);
IntDen          =(0);
CTCF            =(0);
RoiArea         =(0);
I_fish          ={0};
pixelNum         =(0);
MeanFluorescence=(0);
disp('Processing...')
disp(str_positionName)
%--
Time                = round((arr_T(1:num_frames)*SamplingFrequency+MPI_start)-SamplingFrequency);
timepoints      =[60,90,120];
len_n=length(timepoints);
for j=1:len_n
    val             = timepoints(j); % value to find
    tmp             = abs(Time-val);
    [idx idx]       = min(tmp); % idex of closest value
    closest         = Time(idx); %closest value
    tp(j)           = idx;
end;
T_p=cat(2,(1),tp);
%-
iteration       = 200;
everyframe      = (1:num_frames);
everyiteration  = (1:iteration:num_frames);
golditeration   = cat(2,T_p,everyiteration);
Threshold       = 7000;
%for i=(arr_LoopGT) %This was the do the 30 minute interations
for i=everyframe % This is to do every single frame 
    %% Set Parameters
    cd(dir_positionName) %Go back to experiment folder
    tfolder         = d_positionName(i).name; %select Frame folder 
    cd(tfolder) %move to this directory
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
        disp('Determining best z_stack');
        for ii=1:length(f_TexRed)
            figure(ii);imagesc(imread(f_TexRed(ii).name))
        end;
        z_i = input('Whats the best z-stack');
    else
        z_i = 1;
    end;
    everyother  = i;       
    seeMeTR       = f_TexRed(z_i).name;
    seeMeBF     = f_BF(z_i).name;
    I_TR          = imread(seeMeTR); 
    I_BF        = imread(seeMeBF);
    I_BF        = imadjust(I_BF);
    happy       = 'n';
    %% Qualitative
    OriginalImgTR            = I_TR; 
    if ~isempty(find(i==golditeration))%i==1
        disp(T_p)
        disp(i)
        %-Very first iteration to set all the presets
           %%-Thresholding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while happy ~= 'y'
                close all
                A               = I_TR;
                A(A>=Threshold) = 1;
                A(A>1)          = 0;
                imagesc(A);
                save_Threshold  = Threshold;
                
                Threshold       = input(strcat('Set the Low Threshold, [',num2str(Threshold),']'));
                if isempty(Threshold)
                    Threshold = save_Threshold;
                else
                end;
                close all
                A(A>=Threshold) = 1;
                A(A<1)          = 0;
                imagesc(A)
                happy           = input('Are you happy With Threshold','s');
            end;
            %%-REGION Selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if i==1            
                Img1            = I_BF;
                disp('Select ROS Wound Region')
                close all; imagesc(I_BF); colormap('gray');title('ROS REGION');
                RosRegion       = roipoly();
                RosRegionMask   = RosRegion; %added to include in RawPixelIntensity 1/3/2018
                CellPixelNum    = sum(sum(RosRegionMask));
                %------------Background----------------------------------
                disp('Select Background Region')
                boo_bg     = 'y';
                a          = 0;
                temp_Bgs     =(0);
                [x,y]      = size(A);   
                blankimage = zeros(x,y);
                while (isempty(boo_bg)) ||  boo_bg=='y'  
                    a=a+1;
                    if a==1
                        BgRegions           = blankimage;
                    else
                        boo_bg=input('Another Background Region?','s');
                        if boo_bg=='n'
                            break
                        elseif isempty(boo_bg)
                            break
                        end;
                    end;
                    disp('Select Background Region')
                    close all; imagesc(I_BF); title('BACKGROUND REGION');
                    BackgroundRegion         = roipoly();
                    OriginalImgTR            = I_TR; 
                    BgRegions                = BgRegions & BackgroundRegion;
                    BackgroundPixelNum       = sum(sum(BackgroundRegion));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/
                    BackgroundPixelIntensity = sum(sum(uint16(BackgroundRegion).*I_TR))/BackgroundPixelNum;
                    temp_Bgs(a)=BackgroundPixelIntensity;
                    disp(temp_Bgs)
                    close all;


                end;
                close all;
            else
                disp(strcat('Keeping Previous Background Region on frame: ',num2str(i)))
                close all; 
                temp_outlinebg  = bwperim(I_brs{i-1});
                temp_outline    = bwperim(I_rr{i-1});
                temp_segout     = imoverlay(I_TR,temp_outline,[1 0 0]);
                temp_segouts    = imoverlay(temp_segout,temp_outline,[1 1 0]);
                temp_segout_BF  = imoverlay(I_BF,temp_outline,[1 0 0]);
                close all;
                figure();imagesc(temp_segout);
                figure();imagesc(temp_segout_BF);colormap('gray'); title('ROS REGION');
                rr_input=input('Do ROS again?','s');
                if rr_input=='y'
                    RosRegion       = roipoly();
                else
                end;             
            end;
%     elseif ~isempty(find(i==everyiteration))
        %-Will need to check to see if everything is in order       
    else
        %-These are automated processing.
            %-Thresholding
                A               = I_TR;
                A(A>=Threshold) = 1;
                A(A>1)          = 0;
            %- Region Selection
                RosRegion       = I_rr{i-1};
                BackgroundRegion= I_br{i-1};
                BgRegions       = I_brs{i-1};
                temp_Bgs        = MeanBgs(i-1);
                
    end;
    RawPixelIntensity           = sum(sum(uint16(RosRegionMask).*OriginalImgTR))/CellPixelNum;
    B         = A;
    I_fish{i} = I_TR;
    I_rr{i}   = RosRegion;
    I_br{i}   = BackgroundRegion;
    I_brs{i}   = BgRegions;
    MeanBgs(i)= mean(temp_Bgs);        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %% Quantitation
    arr_Bkgd  =(0);
    arr_Raw   =(0);
    arr_Diff  =(0);
    for j=1:length(f_TexRed)
        seeMe_TR1               = f_TexRed(j).name;
        I_TR1                   = imread(seeMe_TR1);  %read image to beMe
%         B=I;
%         B(B>=Threshold)=1;
%         B(B<1)=0;B=~B;B=~B;  
        se90                = strel('line', 3, 90); se0 = strel('line', 3, 0);            
        Bdil                = imdilate(B, [se90 se0]); %figure, imshow(Adil), title('dilated gradient mask');
        BWdfill             = imfill(Bdil, 'holes'); %figure, imshow(BWdfill); title('binary image with filled holes')             
        bw                  = bwareaopen(BWdfill, 5000);%Eliminating objects fewer than 5000 Pixels
        bw2                 = bwmorph(bw,'clean'); %clean not sure, can't hurt
        bw2                 = bwmorph(bw2,'thicken'); %thickens by just a little
        bw                  = bw2;
        bwO                 = bwperim(bw);
        bwO                 = imdilate(bwO,[se90,se0]);
        bw_color_realsegout = imoverlay(I_TR1,bwO,[1 0 0]);%high is red
        %%%Converse is the Background
        se90bck             = strel('line', 20, 90);
        se0bck              = strel('line', 20, 0);
        bw_background       = imdilate(BWdfill, [se90bck se0bck]);
        bw_background       =~bw_background;
            
        %% To see the outline over the ROS Profile
        RR                          = RosRegion&bw;
        RR_outline                  = bwperim(RR);
        RBg_outline                 = bwperim(BackgroundRegion);
        SegoutIntensity             = I_TR1;
        SegoutIntensity(RR_outline) = BPP;
        Fish                        = imadjust(I_TR1, stretchlim(I_TR1), []);
        I_fishRR                    = imoverlay(Fish,RR_outline,[1 0 0]);
        I_fishROIs                  = imoverlay(I_fishRR,RBg_outline,[1 1 0]);
        imshow(I_fishROIs);  
        text( 5, 5, 'Wound Region','Color','red')
        text( 5,15,'Background Region','Color','yellow')
        %% CTCF Calculation
        cc_rr     = regionprops(RR);
        temp_area = [cc_rr.Area];
        RoiArea(i)= sum(temp_area);
        IntDen(i) = RoiArea(i) * RawPixelIntensity; %Integrated Density
        CTCF(i)   = IntDen(i) - (RoiArea(i) *MeanBgs(i));%Corrected Total Cell Fluorescence 
        iterations=4;
        diameter=5;     
        if i==1
            [MMPI, mask_combo] = fortycircles(OriginalImgTR,bw,diameter,iterations);
            %MeanFluorescence(i)= MMPI;
        else
            %alternative possibility.  Get each individual object and do
            %average 
        end
        c40           = mask_combo & bw;
        pixelNum(i)       = sum(sum(c40));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/
        pixelIntensity = sum(sum(uint16(c40).*I_TR))/pixelNum(i);
        MeanFluorescence(i)            = pixelIntensity;
        %% 
        %-PreSet
        OriginalImgTR               = I_TR1;
        RosRegionMask               = RR;
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
        Arr_Ifish{i}             = I_fishROIs;
        Arr_Raw{i}               = arr_Raw;
        Arr_Bkgd{i}              = arr_Bkgd;
        Arr_Difference{i}        = arr_Diff;
        Bg(i,1)                   = arr_Bkgd(z_i);
        Df(i,1)                   = arr_Diff(z_i);
        Rw(i,1)                   = arr_Raw(z_i);
        disp(strcat('Completed: ',num2str((i/num_frames)*100),'% [',num2str(i),'] out of [',num2str(num_frames),']'))
        %     copyfile(combPh2img.name,combinePh2cd_s*tr)
      
end;
disp('END: ROS Detection Injury Loop')
%% --Save Values
disp(strcat('Background: ',num2str(Bg)))
disp(strcat('Raw: ',num2str(Rw)))
disp(strcat('Difference: ',num2str(Df)))
PixValues           = 0;
str_PixValues       = ['Background','Raw','Difference'];
PixValues           = [Df,Rw,Bg];
Time                = round((arr_T(1:num_frames)*SamplingFrequency+MPI_start)-SamplingFrequency);
% Adjust Values
adjBg=Bg;
idx=find(Bg>(std(Bg)+mean(Bg)));
adjBg(idx)=mean(Bg);
len_Bgnoise=length(idx);
if len_Bgnoise>6
    disp('Lenghth of Bg noise higher than SD is more than 5%')
elseif len_Bgnoise>12
    disp('Lenghth of Bg noise higher than SD is more than 10%')
end;
%-
adjDf=Df;
adjRw=Rw;
adjDf=adjRw-adjBg;
adjPixValues=[adjDf,adjRw,adjBg];
%-
dBaseline=((adjDf-num_PreDf)/num_PreDf)*100; %percent change from baseline
%Remember arr_T is the 1:frames
disp('END: After Injury ROS Calculation')
%%
disp('Want to save???????')
close all
plot(MeanFluorescence);
pause()
cd('C:\Users\AndreDaniel\OneDrive\PhD Data\m_ROS');
namename        =strcat(name,str_positionName);
save(namename,'dBaseline','adjPixValues','str_PixValues','PixValues','Time','frames','T_p','Prex','name5','name','t_plate','ti_d','golditeration','Img1');
disp('END:Saved m_Ros Values')
%% CTCF No wound
nowound_IntDen  = (pixelNum).* MeanFluorescence; %Integrated Density
nowound_CTCF            = nowound_IntDen - ((pixelNum).*MeanBgs);%Corrected Total Cell Fluorescence 
canvas          = imoverlay(I_TR,mask_combo,[1 0 0]);
mmpi0=MeanFluorescence;
figure();imshow(canvas);
disp('Load the MeanFlourescence of Non Wound')
cd('C:\Users\AndreDaniel\OneDrive\PhD Data\m_ROS\CTCF');
nowoundname    =strcat(namename,'noWound');
save(nowoundname,'mmpi0','arr_T','MeanFluorescence','ixelNum','nowound_IntDen','MeanBgs','nowound_CTCF','T_p','Time','frames','name5','name','t_plate','ti_d','golditeration','Img1');
disp('END: Saved noWound CTCF');
