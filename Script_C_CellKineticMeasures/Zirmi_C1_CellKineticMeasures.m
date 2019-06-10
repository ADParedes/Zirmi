%% Zirmi Package - Module (3-C1) Computing Cell Migration based Metrics
% Version 1.4 2019-06-03
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% (A) Is Initial Formatting of Pre-processing Data
%- This is intended to take data Organized in manner outlined and place in
%a directory that the Zirmi program can use to throughput computations
%% (B) Is PhagoSight Based Scripts and initial Processing of Fish
%- Contains the following Scripts
%- S0, Loading D
%% (C) Is Computing Cell Migration based Metrics 
%- PLEASE NOTE THE FOLLOWING:
% 1.  Golden times refers to doing 60,90 and 120 time points.  This is not
% really editable in this script (should be in Version 2.0).
% 2.  Track Inclusion Criteria 1A is select tracks that are % for all the
% "interval" outcome measures
% 3.  Track Inclusion Criteria 1B is to select tracks of those initial % of
% each timepoint - so will have less cells, if they arn't present in later
% time point tracks.  
% "SM" outcome measures and "Cumulative" outcome measures so that we can
% do paired analysis.  

% that are also 90% in each time frame.  
%% (D) Is to Quantitate Corrective Total Fluorescence Quantitation

%% Summary
% NOTE WE TAKE NODE NETWORK and Calculate all the Metrics in Categories
% field that I do not use (ANDRE) do i use phagosight handles metrics
% I take the positions and recalculate velocity, dist travel all the
% categories below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This will break up int 60,90,120,150,180 and maximum frame based
%information
% Two categories  #1 All of fish  #2 In Wound
%Category #1
% We will determine distance traveled, velocity (* non phagosight
% based),meandering ratio, tortuosity, and static ratio
%Category #2
%We determine if inWound from new made woundRegion
%NOTE:  nanmean([ ]) will average ignoring nan values!!!!
%% Warning
disp                    ('Please Note the Following Warnings : ')
warning                 ('Zirmi A & B0 needs to be ran on experiment directory before proceeding')
warning                 ('Zirmi: Need to have handles.distmap and accurate spatial domains 1 & 2')%Spatial Domain 1 is Wound GAP Spaital domain 2 is Tip of notochord
warning                 ('Zirmi: If Removing Movement remember to save NEW handles')
PARAMETERS.ParameterS  = 100;      %Leukocyte Spatial Interval (150um)
display                 ('Velocity is represented by pixel/frame!!!!'); %there is Utility in correcting this (12/13/2017)
pause                   (2);
%% In case clauses to ensure everything that needs to be loaded was loaded
switch exist('POI') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 1: PASS! POI exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('No Registered POI  variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run Zirmi B0 before proceeding')
        display                   ('Module 3 DISCONTINUED') 
        return
end;
switch exist('PARAMETERS') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 2: PASS! PARAMETERS exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('No Registered  PARAMETER variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run Zirmi B0 before proceeding')
        display                   ('Module 3 DISCONTINUED') 
        return
end;

switch exist('PhagoSight') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 3: PASS! PhagoSight exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('No Registered  PhagoSight variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('Run B0, B1, and PhagoSight before proceeding')
        display                   ('Module 3 DISCONTINUED') 
        return
end;

switch exist('handles') %#ok<EXIST>
    case{1} % POI does exist
        display                     (strcat('CHECK 4: PASS! handles exists'))
    otherwise %POI Does NOT exist - User did not follow execution procedure
        display                   ('No Registered  PhagoSight variables')
        warning                   ('Zirmi Cannot Proceed if not Metadata is not recognized')
        warning                   ('You need to run PhagoSight before proceeding')
        display                   ('Module 3  DISCONTINUED') 
        return
end;
%% Verify Parameters
prompt                  = {'Module 3: Track Inclusion Criteria 1A (based on centriod detection in ALL FRAMES; 1.0- 0.0 )'...     % 1
                            ,'Module 3: Track Inclusion Criteria 1B (based on centriod detection per Time Domain; 1.0- 0.0 )'... % 2
                            ,'Module 3: Track Inclusion Criteria 2 Time Domain T1,T2, and T3 (minutes)'...                       % 3
                            ,'Module 3: Track Inclusion Criteria 3 Space domain S2,S3, and S4  distance (um ; micrometer)'...    % 4
                            ,'Enter Parameter 6: StaticLimit (um ; micrometer)'};                                                % 5
dlg_title               = 'confirm Module 3 customizable algorithm parameters';
num_lines               = 1;
defaultans              = {num2str(PARAMETERS.ParameterA),num2str(POI.Parameter_gtA),num2str(30),num2str(PARAMETERS.ParameterS),num2str(PARAMETERS.ParameterB)};
answer                  = inputdlg(prompt,dlg_title,num_lines,defaultans);

PARAMETERS.ParameterA              = str2num(answer{1}); % Trackability all frames
POI.Parameter_gtA                  = str2num(answer{2}); % Trackability per time domain
PARAMETERS.ParameterB              = str2num(answer{5}); % Static Ratio
PARAMETERS.ParameterS              = str2num(answer{4});
%% Load PARAMETERS, POI, PhagoSight & ADP  Needed 
%Note: Phagosight DNE PhagoSight (two separate variables)
ParameterA              = PARAMETERS.ParameterA;   % Majority Track Percent (MTP)
ParameterB              = PARAMETERS.ParameterB;   % Staticity Threshold
ParameterS              = PARAMETERS.ParameterS;   % Spatial Threshold - 100um is default

micronsPerPixel         = PARAMETERS.Parameter2;   % LateralPixelResolution  
micronsPerStack         = PARAMETERS.Parameter3;   % micronsPerStack   
SamplingFrequency       = PARAMETERS.Parameter4;   % Sampling Frequency
MPI_start               = PARAMETERS.Parameter5;   % Minutes post injury (Start of Imaging) ; old: t_plate
tplate                  = MPI_start;  %stored original t_plate as tplate incase I need to refer back to it

Parameter_gtA           = POI.Parameter_gtA;       % MTP for all each Temporal segment
dir_EOI                 = POI.Parameter10c;        % Experiment Directory (E.g. AB020 within Confocal)
str_FOI                 = POI.Parameter11c;        % ExperimentFish of Interest (E.g. P2 within AB020) 
valsPh2                 = POI.Parameter13;      

woundRegion             = PhagoSight.woundRegion;
wR1                     = PhagoSight.wR1;          % Wound Gap Epicenter
wR2                     = PhagoSight.wR2;          % Notochord Tip Epicenter
disp('End: Paramters load')
%% Simple Calculations (Pixel based information)
coe                     =(micronsPerPixel/SamplingFrequency); %microns per minute
zstack                  = length(valsPh2);
ZstepPixel              = micronsPerStack/micronsPerPixel;
%Example of the equation (minum is 1 maximum is 3)  so 1 is 0 and 3 is 44.
zval                    = zstack*ZstepPixel-ZstepPixel; %this is just an example, so just take z and multiply by ZstepPixel
disp('chose zstack based on BF/Phase Contrast'); %..not always accurate - Fixed.
%% Clearing the Important Structure Arrays
clear Phagosight Andre
%% Calculate Frame positions for Time 60,90,120,150,180(max), full.  
%Access the Position Folder (E.g. P1 or P2..) to determine Frames used
cd(dir_EOI)
cd(str_FOI)
%Save the directory and access folder names to determine Frames
dname5                  = dir;
len_dname               = length(dname5);
a=0;
clear arr_T
clear tfolder
for i=(3:len_dname) %for MAC is at 5.  
    a       = a+1;
    tfolder = dname5(i).name;
    arr_T(a)= str2num(tfolder(2:end));
    min_t   = min(arr_T);
    max_t   = max(arr_T);
end;
%min to add to t_plate so that deleted frames are considered in time 
adj_tplate              =(min_t-1)*SamplingFrequency+MPI_start;
tplate                  = adj_tplate;
%max is just arr_T length
frames                  = length(arr_T); % or frames=handles.numFrames;
%handles based frame and track
[frameinWound,track]    = size(handles.inWound);
%% Determine the time points after wounding.  
% GT_45   = round((45-MPI_start)/SamplingFrequency); %%%%%%%%%%%%%%%%%%%%CHANGING
GT_45   = 3; %Good So I start after 3 ---> This is only for macrophages - neutrophils I fixed this
GT_60   = round((60-tplate)/SamplingFrequency);
GT_90   = round((90-tplate)/SamplingFrequency);
GT_120  = round((120-tplate)/SamplingFrequency);
GT_150  = round((150-tplate)/SamplingFrequency);  %GT_150=golden(2hours and 30 minutes)
GT_180  = round((180-tplate)/SamplingFrequency);
boo_60  = 1;
boo_90  = 1;
boo_120 = 1;
boo_150 = 1;
boo_180 = 1;
boo_max = 1;
display('_')
if GT_60>frames
    display(strcat(num2str(GT_120),'_ GT_60 exceeds ',' _', num2str(frames)))
    boo_60      =0;
    boo_90      =0;
    boo_120     =0;
    boo_150     =0;
    boo_180     =0;
elseif GT_90>frames
    display(strcat(num2str(GT_120),'_ GT_90 exceeds ',' _', num2str(frames)))
    boo_90      =0;
    boo_120     =0;
    boo_150     =0;
    boo_180     =0;
elseif GT_120>frames
    display(strcat(num2str(GT_120),'_ GT_120 exceeds ',' _', num2str(frames)))
    boo_120     =0;   
    boo_150     =0;
    boo_180     =0;
elseif GT_150>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_150     =0;
    boo_180     =0;
elseif GT_150>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_150     =0;
    boo_180     =0;
elseif GT_180>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_180     =0;
else
    display('you have over 3 hours of total experiment.  Wow')
end;  
arr_boo_GT              =[boo_60,boo_90,boo_120,boo_150,boo_180,boo_max];
arr_GT                  =[GT_60,GT_90,GT_120,GT_150,GT_180,frames]; GT_max=frames;
str_tplate              = num2str(MPI_start);
str_GT                  =[str_tplate,strcat('60min:Frame ',num2str(GT_60)),strcat('90min:Frame ',num2str(GT_90)),strcat('120min:Frame ',num2str(GT_120)),strcat('150min:Frame ',num2str(GT_150)),strcat('180min:Frame ',num2str(GT_180)),strcat('Frame:Frame ',num2str(frames))];
arr_Time                = round((arr_T(1:frames)*SamplingFrequency+MPI_start)-SamplingFrequency)';
disp('End: Times')
%% GT in wound Count
gtinwoundcount=0;
beat=0;
for I=1:length(arr_boo_GT)
    check_boo_GT=arr_boo_GT(I);
    gtinwoundcount=0;
    gt_spot=arr_GT(I);
    if check_boo_GT==1
        for i=1:track
           if length(1:gt_spot)<frameinWound        
                place=handles.inWound(1:gt_spot,i)';
               beat=sum(place);
               if beat>0
                   gtinwoundcount=gtinwoundcount+1;
               else
               end;
           else
               % incase first frame is beyond in wound frame
               if I==1
                   gtinwoundcount=nan;
               else
                    gtinwoundcount=arr_GTinWound(I-1);
               end;
           end;
        end;
    else
        gtinwoundcount=nan;
    end;
    arr_GTinWound(I)=gtinwoundcount;
end;
CTTM=(SamplingFrequency*handles.numFrames)+tplate;%Corrected_Total_Time_Max
imaginglength=(SamplingFrequency*handles.numFrames);
P=str_FOI;
disp('End GT Calculations')
%% GT NodeNetwork
HNN={handles.nodeNetwork}; %Node Network {}
Hnn=handles.nodeNetwork;   %Node Network ()
[m,n]=size(Hnn);
[~,~]=size(handles.finalNetwork); %Num_TrackMax is  the same as mp
Col06=handles.nodeNetwork(:,6); %Unique ID
Col14=HNN{1}(:,13); %Track  
Col05=HNN{1}(:,5); % Time Frame
Col07=HNN{1}(:,7);%Parent
Col09=HNN{1}(:,9);% velocity
Col10=HNN{1}(:,10); %Volume of MP
%Array limit is set by 'arr_boo_GT'
for i=1:length(arr_GT)
    frame_limit=arr_GT(i);
    [idx_over]=find(Col05>frame_limit);
    [idx_in]=find(Col05<frame_limit);
    Arr_hnn_idxgt{i}=idx_in;
    Arr_HNN_GT{i}=Hnn(idx_in,1:n);
end;

disp('End GT NodeNetwork')
%% Phagosight Unique ID Track,Frame & hop per Track, Recalculate HOP
%-Structure Array: Phagosight
Num_TrackMax    = max(Col14); %This is the maximum number of tracks in handles
Num_FrameMax    = handles.numFrames; %phagosight does 1 over, so its frames+1
Arr_uIDtrack    ={0};%Unique ID per Track derived from handles.finalNetwork
NN_uIDtrack     ={0}; %Unique ID per Track derived from handles.nodeNetwork
NN_uIDframe     ={0};%Unique ID per TimeFrame derived from handles.nodeNetwork
Arr_DhopTrack   ={0};%DistancePer"Hop" per Track derived from phagosight stats
% manually determining from node network the hops that land in each track
% absolete since it ends up being the same thing.
clear Zall_disthops
for i=1:Num_TrackMax
    idx_uid=find(Col14==i);
    NN_idx_uidframe{i}=Col05(idx_uid);
    NN_uIDtrack{i}=idx_uid; %based on handles.nodenetwork (NN) 
    clear zdisthops
        for j=2:length(idx_uid) 
            if j==1
            else
                n1=NN_uIDtrack{i}(j-1);
                n2=NN_uIDtrack{i}(j);
                x_node1=handles.nodeNetwork(n1,1);
                y_node1=handles.nodeNetwork(n1,2);
                z_node1=handles.nodeNetwork(n1,3)*ZstepPixel-ZstepPixel;
                x_node2=handles.nodeNetwork(n2,1);
                y_node2=handles.nodeNetwork(n2,2);
                z_node2=handles.nodeNetwork(n2,3)*ZstepPixel-ZstepPixel;
                zdisthops(j)=pdist([x_node1,y_node1,z_node1; x_node2,y_node2,z_node2],'Euclidean');
            end;
        end;
        Zall_disthops{i}=zdisthops';
end;
Phagosight.Zdistofhops=Zall_disthops;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%using Phagosights arrays for distance per hop and hops in each track
%%%%%%%%%%%%%%%%%%%%final entwork and distmap
for i=1:Num_TrackMax
    a=handles.finalNetwork(:,i); %Sort UinqueID by Track
    d=handles.distMaps.absDistPerHop(:,i); %abs Distance per Hop
    a(a==0)=[];
    len_a=length(a);%     d(d==0)=[]; 
    Arr_uIDtrack{i}=a;
    Arr_DhopTrack{i}=d(1:(len_a-1)); %Note if 0 value is at some ends, then final frame is duplicated --ERROR
    clear a
    clear d
end;
%sizing Both up (so far so good)
check_uIDtrack  = [NN_uIDtrack',Arr_uIDtrack']; %if they have same index I trust matlab's final network

for i=1:Num_FrameMax
        b               = find(Col05==i); %Sort UniqueID by Time Frame
        NN_uIDframe{i}  = b;
        clear b
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('End Phagosight Unique ID Track,Frame & hop per Track and REcalculated HOP!!! :)')
%% Andre: By Time Frame and Volume
%-Structure Array: Andre 
%NN_framesinTrack
NN_framesinTrack    ={0}; % The frames in Each Track, therefore {1} will give frames in track 1. 
arr_lensFrame       = 0; %Total Time Length (Max-Min)
NN_VolumeinTrack    ={0};%Volume per Track (MP)
arr_meanVolume      =(0);
for J=1:Num_TrackMax
    for j=Arr_uIDtrack{J}
        c                       = HNN{1}(j,5);
        v                       = HNN{1}(j,10);
    end;
    minC                    = min(c);
    maxC                    = max(c);
    NN_framesinTrack{J}     = c;  %Indivual Frames per Unique ID segmented in tracks
    NN_VolumeinTrack{J}     = v;
    e                       = maxC-minC;
    arr_lensFrame(J)        = e; %Andre Real #tot frames in track
    arr_meanVolume(J)       = mean(v);
    V_std(J)                = std(v);
    clear v
    clear c
    V_radius(J)             =((arr_meanVolume(J)-V_std(J))*(3/4)*(1/pi))^(1/3);
end;
%%%%%%%%%%%%%%%%ANDRE%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Andre.meanVolPerMp      = arr_meanVolume;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arr_lensHop             = handles.distanceNetwork.numHops';
arr_lensFrame           = arr_lensFrame';
Andre.lengthframeofTrack= arr_lensFrame;
disp('End Phagosight Frame and Volume & Struct Andre')
%% ----------------------------------------Pass1-------------------------------------------------------
%-More Structure Array: Phagosight
arr_timediff        =0;
Arr_TimeDiscrepency ={0};
Arr_absVelocity     ={0};
ARR_TrackVel        =0;
ARR_ComboCompare    ={0};
for K=1:Num_TrackMax
    display(strcat('Track:',num2str(K)))
    len_C                   = length(NN_framesinTrack{K});
    minC                    = min(NN_framesinTrack{K});
    maxC                    = max(NN_framesinTrack{K});
    adp_timediff            =(minC:maxC);
    len_adp                 = length(adp_timediff);  
    for k=1:length(NN_framesinTrack{K})
        if k==1
%             val_timediff=0;
        else
            D_cat               = NN_framesinTrack{K}(k-1);
            cur                 = NN_framesinTrack{K}(k);
            val_timediff        = cur-D_cat;
            arr_timediff(k-1)   = val_timediff;
        end;
    end;
    zdistanceperhop         = Zall_disthops{K}(2:end);
    distanceperhop          = Arr_DhopTrack{K};
    Arr_TimeDiscrepency{K}  = arr_timediff';
    Arr_absVelocity{K}      = distanceperhop./arr_timediff';   
    Arr_absVelocityZ{K}     = zdistanceperhop./arr_timediff';  %%This is the new handles
    ARR_ComboCompare{K}     =[NN_framesinTrack{K}(2:len_C),Arr_TimeDiscrepency{K},Arr_absVelocity{K},Arr_DhopTrack{K},NN_VolumeinTrack{K}(2:len_C)]; 
    arr_timediff            = 0;
    ARR_TrackVel(K)         =(sum(distanceperhop)/arr_lensFrame(K)); %#ok<SAGROW>
    numHopsasd(K)           = len_C;
end;
Phagosight.absVelocityZ=Arr_absVelocityZ;
disp('End Determining Frame and Hop Discrepency- Arr_TimeDiscrepency')
%% Create a new handles.finalNetwork
%handles.finalNetwork Without Zstep Consideration
U=[0];
F={0};
H={0};
for L=1:Num_TrackMax
    F{L}=vertcat(U,Arr_absVelocity{L});
    H{L}=vertcat(U,Arr_TimeDiscrepency{L});
    minC=min(NN_framesinTrack{L});
    maxC=max(NN_framesinTrack{L});
    nodecount=1;
    cc=minC-1;
    for l=1:length(NN_framesinTrack{L})
        cc=cc+1;
        nodecount=nodecount+1;
        aa=find(NN_framesinTrack{L}==cc);
        bb=isempty(aa);
        if cc>maxC
            Andre.distmap1(cc,L)=0;
        elseif isempty(aa)
            O=ones(H{L}(l),1);
%             display('empty')
            for k=1:H{L}(l)
                Andre.distmap1(cc,L)=F{L}(l);
                cc=cc+1;
                nodecount=nodecount+1;
            end;
            cc=cc-1;
            nodecount=nodecount-1;
            Phagosight.distmap1(cc,L)=Arr_DhopTrack{L}(l-1);
            distmap.Cheng1(cc,L)=F{L}(l); 
        else %is not empty
            var_catch=F{L}(l);
            Andre.distmap1(cc,L)=var_catch; %Fills skipped hops with the average hop distance over skip
            distmap.Jun1(cc,L)=var_catch; %MOVES TRACE OF ANY DISTANCE TRAVELED ALL TOGETHER
            Phagosight.distmap1(cc,L)=var_catch; %Provides the average of One of the missing spots ignoring rest.
            distmap.Cheng1(cc,L)=var_catch; %Provides the big jump of that frame gap but aknoweldging the rest

        end; 
    end;
    lastNodeFrame(L)=nodecount;
end;
Andre.lengthframeofhops=numHopsasd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2016-09-22 with Z step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zdofp={0};
U=[0];
F={0};
H={0};
for L=1:Num_TrackMax
    Zdofp{L}=Phagosight.Zdistofhops{L}(2:end); %includes Z step, that a boiiiiiiii
    F{L}=vertcat(U,Arr_absVelocityZ{L});
    H{L}=vertcat(U,Arr_TimeDiscrepency{L});
    minC=min(NN_framesinTrack{L});
    maxC=max(NN_framesinTrack{L});
    nodecount=1;
    cc=minC-1;
    for l=1:length(NN_framesinTrack{L})
        cc=cc+1;
        nodecount=nodecount+1;
        aa=find(NN_framesinTrack{L}==cc);
        bb=isempty(aa);
        if cc>maxC
            Andre.distmap(cc,L)=0;
        elseif isempty(aa)
            O=ones(H{L}(l),1);
%             display('empty')
            for k=1:H{L}(l)
                Andre.distmap(cc,L)=F{L}(l);
                cc=cc+1;
                nodecount=nodecount+1;
            end;
            cc=cc-1;
            nodecount=nodecount-1;
            Phagosight.distmap(cc,L)=Zdofp{L}(l-1);
            distmap.Cheng(cc,L)=F{L}(l); 
        else %is not empty
            var_catch=F{L}(l);
            Andre.distmap(cc,L)=var_catch; %Fills skipped hops with the average hop distance over skip
            distmap.Jun(cc,L)=var_catch; %MOVES TRACE OF ANY DISTANCE TRAVELED ALL TOGETHER
            Phagosight.distmap(cc,L)=var_catch; %Provides the average of One of the missing spots ignoring rest.
            distmap.Cheng(cc,L)=var_catch; %Provides the big jump of that frame gap but acknoweldging the rest

        end; 
    end;
    lastNodeFrame(L)=nodecount;
end;
Andre.lengthframeofhops=numHopsasd;
distmap.JunNan=distmap.Jun;
distmap.JunNan(distmap.JunNan==0)=nan;
disp('End Create distmaps distmap.Jun)')
%% Velocity of Andre and Phagosight based distmaps and GTs
%%%-Phagosight excludes the skips including the final determination point
%%%-Andre.distmap gives distance traveled, including assumped average as
%%%-distmap.Jun ignores all skips 
%This is creates velocity based on total time
clear GTs 
clear Jhopdisttrack 
clear J_hopdisttrack
jhdt=0;
Arr_VelCompare=0;
for J=1:length(arr_GT)
    len_gt=arr_GT(J); %length of frame to use
    % Update 2019-06-03 Frames lenght cannot exeed distmap length
    if J==6
        if len_gt > length(Andre.distmap)
            disp('Projected Time segment is beyond the length of imaging')
            len_gt = length(Andre.distmap);
            arr_GT(J) = len_gt;  %This should only be true for J =6
        end;
    end;
        
        if arr_boo_GT(J)==1             
             for I=1:Num_TrackMax           
                temp_adm=Andre.distmap(1:len_gt,I);
                temp_adm(temp_adm==0)=[];
                Andre.DistPerMp(I)=sum(temp_adm);
                Andre.VelPerMp(I)=mean(temp_adm);
                temp_pdm=Phagosight.distmap(1:len_gt,I);
                temp_pdm(temp_pdm==0)=[];
                Phagosight.DistPerMp(I)=sum(temp_pdm);
                Phagosight.VelPerMp(I)=mean(temp_pdm);
                temp_jdm=distmap.Jun(1:len_gt,I);
                temp_jhdt=temp_jdm;
                Jhopdisttrack{I}=temp_jdm;
                temp_jdm(temp_jdm==0)=[];
                distmap.JunDistPerMp(I)=sum(temp_jdm);
                distmap.JunVelPerMp(I)=mean(temp_jdm);
                %ajsdf;alsdkjf;lasdjf;lk
                 if I==1
                    jhdt=temp_jhdt;
                 else
                    jhdt=[jhdt temp_jhdt];
                 end;
             end;
             arr_jhdt{J}=jhdt;
             J_hopdisttrack{J}=Jhopdisttrack;
        else
                arr_jhdt{J}=nan;
                J_hopdisttrack{J}=nan;
                Phagosight.DistPerMp=nan(1,Num_TrackMax);
                Phagosight.VelPerMp=nan(1,Num_TrackMax);
                distmap.JunDistPerMp=nan(1,Num_TrackMax);
                distmap.JunVelPerMp=nan(1,Num_TrackMax);
        end;
    
    GTs.Arr_distpermp{J}=distmap.JunDistPerMp;
    GTs.Arr_velpermp{J}=distmap.JunVelPerMp';
    Phagosight.GTs_ArrDistperMp{J}=Phagosight.DistPerMp';
   if J==1
        phag_gtsdist= Phagosight.DistPerMp';
        phag_gtsvel=Phagosight.VelPerMp';
        gtsdistpermp=distmap.JunDistPerMp';
        gtsvelpermp=distmap.JunVelPerMp';
   else
        phag_gtsdist=[phag_gtsdist Phagosight.DistPerMp'];
        phag_gtsvel=[phag_gtsvel Phagosight.VelPerMp'];
        gtsdistpermp=[gtsdistpermp distmap.JunDistPerMp'];
        gtsvelpermp=[gtsvelpermp distmap.JunVelPerMp'];
    end;
end;
GTs.distmap.distmap=arr_jhdt;
GTs.Arr_distmap=J_hopdisttrack;
GTs.arr_velpermp=gtsvelpermp;
GTs.arr_distpermp=gtsdistpermp;
Phagosight.GTs_arrdistpermp= phag_gtsdist;
Phagosight.GTs_velpermp=phag_gtsvel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This finds velocity based on 30minute iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for J=1:length(arr_GT)
    len_gt=arr_GT(J);
    if arr_boo_GT(J)==1 
        if J==1   
            for I=1:Num_TrackMax
                temp_adm=Andre.distmap(1:len_gt,I);
                temp_adm(temp_adm==0)=[];
                Andre.DistPerMp30min(I)=sum(temp_adm);
                Andre.VelPerMp30min(I)=mean(temp_adm);
                temp_pdm=Phagosight.distmap(1:len_gt,I);
                temp_pdm(temp_pdm==0)=[];
                Phagosight.DistPerMp30min(I)=sum(temp_pdm);
                Phagosight.VelPerMp30min(I)=mean(temp_pdm);
                temp_jdm=distmap.Jun(1:len_gt,I);
                temp_jdm(temp_jdm==0)=[];
                distmap.JunDistPerMp30min(I)=sum(temp_jdm);
                distmap.JunVelPerMp30min(I)=mean(temp_jdm);
            end;  
            gtsvelper30min_1=distmap.JunVelPerMp30min;
        elseif J==length(arr_boo_GT) 
            for I=1:Num_TrackMax
                temp_adm=Andre.distmap(1:len_gt,I);
                temp_adm(temp_adm==0)=[];
                Andre.DistPerMp30min(I)=sum(temp_adm);
                Andre.VelPerMp30min(I)=mean(temp_adm);
                temp_pdm=Phagosight.distmap(1:len_gt,I);
                temp_pdm(temp_pdm==0)=[];
                Phagosight.DistPerMp30min(I)=sum(temp_pdm);
                Phagosight.VelPerMp30min(I)=mean(temp_pdm);
                temp_jdm=distmap.Jun(1:len_gt,I);
                temp_jdm(temp_jdm==0)=[];
                distmap.JunDistPerMp30min(I)=sum(temp_jdm);
                distmap.JunVelPerMp30min(I)=mean(temp_jdm);
            end;  
        else
              for I=1:Num_TrackMax
                    temp_adm=Andre.distmap(arr_GT(J-1)+1:len_gt,I);
                    temp_adm(temp_adm==0)=[];
                    Andre.DistPerMp30min(I)=sum(temp_adm);
                    Andre.VelPerMp30min(I)=mean(temp_adm);
                    temp_pdm=Phagosight.distmap(arr_GT(J-1)+1:len_gt,I);
                    temp_pdm(temp_pdm==0)=[];
                     Phagosight.DistPerMp30min(I)=sum(temp_pdm);
                    Phagosight.VelPerMp30min(I)=mean(temp_pdm);
                    temp_jdm=distmap.Jun(arr_GT(J-1)+1:len_gt,I);
                    temp_jdm(temp_jdm==0)=[];
                    distmap.JunDistPerMp30min(I)=sum(temp_jdm);
                    distmap.JunVelPerMp30min(I)=mean(temp_jdm);
              end; 
        end;
    else
        distmap.JunDistPerMp30min=nan(1,Num_TrackMax);
        distmap.JunVelPerMp30min=nan(1,Num_TrackMax);
    end;
    
    if J==1
        gtsdistper30min=distmap.JunDistPerMp30min';
        gtsvelper30min=distmap.JunVelPerMp30min';
    else
        gtsdistper30min=[gtsdistper30min distmap.JunDistPerMp30min'];
        gtsvelper30min=[gtsvelper30min distmap.JunVelPerMp30min'];
    end;
    GTs.Arr_distpermp30min{J}=distmap.JunDistPerMp30min';
    GTs.Arr_velpermp30min{J}=distmap.JunVelPerMp30min';
end;
    GTs.arr_velpermp30min=gtsvelper30min;
    GTs.arr_distpermp30min=gtsdistper30min;
% Andre.VelPerMp=Andre.VelPerMp';%Note this Velocity averages the skipped hops velocity
% Phagosight.VelPerMp=Phagosight.VelPerMp';%Note this Velocity takes skipped hops as 1 skewing velocity
% distmap.JunVelPerMp=distmap.JunVelPerMp';%Note this Velocity ignores all skips completely.
Arr_StrVelCompare=['Andre', ' Phagosight', ' Jun', ' distanceNetwork'];
Arr_VelCompare=[Andre.VelPerMp',Phagosight.VelPerMp',distmap.JunVelPerMp',handles.distanceNetwork.absVelocity'];
%M_mean=mean(M);
% distmap.Jun

disp('End Different Types of Cummulative Velocities (distmap.JunVelPerMp) and their Iterations ')
%% Velocity to Volume Correlation
%Velocity Correlation to Volume (none shown)
% % % % % % % % % % % % % % % % % % % % figure();
% % % % % % % % % % % % % % % % % % % % x=arr_meanVolume';
% % % % % % % % % % % % % % % % % % % % y=Andre.VelPerMp';
% % % % % % % % % % % % % % % % % % % % format long
% % % % % % % % % % % % % % % % % % % % b1=x\y;
% % % % % % % % % % % % % % % % % % % % yCalc1 = b1*x;
% % % % % % % % % % % % % % % % % % % % scatter(x,y)
% % % % % % % % % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % % % % % % % plot(x,yCalc1)
% % % % % % % % % % % % % % % % % % % % % ax=plot(arr_meanVolume,Andre.VelPerMp,'O','LineWidth',2,...
% % % % % % % % % % % % % % % % % % % % %            'MarkerEdgeColor','k',...
% % % % % % % % % % % % % % % % % % % % %            'MarkerFaceColor','g',...
% % % % % % % % % % % % % % % % % % % % %            'MarkerSize',10);
% % % % % % % % % % % % % % % % % % % %     title('Linear Regression Relation Between Velocity vs Volume','FontSize',16)
% % % % % % % % % % % % % % % % % % % %     ylabel('Velocity (frame/pixel)','FontSize',14)
% % % % % % % % % % % % % % % % % % % %     xlabel('Volume (pixel^3)','FontSize',14)
% % % % % % % % % % % % % % % % % % % %     grid on
% % % % % % % % % % % % % % % % % % % % % Improve the fit by including a y-intercept $\beta_0$ in your model as $y =
% % % % % % % % % % % % % % % % % % % % % \beta_0 + \beta_1x$. Calculate $\beta_0$ by padding |x| with a
% % % % % % % % % % % % % % % % % % % % % column of ones and using the |\| operator.    
% % % % % % % % % % % % % % % % % % % % X = [ones(length(x),1) x];
% % % % % % % % % % % % % % % % % % % % b = X\y;
% % % % % % % % % % % % % % % % % % % % % Visualize the relation by plotting it on the same figure.
% % % % % % % % % % % % % % % % % % % % yCalc2 = X*b;
% % % % % % % % % % % % % % % % % % % % plot(x,yCalc2,'-')
% % % % % % % % % % % % % % % % % % % %     legend('Data','Slope','Slope & Intercept','Location','best');
% % % % % % % % % % % % % % % % % % % % Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2); %Slope 1  
% % % % % % % % % % % % % % % % % % % % Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2); %Slope 2 (including y-intercept)
% % % % % % % % % % % % % % % % % % % % format short
disp('Linear Regression Relation Between Speed & Volume')
%% Velocity per Frame
N=0;
for I=1:Num_FrameMax
    temp_adf=Andre.distmap(I,:);
    temp_adf(temp_adm==0)=[];
    Andre.VelPerFram(I)=mean(temp_adf);
    temp_pdf=Phagosight.distmap(I,:);
    temp_pdf(temp_pdf==0)=[];
    Phagosight.VelPerFram(I)=mean(temp_pdf);
end;
N=[Andre.VelPerFram',Phagosight.VelPerFram'];
N_mean=mean(N);
disp('End Andre and Phagosight Velocity per Frame')
%% Static Ratio
clear Arr_StaticMap
clear chum
Bigger=[0];
for I=1:Num_TrackMax
%     temp_arr=find(Phagosight.distmap(:,I)>0);
    len=length(Arr_DhopTrack{I});
    rad=V_radius(I);
    for i=1:len
        val=Arr_DhopTrack{I}(i);
        if val>=rad;
%             display('Bigger')
            chum(i)=0;
        else
%             display('Smaller')
            chum(i)=1;
        end;
%         Bigger(a)=Phagosight.distmap(i,I)>V_radius(I)
    end;
%     display(num2str(I))
%     display(chum)
    Phagosight.MpStaticMap{I}=chum;
    Phagosight.MpStaticRatio(I)=mean(chum);
     clear chum
end;
%%Per Time 
chum=0;
for I=1:Num_FrameMax
    temp_arr=find(Phagosight.distmap(I,:)>0);
    len=length(temp_arr);
    if isempty(temp_arr)
        chum=0;
    elseif temp_arr==0 
        chum=0;
    else
        for i=1:len
            j=temp_arr(i);
                rad=V_radius(j);
            val=Phagosight.distmap(I,i);
            if val>=rad;
%                 display('Bigger')
                chum(i)=0;
            else
%                 display('Smaller')
                chum(i)=1;
            end;
    %         Bigger(a)=Phagosight.distmap(i,I)>V_radius(I)
        end;
    end;
%     display(num2str(I))
%     display(chum)
    Phagosight.FramStaticMap{I}=chum;
    Phagosight.FramStaticRatio(I)=mean(chum);
     clear chum
end;
disp('End Phagosight Frame and Track Static Ratio')
%% Calculate num of Hops - Andreized.  
%get the first node
%%note current $$ handles.distanceNetwork.numHops doesn't count first frame
distmap.timediscrepency=Arr_TimeDiscrepency;
distmap.distdiscrepency=Phagosight.distmap-distmap.Jun;
%arr_ddsum=sum(distmap.distdiscrepency); %the distance discrepency sum
arr_ddsum=Phagosight.GTs_arrdistpermp-GTs.arr_distpermp;
% Andre.lengthframeofhops; vs Andre.lengthframeofTrack
NN_uIDtrack; %Unique ID per Track derived from handles.nodeNetwork
NN_framesinTrack;
%%%%%%%%%%%%%%%%%%
numTracks=max(handles.finalLabel);
numHops=(0);
P_tat=0;J_tat=0; %don't need to because its essentiall reseting the array
P_totaldist=(0);J_totaldist=(0);
J_absVel=(0); P_absVel=(0);
dist_nodes=(0);
P_tortuosity=(0); J_tortuosity=(0);
P_meanderRatio=(0); J_meanderRatio=(0);
%%%%%%%%%%%%%% FY NODE NETWORK IS TWISTED UP [Y X Z] i.e 2 1 3
staticRatio=(0); 
nodeNetwork=handles.nodeNetwork(:,[1 2 5]);
nodeNetworkX = (nodeNetwork(:,2));
nodeNetworkY = (nodeNetwork(:,1));
for J=1:length(arr_GT)
    lensgt=arr_GT(J);
    GTsuID={0};
    gt_nodxy={0};
    gt_xyuIDwRF={0};
    in_Wound={0};
    in_WoundCount=0;
    frame_node={0};
    GT_hops={0};
    if arr_boo_GT(J)==1 
            for currTrack=1:numTracks
                ft=NN_framesinTrack{currTrack};
                if ft(1)<= lensgt
                    if J==1
                        temp_idx=find(ft<=lensgt);
                        temp_uID=NN_uIDtrack{currTrack}(temp_idx);
                        gt_node=temp_idx(end);
                        node1=NN_uIDtrack{currTrack}(1);%handles.finalNetwork(1,currTrack);  %start
                        node2=NN_uIDtrack{currTrack}(gt_node);
                    elseif J==length(arr_GT)
%                         orig_numHp=Andre.lengthframeofhops(currTrack); %updated
                        temp_uID=NN_uIDtrack{currTrack};
                        gt_node=NN_uIDtrack{currTrack}(end);
                        node1=NN_uIDtrack{currTrack}(1);%handles.finalNetwork(1,currTrack);  %start
                        node2=gt_node;%handles.finalNetwork((orig_numHp),currTrack);
                    else
                        temp_idx=find(ft<=lensgt);
                        temp_uID=NN_uIDtrack{currTrack}(temp_idx);
                        gt_node=temp_idx(end);
                        node1=NN_uIDtrack{currTrack}(1);%handles.finalNetwork(1,currTrack);  %start
                        node2=NN_uIDtrack{currTrack}(gt_node);
                    end;
%                          disp(temp_uID')
%                          pause() 
                    frame1=handles.nodeNetwork(node1,5);
                    frame2=handles.nodeNetwork(node2,5);
                    gotcha=diff([frame1 frame2]);
                    numHops(currTrack)=gotcha;            
               %%%%Method#1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Calculate Distance using handles based metrics%%%%%%%%%%
    %                 tat=handles.distMaps.absDistPerHop(1:end,currTrack); %handles.distMaps.absDistPerHop(1:numHops(currTrack),currTrack)
    %                 tat(tat==0)=[];
    %                 totalDist(currTrack)=sum(tat);
               %%%%Method#2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
               %Calculate Distance using Andreized Metrics%%%%%%%%%%%%%%
%                    P_tat=Phagosight.distmap(1:end,currTrack);
%                    P_tat(P_tat==0)=[];
%                    P_totaldist(currTrack)=P_tat;%so this is
%                    outdated..hence bellow
                   P_totaldist(currTrack)=Phagosight.DistPerMp(currTrack);
                   J_totaldist(currTrack)=GTs.Arr_distpermp{J}(currTrack); %This takes into consideration the GTs distances
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   P_absVel(currTrack)=P_totaldist(currTrack)/numHops(currTrack);%or Andre.lengthframeofTrack (Phagosight averages skip)
                   J_absVel(currTrack)=J_totaldist(currTrack)/(Andre.lengthframeofhops(currTrack)-2);%numHops(currTrack); (Jun Ignore skips)
               %%%tortosity and meanduring ratio%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    x_node1=handles.nodeNetwork(node1,1);
                    y_node1=handles.nodeNetwork(node1,2);
                    z_node1=handles.nodeNetwork(node1,3)*ZstepPixel-ZstepPixel;
                    x_node2=handles.nodeNetwork(node2,1);
                    y_node2=handles.nodeNetwork(node2,2);
                    z_node2=handles.nodeNetwork(node2,3)*ZstepPixel-ZstepPixel;
                    dist_nodes(currTrack)=pdist([x_node1,y_node1,z_node1; x_node2,y_node2,z_node2],'Euclidean');
                    P_dist_nodes(currTrack)=dist_nodes(currTrack);
                    j_checkdistnodes=dist_nodes(currTrack)-arr_ddsum(currTrack,J);
                    if j_checkdistnodes>0
                        j_check=j_checkdistnodes;
                    else
                        j_check=1;
                    end;
                    J_dist_nodes(currTrack)=j_check;
                    P_tortuosity(currTrack)=P_totaldist(currTrack)/dist_nodes(currTrack);
                    P_meanderRatio(currTrack)=1/P_tortuosity(currTrack); 
                    J_tortuosity(currTrack)=J_totaldist(currTrack)/J_dist_nodes(currTrack);
                    J_meanderRatio(currTrack)=1/J_tortuosity(currTrack);
                %%%%%%%%%%Static Ratio%%%%%%%%%%%%%%%%%%%%%%%%
                    % mtotdist=mean(arr_totDist);
                    % stdtotdist=std(arr_totDist);
                    fudgefactor=1/1;%-stdtotdist;
                    Mav=mean(distmap.JunVelPerMp);
                    Std=std(distmap.JunVelPerMp);
                    if fudgefactor<0
                        fudgefactor=1/1;%(pixel/frame) i.e. (1.55/1.9)
                    else
                    end;  
%                     num_staticMps=length(find(handles.distMaps.absDistPerHop(1:numHops(currTrack),currTrack)<=fudgefactor));
%                     staticRatio(currTrack)=num_staticMps/numHops(currTrack);                
                    J_tat=GTs.Arr_distmap{J}{currTrack}; %Took a distmap at that time GT only
                    J_tat(J_tat==0)=[];
                    lens_jtat=length(J_tat);
                    num_staticMps=length(find(J_tat<=fudgefactor));
                    staticRatio(currTrack)=num_staticMps/lens_jtat;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% In determine in Wound or Nt
                clear gtnodxy
                GTsuID{currTrack}=temp_uID;
                lens_uId=length(temp_uID);
%                 figure(currTrack);
%                     title(strcat('Track : ',num2str(currTrack)));
% %                     imshow(woundRegion)
 %                    imagesc(valsPH2{1}{(handles.numFrames-2)});colormap('pink');
%                 hold on
                inWound=(0);
                framenode=(0);
                for ii=1:lens_uId
                    idx=temp_uID(ii);
                    i_nodx=handles.nodeNetwork(idx,1);
                    i_nody=handles.nodeNetwork(idx,2);
                    txt1 = strcat('\rightarrow ',num2str(ii));
%                         text(i_nody,i_nodx,txt1)
%                         plot(i_nody,i_nodx,'.','MarkerSize',20);     
                    gtnodxy(ii,1)=i_nodx;
                    gtnodxy(ii,2)=i_nody;
                    switch isnan(i_nodx)
                        case{1}
                            inWound(ii)=NaN;
                        otherwise
                            inWound(ii)=woundRegion(round(i_nodx),round(i_nody));
                    end;

                end;                               
                framenode=handles.nodeNetwork(temp_uID,5);
                in_Wound{currTrack}=inWound';
                in_WoundCount(currTrack)=sum(inWound)>0; %just tells us if its wound
                frame_node{currTrack}=framenode;
                gt_nodxy{currTrack}=gtnodxy;
                %%%%%%%%%%%%REDO HOP DISTANCES
                %%I'm finding a discorrelation between distmap and
                %%handles.nodeNetwork hops
                %%This section is taking every ID from the track and
                %%listing its hops---hopefully we are meant to ignore the
                %%ones that distmaps is ignoring.  If not ill need to redo
                %%velocity. 
                GThops=(0);
                for ii=2:lens_uId
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    idx1=temp_uID(ii-1);
                    i_nodx1=handles.nodeNetwork(idx1,1);
                    i_nody1=handles.nodeNetwork(idx1,2);
                    i_nodz1=handles.nodeNetwork(idx1,3)*ZstepPixel-ZstepPixel;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    idx=temp_uID(ii);
                    i_nodx=handles.nodeNetwork(idx,1);
                    i_nody=handles.nodeNetwork(idx,2);
                    i_nodz=handles.nodeNetwork(idx,3)*ZstepPixel-ZstepPixel;
                    pgthops=pdist([i_nodx1,i_nody1,i_nodz1; i_nodx,i_nody,i_nodz],'Euclidean');
                    GThops(ii)=pgthops;
                end;
                GT_hops{currTrack}=GThops';                   
                xyuIDwRF=horzcat(temp_uID, gtnodxy(:,1),gtnodxy(:,2),GThops',inWound',framenode); 
                gt_xyuIDwRF{currTrack}=xyuIDwRF;
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                A=distmap.Jun;
                A(A==0)=nan;
                AB=nanmean(A);
                ABC=mean(AB); %gives the average hop ---but actually velocity is more accurate
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    GT_hops{currTrack}=nan;
                    frame_node{currTrack}=nan;
                    in_Wound{currTrack}=nan;
                    in_WoundCount(currTrack)=nan;
                    gt_xyuIDwRF{currTrack}=nan;
                    gt_nodxy{currTrack}=nan;
                    GTsuID{currTrack}=nan;
                    staticRatio(currTrack)=nan;
                    numHops(currTrack)=nan;  
                    P_totaldist(currTrack)=nan;
                    J_totaldist(currTrack)=nan;
               
                    P_absVel(currTrack)=nan;
                    J_absVel(currTrack)=nan;
                    dist_nodes(currTrack)=nan;
                    P_dist_nodes(currTrack)=nan;
                    J_dist_nodes(currTrack)=nan;
                    P_tortuosity(currTrack)=nan;
                    P_meanderRatio(currTrack)=nan;
                    J_tortuosity(currTrack)=nan;
                    J_meanderRatio(currTrack)=nan;
                end;      
            end;
           %%%%%%%%%%%%%%%%Static Ratio By Frame%%%%%%%%%%%%%%%%%%%%%%%%
           %%if you wanted you can do velocity by frame or dist by frame
            GTdistmap=GTs.distmap.distmap{J};
            [lensx,lensy]=size(GTdistmap);    
            mpFrame=0;
            frame_staticRatio=0;
            cf=0;
            for currFrame=2:lensx
                cf=cf+1;
                mpFrame=GTdistmap(currFrame,:);
                mpFrame(mpFrame==0)=[];
                lens_mpFrame=length(mpFrame);
                num_framstaticMps=length(find(mpFrame<=fudgefactor));
                frame_staticRatio(cf)=num_framstaticMps/lens_mpFrame;
            end;
    else 
        frame_staticRatio=nan;
        staticRatio=nan;
        J_tortuosity=nan;
        J_meanderRatio=nan;
        J_totaldist=nan;
    end;    
    GTs.hops{J}=GT_hops;
    GTs.frameNode{J}=frame_node;
    GTs.inWound{J}=in_Wound;
    GTs.inWoundCount{J}=in_WoundCount;
    GTs.xyuIDwR{J}=gt_xyuIDwRF;
    GTs.nodeXY{J}=gt_nodxy;
    GTs.uID{J}=GTsuID;
    GTs.frame_staticRatio{J}=frame_staticRatio;
    GTs.staticRatio{J}=staticRatio;
    GTs.tortuosity{J}=J_tortuosity;
    GTs.meanderRatio{J}=J_meanderRatio;
    GTs.totaldist{J}=J_totaldist;
end;
Phagosight.tortuosity=P_tortuosity; 
%     FigHandle = figure;
%     set(FigHandle, 'Position', [100, 100, 350, 250]); %[x1, y1, x2, y2]
%      histogram(Phagosight.tortuosity,Num_TrackMax);
%      ylabel('Phagosight Tortuosity')
%     title(strcat('Histogram of Tortuosity  30-',str_ArrGT{3}))
Phagosight.meander=P_meanderRatio;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [100, 100, 350, 250]); %[x1, y1, x2, y2]
%      histogram(Phagosight.meander,50);
%      ylabel('Phagosight Meandering Index')
%     title(strcat('Histogram of Meandering Index  30-',str_ArrGT{3}))
distmap.JunTortuosity=J_tortuosity;
distmap.JunMeander=J_meanderRatio;
distmap.JunframeStaticRatio=frame_staticRatio;
distmap.JunStaticRatio=staticRatio;
disp('End tortuosity, meanderRatio, & Static Ratio')

%% Paredes for excel with Golden Times Cummulative 
%%%Missing Intervals of Golden Times (non cummulative)
clear Paredes
Paredes.totalDist=J_totaldist;
Paredes.absVel=distmap.JunVelPerMp;
Paredes.tortuosity=J_tortuosity;
Paredes.meanderRatio=J_meanderRatio;
Paredes.staticRatio=staticRatio;    
Paredes.frame_staticRatio=frame_staticRatio;
Paredes.twist=sum(J_totaldist)/sum(dist_nodes);
Paredes.mr=sum(dist_nodes)/sum(J_totaldist);
Tweeker=[Paredes.twist;Paredes.mr];
disp('End Paredes Metrics and Tweeker to Check')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New with GT things to save In Excel Sheet%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tasks 1
Dr.totalDist=GTs.Arr_distpermp; %horizontal
Dr.absVel=GTs.Arr_velpermp; %vertical
Dr.staticRatio=GTs.staticRatio;
Dr.frame_staticRatio=GTs.frame_staticRatio;
Dr.tortuosity=GTs.tortuosity;
Dr.meanderRatio=GTs.meanderRatio;
disp('End Script Zebrafish_TotalDoctorIt')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Search for Time Placements of Tracks in each GT
if arr_boo_GT(3)
    f0tof120=Andre.distmap(1:GT_120,:);
    distmapf0tof120=GTs.distmap.distmap{3};
    [~,mp]=size(distmapf0tof120);   %GT_120;
else
    GT_120=frames-1;
    arr_boo_GT(3)=1;
    disp('Ending 120 minute at max frames')
    [~,mp]=size(distmap.Jun);   %GT_120;
end;
%%%%%%%%%%%%%%%%
% Determining Time POints
clear time
time.i=(0); %Intial Time point of Mp
time.f=(0); % Final Time point of Mp
time.idx_GT60=(0); % idx of mp in nodeNetwork  at time point 60
time.idx_GT90=(0); % idx of mp in nodeNetwork  at time point 90
time.idx_GT120=(0);% idx of mp in nodeNetwork  at time point 120
time.idx_GT45=(0);% idx of mp in nodeNetwork  at time point 120
time.GT60=(0); % idx of mp in nodeNetwork  at time point 60
time.GT90=(0); % idx of mp in nodeNetwork  at time point 90
time.GT120=(0);% idx of mp in nodeNetwork  at time point 120
time.GT45=(0);% idx of mp in nodeNetwork  at time point 120

for i=1:mp
    [idx_mp]=find(Col14==i);
    time.i(i)=Col05(idx_mp(1));%NOTE:Col05 is Time Frame   
    time.f(i)=Col05(idx_mp(end));
    time.idx_i(i)=(idx_mp(1));
    time.idx_f(i)=(idx_mp(end));    
    if arr_boo_GT(3)==1
        [idx_mpx]=find(Col14==i);
        [idx_whaa45]=find(Col14==i & Col05==GT_45);%NOTE:Col14 Track
        [idx_whaa60]=find(Col14==i & Col05==GT_60);
        [idx_whaa90]=find(Col14==i & Col05==GT_90);
        [idx_whaa120]=find(Col14==i & Col05==GT_120);
        boo_check45=time.i(i)<GT_45 && GT_45<=time.f(i);
        boo_check60=time.i(i)<GT_60 && GT_60<=time.f(i);
        boo_check90=time.i(i)<GT_90 && GT_90<=time.f(i);
%         boo_check120=time.i(i)<GT_120 && GT_120<=time.f(i); %5/11 update
        boo_check120=time.i(i)<GT_120;
        if isempty(idx_whaa45)
             sup45=nan;
                if boo_check45
                    Jun=distmap.Jun(1:GT_45,i);
                   	[non_zero]=find(Jun~=0);
                    tmp=abs(non_zero-GT_45);
                    [idx idx] = min(tmp); % index of closest value
                    closest = non_zero(idx); % closest value
                    [idx_whaa45]=find(Col14==i & Col05==closest);
                     sup45=Col05(idx_whaa45);
                else
                    idx_whaa45=nan;
                end;             
        else 
            sup45=Col05(idx_whaa45);
        end;  
       
        if isempty(idx_whaa60)
                sup60=nan;
                if boo_check60
                    Jun=distmap.Jun(1:GT_60,i);
                   	[non_zero]=find(Jun~=0);
                    tmp=abs(non_zero-GT_60);
                    [idx, idx] = min(tmp); % index of closest value
                    closest = non_zero(idx); % closest value
                    [idx_whaa60]=find(Col14==i & Col05==closest);
                    sup60=Col05(idx_whaa60);
                else
                    idx_whaa60=nan;
                end;             
        else
            sup60=Col05(idx_whaa60);
        end;
        if isempty(idx_whaa90)
            sup90=nan;
                if boo_check90
                    Jun=distmap.Jun(1:GT_90,i);
                   	[non_zero]=find(Jun~=0);
                    tmp=abs(non_zero-GT_90);
                    [idx idx] = min(tmp); % index of closest value
                    closest = non_zero(idx); % closest value
                    [idx_whaa90]=find(Col14==i & Col05==closest);
                    sup90=Col05(idx_whaa90);
                else
                    idx_whaa90=nan;
                end;
        else       
             sup90=Col05(idx_whaa90);
        end;
        
        if isempty(idx_whaa120)
            sup120=nan;
                if boo_check120
                    if time.f(i)>=GT_90
                        Jun             =   distmap.Jun(1:GT_120,i);
                        [non_zero]  	=   find(Jun~=0);
                        tmp             =   abs(non_zero-GT_120);
                        [idx idx]       =   min(tmp); % index of closest value
                        closest         =   non_zero(idx); % closest value
                        [idx_whaa120]   =   find(Col14==i & Col05==closest);
                        sup120          =   Col05(idx_whaa120);
                    else
                        idx_whaa120     =   nan;
                    end;
                else
                    idx_whaa120=nan;
                end;
        else
            sup120=Col05(idx_whaa120);
        end;
        
    time.GT45(i)        =   sup45;%%problem is that it will only capture those that start at 45.  
    time.GT60(i)        =   sup60;
    time.GT90(i)        =   sup90;
    time.GT120(i)       =   sup120;
    time.idx_GT45(i)    =   idx_whaa45;%%problem is that it will only capture those that start at 45.  
    time.idx_GT60(i)    =   idx_whaa60;
    time.idx_GT90(i)    =   idx_whaa90;
    time.idx_GT120(i)   =   idx_whaa120;
    else
        disp('Imaging Does not go Beyond up to 120 minutes after injury')
    end;  
end;
time.GT         ={time.GT60,time.GT90,time.GT120};
time.idx_GT     ={time.idx_GT60,time.idx_GT90,time.idx_GT120};
disp('End: "time" structure array created')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%% Total Distance for Tortosity
% Note to ADP - Inefficently Written.  
clear SM Dizzle Zirmi
str_ArrGT={'60min','90min'...
    ,'120min'};
% Iteration 30 min Intervals
Dizzle.distmap.f45to60=distmap.Jun(GT_45:GT_60,:);
Dizzle.distmap.f45to90=distmap.Jun(GT_45:GT_90,:);
Dizzle.distmap.f45to120=distmap.Jun(GT_45:GT_120,:);
Dizzle.distmap.dm30cum={Dizzle.distmap.f45to60...
    ,Dizzle.distmap.f45to90,Dizzle.distmap.f45to120};%distmap30minIntervals
% Cummulative 30 min
Dizzle.distmap.f45to60=distmap.Jun(GT_45:GT_60,:);
Dizzle.distmap.f60to90=distmap.Jun(GT_60:GT_90,:);
Dizzle.distmap.f90to120=distmap.Jun(GT_90:GT_120,:);
Dizzle.distmap.dm30={Dizzle.distmap.f45to60...
    ,Dizzle.distmap.f45to90,Dizzle.distmap.f45to120};
%%Selected Metrics - Tortuosity, Meandering Ratio, 
% DISTMAP= Andre.distmap
% Iteration 30 min Intervals
SM.distmap.f45to60      =   Andre.distmap(GT_45:GT_60,:);
SM.distmap.f60to90      =   Andre.distmap(GT_60:GT_90,:);
SM.distmap.f90to120     =   Andre.distmap(GT_90:GT_120,:);
SM.distmap.dm30         =   {SM.distmap.f45to60...
    ,SM.distmap.f60to90,SM.distmap.f90to120};%distmap30minIntervals
% Cummulative 30 min
SM.distmap.f45to60      =   Andre.distmap(GT_45:GT_60,:);
SM.distmap.f45to90      =   Andre.distmap(GT_45:GT_90,:);
SM.distmap.f45to120     =   Andre.distmap(GT_45:GT_120,:);
SM.distmap.dm30cum      =   {SM.distmap.f45to60,...
                                SM.distmap.f45to90,SM.distmap.f45to120};
%% Selected Macrophages
%--First Select Macrophage tracks of interest (within a HIGH% of the
%frames)
[~,mp]      =   size(distmap.Jun);
frames      =   handles.numFrames;
 Selected   =   (0); % Selected Metrics
 checkSM    =   (0);
 a          =   0;
 %MajorityTracksPercent=.70; %<--HIGH% selection - ADP
for i=1:mp
    a           =   a+1;
    Awo         =   distmap.Jun(:,i);
    Awo         =   SM.distmap.f45to120(:,i);
    Amo         =   Awo(Awo~=0);
    lens_track  =   length(Amo);
    det_lens    =   lens_track/GT_120;
    if det_lens > ParameterA %MajorityTracksPercent 
        Selected(a)     =   i;
    else
    end;
    checkSM(i)  =   det_lens;
end;
Selected    =   Selected(Selected~=0); %Provides the Tracks that fall within
display('END: SM-Selected Macrophages of Interested')
%% --Second Do it Per 30 minute intervals 
% Parameter_gtA=0.9; %90 percent
gt_SM       =   {0};
for J=1:3%length(arr_GT) %<-- This is the list of frames NOT time points   
    if arr_boo_GT(J)==1
        if J==1
            Bwo     =   SM.distmap.f45to60;
        elseif J==2
            Bwo     =   SM.distmap.f60to90;
        elseif J==3
            Bwo     =   SM.distmap.f90to120;
        end;
    else 
        Bwo     =   SM.distmap.f90to120;
    end;
    [x,y]   =   size(Bwo);
    gtsm    =   (0);
    gtSM    =   0;
    b       =   0;
    Bmo     =   0;
    for i=1:y
        b           =   b+1;
        Cwo         =   Bwo(:,i);    
        Bmo         =   Cwo(Cwo~=0);
        lens_track  =   length(Bmo);
        det_lens    =   lens_track/x;
        if det_lens > Parameter_gtA %Percent within Tracks 
            gtSM(b)         =   i;
        else
        end;
        gtcheckSM(i)=   det_lens;       
    end;
    gtSM                =   gtSM(gtSM~=0);
    gtSelectedMetrics   =   [gtSM,Selected];
    gtSMcombo           =   unique(gtSelectedMetrics);
    %gt_SM{J}           =   gtSM;
    gt_SM{J}            =   gtSMcombo;

end;
%%%%%%%%%%%%%%%%%%-updated 05/11/2017 add Selected to gt_SM since selected
%%%%%%%%%%%%%%%%%%are the values of importance here.  Need to make sure
%%%%%%%%%%%%%%%%%%they get selected.  
disp('END: Finished Calculating Selected Macrophages of Interest')
%% Section 19 - Selected Tortuosity and Meandering, Velocity  && Velocity Per Frame

clear SM.tortuosity SM.meandering SM.staticratio
%--Third Apply this knowledge to getting the Data
%-Structure Array
clear Meandering Tortuosity Velocity StaticRatio
Thismap             = 0;
Thismapcum          = 0;
%-Cell Arrays
d_interval          ={0};
d_cum               ={0};
l_nodes_interval    ={0};
l_nodes_cum         ={0};
%--Prep Static Ratio
Thismap             = distmap.Jun;  
Nanmap              = Thismap;
Nanmap(Nanmap==0)   = nan;
Staticmap           = Nanmap;
[x,y]               = size(Nanmap);
% StaticLimit         = 1/ti_d; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%- ADP
%-Determin the Static Map of the nanmap from distmap.Jun
for i=1:x
    for j=1:y
        if isnan(Nanmap(i,j))
        else
            Staticmap(i,j)=Nanmap(i,j)>=ParameterB;%(1/micronsPerPixel); %micronsPerPixel=1.66 on confocal
        end;
    end;
end;
disp('END: Static Limit')
%--Selected Tortuosity and  Meandering, Velocity  && Velocity Per Frame

for I=1:3 % per GT_60, GT_90, and GT_120
    v_perframe_interval     = 0;
    v_perframe_cum          = 0;
    s_perframe_interval     = 0;
    s_perframe_cum          = 0;
    S_iterations            =(0);
    S_iterations            =[gt_SM{I} Selected];
    S_iterations            = unique(S_iterations);
     if arr_boo_GT(I) %ensures that these frames exist
        for i=S_iterations  %Error gt_SM doesn't include SM, need to combine the 2.  
       
            frame_GT=time.GT{I}(i);
            idx_GT=time.idx_GT{I}(i);
            if isnan(idx_GT) % If Nan then this track can't be done
                warning(strcat('Zirmi Metrics Part1:GT#',num2str(I),'-Selected Track#',num2str(i),' is nan'))
%                 pause();
                Velocity.interval{I}(i)     =nan; %Jun - all zeros discounted
                Velocity.cum{I}(i)          =nan; % Jun -all zeros discounted    
                Tortuosity.interval{I}(i)   =nan;
                Tortuosity.cum{I}(i)        =nan;
                Meandering.interval{I}(i)   =nan;
                Meandering.cum{I}(i)        =nan;
                StaticRatio.interval{I}(i)  =nan;
                StaticRatio.cum{I}(i)       =nan;
  
            else 
            %%%#1-Identifying the node of interest               
                nodei=time.idx_i(i);
                framei=time.i(i);
                %not nan
                  if I==1; %if I=1 therefore up to 60minutes
                      node1=time.idx_i(i); %select initial Point (not 45minutes)
                      frame1=time.i(i);
                      coords           =[nan,nan,nan];
                  else
                      node1=time.idx_GT{I-1}(i); %if not start then select the  30 min interval before current
                      frame1=time.GT{I-1}(i);
                  end;
                node2=time.idx_GT{I}(i); %this will the 30 min interval of now
                frame2=time.GT{I}(i);%this will the 30 min interval of now
            %%%#2-Calculate 'd' the sum of the distances
                %-to account for error if starts not in first 30 minute
                %window
                if isnan(frame1)
                    frame1=framei;
                    node1=nodei;
                else
                end;
                d_interval{I}(i)=sum(Andre.distmap(frame1:frame2,i));%/(frame2-frame1); % Andre- all zeros averaged
                d_cum{I}(i)=sum(Andre.distmap(framei:frame2,i));%/(frame2-framei); % Andre- all zeros averaged
                vi=(Thismap(frame1:frame2,i));
                vc=(Thismap(framei:frame2,i));
                vi=vi(vi~=0);
                vc=vc(vc~=0);
                Velocity.interval{I}(i)=mean(vi); %Jun - all zeros discounted
                Velocity.cum{I}(i)=mean(vc); % Jun -all zeros discounted              
            %%%#3-Vector Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                z_nodei                     =handles.nodeNetwork(nodei,3)*ZstepPixel-ZstepPixel;
                x_nodei                     =handles.nodeNetwork(nodei,1);
                y_nodei                     =handles.nodeNetwork(nodei,2);
                %--
                z_node1                     =handles.nodeNetwork(node1,3)*ZstepPixel-ZstepPixel;
                x_node1                     =handles.nodeNetwork(node1,1);
                y_node1                     =handles.nodeNetwork(node1,2); 
                %--
                x_node2                     =handles.nodeNetwork(node2,1);
                y_node2                     =handles.nodeNetwork(node2,2);
                z_node2                     =handles.nodeNetwork(node2,3)*ZstepPixel-ZstepPixel;
            %%%#4 Vector Distances 
                l_nodes_interval{I}(i)      =pdist([x_node1,y_node1,z_node1; x_node2,y_node2,z_node2],'Euclidean');
                l_nodes_cum{I}(i)           =pdist([x_nodei,y_nodei,z_nodei; x_node2,y_node2,z_node2],'Euclidean');
            %%%#5-tortosity and meanduring ratio%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                Tortuosity.interval{I}(i)   =d_interval{I}(i)/l_nodes_interval{I}(i);
                Tortuosity.cum{I}(i)        =d_cum{I}(i)/l_nodes_cum{I}(i);
                Meandering.interval{I}(i)   =l_nodes_interval{I}(i)/d_interval{I}(i);
                Meandering.cum{I}(i)        =l_nodes_cum{I}(i)/d_cum{I}(i);
            %%%#6-Static Ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %-interval
                srci                        =Staticmap(frame1:frame2,i);
                arr_srci                    =srci(~isnan(srci));
                sr_total_interval           =length(arr_srci); %total length (number of things)
                sr_count_interval           =sum(arr_srci); % only number that are static
                StaticRatio.interval{I}(i)  =(sr_total_interval-sr_count_interval)/sr_total_interval;
                %-cummulative
                srcc                        =Staticmap(framei:frame2,i);
                arr_srcc                    =srcc(~isnan(srcc));
                sr_total_cum                =length(arr_srcc);
                sr_count_cum                =sum(arr_srcc);
                StaticRatio.cum{I}(i)       =(sr_total_cum-sr_count_cum)/(sr_total_cum);
%                 pause()
            end;    
        end;
        Coordinates{I}  = coords;
        %--Per Frames
        framei=1;
        if I==1
            frame1          = 1;
        else
            frame1          = arr_GT(I-1);
        end;
        frame2          = arr_GT(I);
        [col,row]       = size(Thismap);
        if frame2>col
            frame2 = col;
        else
        end;
        %-Velocity per Frame
        v_perframe_interval=Thismap(frame1:frame2,:);
        v_perframe_cum=Thismap(framei:frame2,:);    
        Velocity.frame.interval{I}=v_perframe_interval(v_perframe_interval~=0);
        Velocity.frame.cum{I}=v_perframe_cum(v_perframe_cum~=0); 
        %-Static per Frame
        s_perframe_interval=Staticmap(frame1:frame2,:);
        s_perframe_cum=Staticmap(framei:frame2,:);    
        StaticRatio.frame.interval{I}=s_perframe_interval(~isnan(s_perframe_interval));
        StaticRatio.frame.cum{I}=s_perframe_cum(~isnan(s_perframe_cum)); 
        %Note: Its the length(StaticRa...-sum(StaticRa...)/length(StaticRa...)
    else % for the time points that the frames do not exist
        Velocity.interval{I}=nan;
        Velocity.cum{I}=nan;
        Tortuosity.interval{I}=nan;
        Tortuosity.cum{I}=nan;
        Tortuosity.phagosight=handles.distanceNetwork.tortuosity;
        Meandering.interval{I}(i)=nan;
        Meandering.cum{I}(i)=nan;
        Meandering.phagosight=handles.distanceNetwork.meanderRatio;
        StaticRatio.interval{I}(i)=nan;
        StaticRatio.cum{I}(i)=nan;
        StaticRatio.phagosight=handles.distanceNetwork.staticRatio;
    end;   
    SM.velocity.interval{I}=Velocity.interval{I}(gt_SM{I});
    SM.velocity.intervalsm{I}=Velocity.interval{I}(Selected);
    SM.velocity.cum{I}=Velocity.cum{I}(Selected); % -ERRORED 4/27/2017 Index exceeds matrix dimensions.
    SM.tortuosity.interval{I}=Tortuosity.interval{I}(gt_SM{I});
    SM.tortuosity.intervalsm{I}=Tortuosity.interval{I}(Selected);
    SM.tortuosity.cum{I}=Tortuosity.cum{I}(Selected); %maybe not incorrect
    SM.meandering.interval{I}=Meandering.interval{I}(gt_SM{I});
    SM.meandering.intervalsm{I}=Meandering.interval{I}(Selected);
    SM.meandering.cum{I}=Meandering.cum{I}(Selected); %maybe not incorrect
    SM.staticratio.interval{I}=StaticRatio.interval{I}(gt_SM{I});
    SM.staticratio.intervalsm{I}=StaticRatio.interval{I}(Selected);
    SM.staticratio.cum{I}=StaticRatio.cum{I}(Selected); %maybe not incorrect
    
    Zirmi.velocity.interval{I}=Velocity.interval{I};
    Zirmi.velocity.intervalsm{I}=Velocity.interval{I};
    Zirmi.velocity.cum{I}=Velocity.cum{I}; % -ERRORED 4/27/2017 Index exceeds matrix dimensions.
    Zirmi.tortuosity.interval{I}=Tortuosity.interval{I};
    Zirmi.tortuosity.intervalsm{I}=Tortuosity.interval{I};
    Zirmi.tortuosity.cum{I}=Tortuosity.cum{I}(Selected); %maybe not incorrect
    Zirmi.meandering.interval{I}=Meandering.interval{I};
    Zirmi.meandering.intervalsm{I}=Meandering.interval{I};
    Zirmi.meandering.cum{I}=Meandering.cum{I}(Selected); %maybe not incorrect
    Zirmi.staticratio.interval{I}=StaticRatio.interval{I};
    Zirmi.staticratio.intervalsm{I}=StaticRatio.interval{I};
    Zirmi.staticratio.cum{I}=StaticRatio.cum{I}; %maybe not incorrect
end;
SM.velocity.phagosight=handles.distanceNetwork.absVelocity(Selected); %correct
SM.tortuosity.phagosight=handles.distanceNetwork.tortuosity(Selected); %correct
SM.meandering.phagosight=handles.distanceNetwork.meanderRatio(Selected); % correct
SM.staticratio.phagosight=handles.distanceNetwork.staticRatio(Selected);%correct
% plotTracks(handles,1,find(handles.distanceNetwork.numHops>50));%%
% plotTrackStats(handles,2,15)
disp('END: Selected Metrics and gt_SM'); 
%% Section 20 - Forward and Backward Determinations
%-First determine centroid of Wound center.  
clear dd
dd                          =  regionprops(wR1); % This is the center of the wound gap
y_nodeW                     =  dd.Centroid(1); %note the swap of (1)=x
x_nodeW                     =  dd.Centroid(2); %note the swap of (2)=y
ee                          =  regionprops(wR2); %This is the Center of Notochord
y_nodeN                     =  ee.Centroid(1); %note the swap of (1)=x
x_nodeN                     =  ee.Centroid(2); %note the swap of (2)=y
 %Distance in pixels between Wound Gap epicenter and Notochard marked manually
dist_notochord=(pdist([x_nodeW,y_nodeW; x_nodeN,y_nodeN],'Euclidean')*coe); 
% % % % figure(); 
% % % % imagesc(wR1);colormap(gray); hold on;
% % % % plot(dd.Centroid(1),dd.Centroid(2),'or',...
% % % %     'LineWidth',50);hold on; 
% % % % %-To get hPlotNet variable you need to run 'ZebrafishTracks_2D.m' first 
% % % % %-This plots the (1) track on the same plane as wR1
% % % % numTrack=5;
% % % % ZZ=ZZ*0+2; 
% % % %  hPlotNet.handlePlot(numTrack)=plot3(...
% % % %     YY,XX,ZZ,'LineStyle',':','marker',colorID1{1+rem(counterTrack-1,11)},... %Note the switch between YY and XX to plot on in same references, but is the same true for measuring????
% % % %     'color',colorID3(1+rem(counterTrack-1,20),:),'markersize',...
% % % %     4,'linewidth',4); hold on;
% % % % txt1=strcat('\leftarrow ',num2str(numTrack));
% % % % text(YY(1),XX(1),txt1,'Color','cyan','FontSize',12)
%--Loop
%Clear Structure Arrays
clear Forward Backward FtoB FBratio WoundScoreUm
clear WoundScore1234 WoundStartUm WoundDist
NodeDistDiff_cum            ={0};
NodeDistDiff_interval       ={0};

 for numGT = 1:3 % per GT_60, GT_90, and GT_120
     %Solved: Error gt_SM doesn't include SM, need to combine the 2.  
     S_iterations                =(0);
     S_iterations                =[gt_SM{numGT} Selected];
     S_iterations                = unique(S_iterations);
     wound_dist_nodes_interval   ={0}; 
     wound_dist_nodes_cum        ={0};
     distdiff_cum                ={0};
     distdiff_interval           ={0};
     forwardratio_interval       =(0);
     forwardratio_cum            =(0);
     backwardratio_interval      =(0);
     backwardratio_cum           =(0);
     ftobratio_interval          =(0);
     ftobratio_cum               =(0); 
     wound_um_interval           =(0);
     wound_um_cum                =(0);
     startwound_um_interval      =(0);
     startwound_um_cum           =(0);
     wound1234_cum               =(0);
     wound1234_interval          =(0);
     if arr_boo_GT(numGT) %ensures that these frames exist
         
        for numTrack=S_iterations %This specific 'Selected group' TRACK 
            frame_GT        = time.GT{numGT}(numTrack);
            idx_GT          = time.idx_GT{numGT}(numTrack);
            idx_initial     = time.idx_i(numTrack); %-alternatively time.idx_GT45 %%%%%%%This right now is starting from the literal initial point NoT 45. 
            
            if isnan(idx_GT) % If Nan then this track can't be done
            %Since we are doing gt_SM and Selected some gt_SM are not in
            %early GTs but in late GTs and/or vice versa
                warning(strcat('Zirmi Metrics Part2:GT#',num2str(numGT),'-Selected Track#',num2str(numTrack),' is nan'))
                distdiff_cum{numTrack}              = nan;
                forwardratio_cum(numTrack)          = nan;  
                backwardratio_cum(numTrack)         = nan;
                ftobratio_cum(numTrack)             = nan;
                distdiff_interval{numTrack}         = nan;
                forwardratio_interval(numTrack)     = nan;  
                backwardratio_interval(numTrack)    = nan;
                ftobratio_interval(numTrack)        = nan;
                wound_dist_nodes_cum{numTrack}      = nan;
                wound_dist_nodes_interval{numTrack} = nan;
                wound_um_interval(numTrack)         = nan;
                wound_um_cum(numTrack)              = nan;
                startwound_um_interval(numTrack)    = nan;
                startwound_um_cum(numTrack)         = nan;
                wound1234_cum(numTrack)             = nan;
                wound1234_interval(numTrack)        = nan;
            elseif isnan(idx_initial)
            %This would happen only when we do 45 because it might not exist    
                distdiff_cum{numTrack}              = nan;
                forwardratio_cum(numTrack)          = nan;  
                backwardratio_cum(numTrack)         = nan;
                ftobratio_cum(numTrack)             = nan; 
                distdiff_interval{numTrack}         = nan;
                forwardratio_interval(numTrack)     = nan;  
                backwardratio_interval(numTrack)    = nan;
                ftobratio_interval(numTrack)        = nan;
                wound_dist_nodes_cum{numTrack}      = nan;
                wound_dist_nodes_interval{numTrack} = nan;
                wound_um_interval(numTrack)         = nan;
                wound_um_cum(numTrack)              = nan;
                startwound_um_interval(numTrack)    = nan;
                startwound_um_cum(numTrack)         = nan;
                wound1234_cum(numTrack)                       = nan;
                wound1234_interval(numTrack)                  = nan;                
            else
                track_idxs                  = handles.finalNetwork(:,numTrack);
                track_idxs(track_idxs==0)   =[];
                upperlimit_idx              = find(track_idxs==idx_GT); %This tells where to go up to
                %--distinguish idx lowerlimit for interval\alternatively---
                %--cummulative will have initiallimit_idx for deriving values
                if numGT==1;
                    lowerlimit_idx                      = find(track_idxs==idx_initial);
                else
                    lowerlimit_idx                      = find(track_idxs==time.idx_GT{numGT-1}(numTrack));
                    %SOLVED: Error-what if its only present from current frame on!
                    %its only empty if it can't find that value
                    %if its empty its because its Nan in previous GTs
                    %therefore it starts on current frame
                    if isempty(lowerlimit_idx)
                        lowerlimit_idx                      = find(track_idxs==idx_initial);
                    else
                    end;
                end;                        
                %--for interval
                withzeroeswounddistnodesinterval=0;
                for i=lowerlimit_idx:upperlimit_idx                  
                    x_node1                             =HNN{1}(track_idxs(i),1);
                    y_node1                             =HNN{1}(track_idxs(i),2);
                    %as is now there are zeroes all up in these arrays.
                    withzeroeswounddistnodesinterval(i)=pdist([x_node1,y_node1; x_nodeW,y_nodeW],'Euclidean');   
                end;  
                wound_dist_nodes_interval{numTrack}=withzeroeswounddistnodesinterval(withzeroeswounddistnodesinterval~=0);%By doing this we lose the idx of nodes but shouldn't matter for them.
                %(0)=static|(+)=backward|(-)=forward
                distdiff_interval{numTrack}     = diff(wound_dist_nodes_interval{numTrack});
                lens_interval                   = length(distdiff_interval{numTrack} ); %Total 
                positive_backward_interval      = distdiff_interval{numTrack}(distdiff_interval{numTrack}>0); %greater than 0
                negative_forward_interval       = distdiff_interval{numTrack}(distdiff_interval{numTrack}<0); %less than than 0
                forwardratio_interval(numTrack) = length(negative_forward_interval)/lens_interval;  
                backwardratio_interval(numTrack)= length(positive_backward_interval)/lens_interval;
                ftobratio_interval(numTrack)    = length(negative_forward_interval)/length(positive_backward_interval);
                wound_um_interval(numTrack)     =(diff([wound_dist_nodes_interval{numTrack}(1) wound_dist_nodes_interval{numTrack}(end)]))*coe; %This is in um already NOTE: diff([100 40]) is -60.
                startwound_um_interval(numTrack)=wound_dist_nodes_interval{numTrack}(1);
                %--To Determine The Wound Score 1 2 3 4 by 150um
                if startwound_um_interval(numTrack) == 0 || isnan(startwound_um_interval(numTrack));
                    wound1234_interval(numTrack)              = nan;
                elseif startwound_um_interval(numTrack)<=dist_notochord %if less than distance from notochord Score =1 
                    wound1234_interval(numTrack)              = 1;
                elseif dist_notochord < startwound_um_interval(numTrack)&& startwound_um_interval(numTrack) <=(dist_notochord+ParameterS)
                    wound1234_interval(numTrack)              = 2;
                elseif (dist_notochord+ParameterS) < startwound_um_interval(numTrack)&& startwound_um_interval(numTrack) <=(dist_notochord+ParameterS+ParameterS) 
                     wound1234_interval(numTrack)              = 3;
                elseif (dist_notochord+ParameterS+ParameterS) < startwound_um_interval(numTrack)&& startwound_um_interval(numTrack) <=(dist_notochord+ParameterS+ParameterS+ParameterS+1000) 
                     wound1234_interval(numTrack)              = 4;
                end;
                %--for cummulative
                withzeroeswounddistnodescum=0;
                for i=1:upperlimit_idx
                    x_node1                             = HNN{1}(track_idxs(i),1);
                    y_node1                             = HNN{1}(track_idxs(i),2);
                    withzeroeswounddistnodescum(i)      = pdist([x_node1,y_node1; x_nodeW,y_nodeW],'Euclidean');   
                end; 
                wound_dist_nodes_cum{numTrack}=withzeroeswounddistnodescum(withzeroeswounddistnodescum~=0);
                %(0)=static|(+)=backward|(-)=forward
                distdiff_cum{numTrack}          = diff(wound_dist_nodes_cum{numTrack});
                lens_cum                        = length(distdiff_cum{numTrack} ); %Total
                positive_backward_cum           = distdiff_cum{numTrack}(distdiff_cum{numTrack}>0); %greater than 0
                negative_forward_cum            = distdiff_cum{numTrack}(distdiff_cum{numTrack}<0); %less than than 0
                forwardratio_cum(numTrack)      = length(negative_forward_cum)/lens_cum;  
                backwardratio_cum(numTrack)     = length(positive_backward_cum)/lens_cum;
                ftobratio_cum(numTrack)         = length(negative_forward_cum)/length(positive_backward_cum);
                wound_um_cum(numTrack)          =(diff([wound_dist_nodes_cum{numTrack}(1) wound_dist_nodes_cum{numTrack}(end)]))*coe; %This is in um already - subtracts first node distance from wound from last node statnce from wound
                startwound_um_cum(numTrack)     = wound_dist_nodes_cum{numTrack}(1); % we take the start point only and then create score
                %--To Determine The Wound Score 1 2 3 4 by 150um
                if startwound_um_cum(numTrack) == 0 || isnan(startwound_um_cum(numTrack));
                    wound1234_cum(numTrack)              = nan;
                elseif startwound_um_cum(numTrack)<=dist_notochord %if less than distance from notochord Score =1 
                    wound1234_cum(numTrack)              = 1;
                elseif dist_notochord < startwound_um_cum(numTrack)&& startwound_um_cum(numTrack) <=(dist_notochord+ParameterS)
                    wound1234_cum(numTrack)              = 2;
                elseif (dist_notochord+ParameterS) < startwound_um_cum(numTrack) && startwound_um_cum(numTrack) <=(dist_notochord+ParameterS+ParameterS) 
                     wound1234_cum(numTrack)              = 3;
                elseif (dist_notochord+ParameterS+ParameterS) < startwound_um_cum(numTrack) && startwound_um_cum(numTrack)<=(dist_notochord+ParameterS+ParameterS+ParameterS+1000) %Should I make this till the end??
                     wound1234_cum(numTrack)              = 4;
                else
                    disp(numTrack)
                    disp('Andre you suck')
                end;               
            end; 
%              disp(strcat('Numtrack:',num2str(numTrack),' startwound_um_cum:',num2str(startwound_um_cum(numTrack)), 'wound1234_cum:',num2str(wound1234_cum(numTrack))))           
        end;
        WoundStartUm.interval{numGT}    = startwound_um_interval; %This is proof 08/14/17
        WoundStartUm.cum{numGT}         = startwound_um_cum; % This is proof 08/14/17
        WoundScore1234.cum{numGT}       = wound1234_cum;
        WoundScore1234.interval{numGT}  = wound1234_interval;
        WoundScoreUm.cum{numGT}         = wound_um_cum; %in um  % This is the 6th Metric (difference from start to finish
        WoundScoreUm.interval{numGT}    = wound_um_interval; % in um
        WoundDist.cum{numGT}            = wound_dist_nodes_cum; % This is a metric checker - lists all the distances from node to wound
        WoundDist.interval{numGT}       = wound_dist_nodes_interval;
        NodeDistDiff_cum{numGT}         = distdiff_cum;
        NodeDistDiff_interval{numGT}    = distdiff_interval;
        Forward.cum{numGT}              = forwardratio_cum;
        Forward.interval{numGT}         = forwardratio_interval;
        Backward.cum{numGT}             = backwardratio_cum;
        Backward.interval{numGT}        = backwardratio_interval;
        FtoB.cum{numGT}                 = ftobratio_cum;
        FtoB.interval{numGT}            = ftobratio_interval;
        FBratio.cum{numGT}              = forwardratio_cum - backwardratio_cum;
        FBratio.interval{numGT}         = forwardratio_interval - backwardratio_interval;
     else 
        WoundStartUm.interval{numGT}    = nan;
        WoundStartUm.cum{numGT}         = nan;
        WoundScore1234.cum{numGT}       = nan;
        WoundScore1234.interval{numGT}  = nan;        
        WoundScoreUm.cum{numGT}         = nan;
        WoundScoreUm.interval{numGT}    = nan;
        WoundDist.cum{numGT}            = nan;
        WoundDist.interval{numGT}       = nan;
        NodeDistDiff_cum{numGT}         = nan;
        NodeDistDiff_interval{numGT}    = nan;
        Forward.cum{numGT}              = nan;
        Forward.interval{numGT}         = nan;
        Backward.cum{numGT}             = nan;
        Backward.interval{numGT}        = nan;
        FtoB.cum{numGT}                 = nan;
        FtoB.interval{numGT}            = nan;
        FBratio.cum{numGT}              = nan;
        FBratio.interval{numGT}         = nan;
     end;
     SM.WoundStartUm.interval{numGT}    = WoundStartUm.interval{numGT}(gt_SM{numGT}); %These were just to proof, looks good now 08/14/2017
     SM.WoundStartUm.intervalsm{numGT}  = WoundStartUm.interval{numGT}(Selected);
     SM.WoundStartUm.cum{numGT}         = WoundStartUm.cum{numGT}(Selected);
     SM.WoundScore1234.interval{numGT}  = WoundScore1234.interval{numGT}(gt_SM{numGT});
     SM.WoundScore1234.intervalsm{numGT}= WoundScore1234.interval{numGT}(Selected);
     SM.WoundScore1234.cum{numGT}       = WoundScore1234.cum{numGT}(Selected);
     SM.WoundScoreUm.interval{numGT}    = WoundScoreUm.interval{numGT}(gt_SM{numGT});
     SM.WoundScoreUm.intervalsm{numGT}  = WoundScoreUm.interval{numGT}(Selected);
     SM.WoundScoreUm.cum{numGT}         = WoundScoreUm.cum{numGT}(Selected);
     SM.forward.interval{numGT}         = Forward.interval{numGT}(gt_SM{numGT});
     SM.forward.intervalsm{numGT}       = Forward.interval{numGT}(Selected);
     SM.forward.cum{numGT}              = Forward.cum{numGT}(Selected);
     SM.backward.interval{numGT}        = Backward.interval{numGT}(gt_SM{numGT});
     SM.backward.intervalsm{numGT}      = Backward.interval{numGT}(Selected);
     SM.backward.cum{numGT}             = Backward.cum{numGT}(Selected);
     SM.FtoB.interval{numGT}            = FtoB.interval{numGT}(gt_SM{numGT});
     SM.FtoB.intervalsm{numGT}          = FtoB.interval{numGT}(Selected);
     SM.FtoB.cum{numGT}                 = FtoB.cum{numGT}(Selected);
     SM.FBratio.interval{numGT}         = FBratio.interval{numGT}(gt_SM{numGT});
     SM.FBratio.intervalsm{numGT}       = FBratio.interval{numGT}(Selected);
     SM.FBratio.cum{numGT}              = FBratio.cum{numGT}(Selected);
     
     %Zirmi includes all and respective positions in arrays, but SM cleans
     %it for simple incorporation for Zirmi C2 Script. 
     Zirmi.WoundStartUm.interval{numGT}    = WoundStartUm.interval{numGT}; %These were just to proof, looks good now 08/14/2017
     Zirmi.WoundStartUm.intervalsm{numGT}  = WoundStartUm.interval{numGT};
     Zirmi.WoundStartUm.cum{numGT}         = WoundStartUm.cum{numGT};
     Zirmi.WoundScore1234.interval{numGT}  = WoundScore1234.interval{numGT};
     Zirmi.WoundScore1234.intervalsm{numGT}= WoundScore1234.interval{numGT};
     Zirmi.WoundScore1234.cum{numGT}       = WoundScore1234.cum{numGT};
     Zirmi.WoundScoreUm.interval{numGT}    = WoundScoreUm.interval{numGT};
     Zirmi.WoundScoreUm.intervalsm{numGT}  = WoundScoreUm.interval{numGT};
     Zirmi.WoundScoreUm.cum{numGT}         = WoundScoreUm.cum{numGT};
     Zirmi.forward.interval{numGT}         = Forward.interval{numGT};
     Zirmi.forward.intervalsm{numGT}       = Forward.interval{numGT};
     Zirmi.forward.cum{numGT}              = Forward.cum{numGT};
     Zirmi.backward.interval{numGT}        = Backward.interval{numGT};
     Zirmi.backward.intervalsm{numGT}      = Backward.interval{numGT};
     Zirmi.backward.cum{numGT}             = Backward.cum{numGT};
     Zirmi.FtoB.interval{numGT}            = FtoB.interval{numGT};
     Zirmi.FtoB.intervalsm{numGT}          = FtoB.interval{numGT};
     Zirmi.FtoB.cum{numGT}                 = FtoB.cum{numGT};
     Zirmi.FBratio.interval{numGT}         = FBratio.interval{numGT};
     Zirmi.FBratio.intervalsm{numGT}       = FBratio.interval{numGT};
     Zirmi.FBratio.cum{numGT}              = FBratio.cum{numGT};
 end
SM.forward.phagosight    = handles.distanceNetwork.forwardRatioTot(Selected);
SM.backward.phagosight   = handles.distanceNetwork.backwardRatioTot(Selected);
  
disp('END: Foward and Backward Determinations');
%% Frame (Velocity & Static Ratio)
clear Frame
perframe=(0);
perframe_nonzero=(0);
for i=1:frames
    %-velocity
    perframe=Thismap(i,:);
    perframe_nonzero=perframe(perframe~=0);
    %-static
    s_perframe=Staticmap(i,:);
    s_perframe_nonnan=s_perframe(~isnan(s_perframe));    
    if isempty(perframe_nonzero)
        perframe_nonzero=nan;
        s_perframe_nonnan=nan;
    else
    end;
    Velocity.frame.perframe_nonzero{i}=perframe_nonzero;
    StaticRatio.frame.s_perframe_nonnan{i}=s_perframe_nonnan;
end;
Frame.velocity=Velocity.frame;
Frame.staticratio=StaticRatio.frame;
disp('End: Frame')
%% For Plotting Worthy 
%need to make arrays similar sizes 
%will fill in with
clear Worthy
MaxTracks   = mp;

lens        = length(str_ArrGT);

hops        = frames*mp;
for i=1:lens
    %--Interval Lengths (They should all be the same)
    len_1                   = mp-length(SM.velocity.interval{i});       % Metric #1 Velocity
    len_2                   = mp-length(SM.staticratio.interval{i});    % Metric #2 Static Ratio
    len_3                   = mp-length(SM.meandering.interval{i});     % Metric #3 Meandering Index
    len_4                   = mp-length(SM.tortuosity.interval{i});     % Metric #4 Tortuosity
    len_5                   = mp-length(SM.forward.interval{i});        % Metric #5 Forward
    len_6                   = mp-length(SM.FBratio.interval{i});        % Metric #6 FBratio
    len_7                   = mp-length(SM.WoundScoreUm.interval{i});   % Metric #7 Start-End Toward Wound
    len_8                   = mp-length(SM.WoundScore1234.interval{i}); % Metric #8 Assigned Scores           
    %--------------------------------------------------    
    lenc_1                  = mp-length(SM.velocity.cum{i});
    lenc_2                  = mp-length(SM.staticratio.cum{i});   
    lenc_3                  = mp-length(SM.meandering.cum{i});
    lenc_4                  = mp-length(SM.tortuosity.cum{i});
    %---------------------------------------------------
    Worthy.SM.velocity.interval{i}      =    horzcat(SM.velocity.interval{i},nan(1,len_1));     %Metric #1 velocity
    Worthy.SM.velocity.cum{i}           =    horzcat(SM.velocity.cum{i},nan(1,lenc_1));         %Metric #1 
    Worthy.SM.staticratio.interval{i}   =    horzcat(SM.staticratio.interval{i},nan(1,len_2));  %Metric #2 static ratio
    Worthy.SM.staticratio.cum{i}        =    horzcat(SM.staticratio.cum{i},nan(1,lenc_2));      %Metric #2  
    Worthy.SM.meandering.interval{i}    =    horzcat(SM.meandering.interval{i},nan(1,len_3));   %Metric #3 Meandering
    Worthy.SM.meandering.cum{i}         =    horzcat(SM.meandering.cum{i},nan(1,lenc_3));       %Metric #3  
    Worthy.SM.tortuosity.interval{i}    =    horzcat(SM.tortuosity.interval{i},nan(1,len_4));   %Metric #4 Tortuosity
    Worthy.SM.tortuosity.cum{i}         =    horzcat(SM.tortuosity.cum{i},nan(1,lenc_4));       %Metric #4 
    Worthy.SM.forward.interval{i}       =    horzcat(SM.forward.interval{i},nan(1,len_5));      %Metric #4 Forward
    Worthy.SM.forward.cum{i}            =    horzcat(SM.forward.cum{i},nan(1,lenc_1));          %Metric #4 
    Worthy.SM.FBratio.interval{i}       =    horzcat(SM.FBratio.interval{i},nan(1,len_6));      %Metric #4 FBratio
    Worthy.SM.FBratio.cum{i}            =    horzcat(SM.FBratio.cum{i},nan(1,lenc_1));          %Metric #4 
    Worthy.SM.WoundScoreUm.interval{i}  =    horzcat(SM.WoundScoreUm.interval{i},nan(1,len_7)); %Metric #4 WoundScoreUm
    Worthy.SM.WoundScoreUm.cum{i}       =    horzcat(SM.WoundScoreUm.cum{i},nan(1,lenc_1));     %Metric #4 
    Worthy.SM.WoundScore1234.interval{i}=    horzcat(SM.WoundScore1234.interval{i},nan(1,len_7)); %Metric #4 WoundScoreUm
    Worthy.SM.WoundScore1234.cum{i}     =    horzcat(SM.WoundScore1234.cum{i},nan(1,lenc_1));     %Metric #4 
   %--------------------------------------------------
    len_Fveli    = hops-length(Frame.velocity.interval{i});
    len_Fstai    = hops-length(Frame.staticratio.interval{i});    
    len_Fvelc    = hops-length(Frame.velocity.cum{i});
    len_Fstac    = hops-length(Frame.staticratio.cum{i}); 
    %------
    Worthy.Frame.velocity.interval{i}      =    horzcat(Frame.velocity.interval{i}',nan(1,len_Fveli));
    Worthy.Frame.velocity.cum{i}           =    horzcat(Frame.velocity.cum{i}',nan(1,len_Fvelc));
    Worthy.Frame.staticratio.interval{i}   =    horzcat(Frame.staticratio.interval{i}',nan(1,len_Fstai));
    Worthy.Frame.staticratio.cum{i}        =    horzcat(Frame.staticratio.cum{i}',nan(1,len_Fstac));
end;

display('End: Worthy')
%% Histogram of Hops
% %All this is is a histogram of all HOPS
% a=Frame.velocity.cum{3};
% [sortedValues,~] = sort(a,1,'descend');
% exclude = sortedValues > 12;
% sortedValues(exclude) = [];
% highestValue = sortedValues*coe;
% 
% 
%     FigHandle = figure;
%     set(FigHandle, 'Position', [100, 100, 350, 250]); %[x1, y1, x2, y2]
%      histogram(highestValue,50);
%      ylabel('Absolute Velocity(um/min)')
%     title(strcat('Histogram of velocity(um/min)  30-',str_ArrGT{3}))
%%  Clear vars
Worthy.time      = time;
ID.interval      = gt_SM;
ID.intervalsm    = Selected;
ID.cum           = Selected;
ID.phagosight    = (1:mp);
Worthy.ID        = ID;
%Note:  Deleting Structure Arrays: 1.Andre  2.Phagosight 3.distmap 4.Jun
 clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 time Zirmi equationCTF WoundScoreUm
display ('FINISHED: Zirmi_C - Compute Metrics')    