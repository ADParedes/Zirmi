function [arr_T,arr_boo_GT,arr_GT,arr_LoopGT,frames ] = function_calculateGT( d_positionName,SamplingFrequency,MPI_start)
%% Summary of function_calculateGT
% This function is Used in Zirmi_D Scripts for CTF calculations
% The following time points are determined in the arrays
    %index 1) GT_60 
    %index 2) GT_90
    %index 3) GT_120
    %index 4) GT_150 
    %index 5) GT_180
    %index 6) GT_max
%Outputs
%       0: arr_T       - this is important; it lists every Time Folder - so if you delete a time folder there won't be an error :) 
%       1: arr_boo_GT  - this is a boo array of a present frame from each
%       index i.e. 1-6 listed above
%       2: arr_GT      - This is the Array of Frame# for this unique fish 
%       3: arr_LoopGT  - Includes frame numbers that are only REAL or
%       present for this fish %--this is most useful.
%       3: frames      - Total number of frames in this specific fish
%%
len_dpn                             = length(d_positionName);
% [temp ]                             = calculateGT(d_positionName);
for i=(1:len_dpn) %for MAC is at 5.  
    dir_time                        = d_positionName(i).name;
    arr_T(i)                        = str2num(dir_time(2:end));
    min_t                           = min(arr_T);
    max_t                           = max(arr_T);
end;
adj_tplate          = (min_t-1)*SamplingFrequency+MPI_start;%min to add to t_plate so that deleted frames are considered in time 
tplate              = adj_tplate;
%--max is just arr_T length
frames              = length(arr_T); 
disp(frames)
arr_time                = round((arr_T(1:frames)*SamplingFrequency+MPI_start)-SamplingFrequency);
% frames=99;
%--determine the time points after wounding.  
GT_60               = round((60-tplate)/SamplingFrequency);
GT_90               = round((90-tplate)/SamplingFrequency);
GT_120              = round((120-tplate)/SamplingFrequency);
GT_150              = round((150-tplate)/SamplingFrequency);  %GT_150=golden(2hours and 30 minutes)
GT_180              = round((180-tplate)/SamplingFrequency);
boo_60              = 1;
boo_90              = 1;
boo_120             = 1;
boo_150             = 1;
boo_180             = 1;
boo_max             = 1;
display('_')
if GT_60>frames
    display(strcat(num2str(GT_120),'_ GT_60 exceeds ',' _', num2str(frames)))
    boo_60  = 0;
    boo_90  = 0;
    boo_120 = 0;
    boo_150 = 0;
    boo_180 = 0;
elseif GT_90>frames
    display(strcat(num2str(GT_120),'_ GT_90 exceeds ',' _', num2str(frames)))
    boo_90  = 0;
    boo_120 = 0;
    boo_150 = 0;
    boo_180 = 0;
elseif GT_120>frames
    display(strcat(num2str(GT_120),'_ GT_120 exceeds ',' _', num2str(frames)))
    boo_120 = 0;   
    boo_150 = 0;
    boo_180 = 0;
elseif GT_150>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_150 = 0;
    boo_180 = 0;
elseif GT_150>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_150 = 0;
    boo_180 = 0;
elseif GT_180>frames
    display(strcat(num2str(GT_150),'_GT_150 exceeds ',' _', num2str(frames)))
    boo_180 = 0;
else
    display('you have over 3 hours of total experiment.  Wow')
end;  
disp('_')
arr_boo_GT      = [boo_60,boo_90,boo_120,boo_150,boo_180,boo_max];
arr_GT          = [GT_60,GT_90,GT_120,GT_150,GT_180,frames];
gtinwoundcount  = 0;
beat            = 0;
%-Update to only do these time points
boo_max         = 0;
arr_LoopGT      = arr_GT(1==(arr_boo_GT)); %includes only the values that arn't zero


end

