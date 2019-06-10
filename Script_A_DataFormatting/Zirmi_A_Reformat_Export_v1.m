%% Zirmi Package - Module (1) Reformatting Raw Mirosope images
% Version 1.3 11-20-17
% Files needed be exported in diretories  by channel (BF, GFP, TexRed)
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%%-Zirmi Software Workflow:
    %%-Module (1) Is Initial Formatting of Pre-processing Data
    %- This is intended to take data Organized in manner outlined and place in
    %  a directory that the Zirmi program can use to throughput computations
    %%-Module (2)Is to Quantitate Corrective Total Fluorescence Quantitation
    %%-Module (3)Is Computing Cell Migration based Metrics 
    %  Is PhagoSight Based Scripts and initial Processing of Fish 
    %%-Module (4) Export to Exel (ammendable file) 
%% Close all/Clear Variables/Warning
clear ; 
clc; 
close all;
boo_pc                     =   ispc;

switch      boo_pc
    case {1}
        disp('Zirmi is Compatible with your PC')
        disp('NOTE: Zirmi is not able to process more than 99 positions')
        pause(2)
    otherwise
        disp('Zirmi is not yet compatible with your MAC...YET.')
        disp('Look for Future updates')
        display                   ('Module 1 DISCONTINUED')
        return
        
end;  
%% Computer System Identification 
!hostname
    hostname                    =   char( getHostName( java.net.InetAddress.getLocalHost ) );
% hostname = 'A';
ip                          =   char( getHostAddress( java.net.InetAddress.getLocalHost ) );
user                        =   getenv('UserName');    
dir_desktop                 = strcat('C:\Users\',user,'\Desktop');
cd                          (dir_desktop)
%% Warning
disp                    ('Please Note the Following Warnings : ')
warning                 ('This is Module - Zirmi A (Module 1 Part 1) and needs to be ran reformatt Files In a callable form')
check_functions         = exist('zirmiChVal'); %#ok<EXIST>
switch check_functions
    case {0} 
        disp            ('You havent loaded Zirmi functions to matlab, the following is a temporary solution')
        warning         ('You need set path to Zirmi functions!!')
        
        disp            ('Select Directory of Zirmi Directory: "Zirmi"')
        dir_Zirmi       = uigetdir();
        dir_zirmiFuns   = strcat(dir_Zirmi,'\','Script_E_Functions');
        addpath          (dir_zirmiFuns);%add path to *.mfile folder and *.tif foldersss
    otherwise
end;
%% Save Time part
input_clearall                  =   input('Should I clear all?  y/n || need "y" if first go ','s');
if input_clearall               ==  'y'   
    %% Define Directory Part I
    cd                          (dir_desktop)
    str_dirZirmiFormat          = 'Zirmi_Data\Data_Formatted';
    mkdir                        (str_dirZirmiFormat)
    dir_default                 = strcat(dir_desktop,'\',str_dirZirmiFormat); 
    if strcmp(hostname,'ADP-Chicago'); %this is the desktop
        disp                        ('Did NOT select Labtop Directory')
        boo_desktop                 =   1; %manual save logic of computer type | "1" is Yes for this is desktop
    elseif strcmp(hostname,'ADPAREDES') %this is the labtop
        disp                        ('Hello ANDRE: Selected Labtop Directory')

        dir_ZirmiFormattedData      =   strcat('C:\Users\',user,'\OneDrive\MrBigTest\Confocal Data\Processed Data'); 
        boo_desktop                 =   0;
    else
        clc
        disp                        (strcat('Hello New User: ',hostname,': We Need to locate TWO Directories'))
        warning                     ('We need to locate Your (1) Where you want to save reformatted files  (2) locate your specific batch of Raw GreyScale TIFFs')
        disp                        ('(1) Where would you like to save: Current Default Saves to Desktop')
        pause                       (2)
        dir_RawBatchData            = strcat('C:\Users\',user,'\Desktop'); %Experiment Folder- Starting Point to find Experiment to Analyze
        dir_ZirmiFormattedData      = uigetdir();
        boo_desktop                 = 2;
    end;
    %% Defining Directories Part II
    disp                        ('(2) Where are your Microscope exported RAW greyscale exported TIFFs (needs to umbrella Channel Directories)')
    dir_exportedBatchBranch         =   uigetdir();   
    [pa, name, ext]                 =   fileparts(dir_exportedBatchBranch); %break up file name
    figurename                      =   name;
    cd                              (dir_exportedBatchBranch);%change the current folder    
    %% Determine Channel Names 
    clear channels
    prompt                      = {'Imaging Channel #1 : BrightField or DIC or Ph2 based channel:'...
                                ,'Imaging Channel #2 : Fluorescence of Probe'...
                                ,'Imaging Channel #3 : Fluorescence of cell'};                         
    dlg_title                   = 'Please input the name of Imaging Channel #1-3 as labeled on the directory e.g.(BF,GFP,TexRed)';
    num_lines                   = 1;
    defaultans                  = {'BF','TexRed','GFP'};
    answer                      = inputdlg(prompt,dlg_title,num_lines,defaultans);
    channels.BF                 = answer{1}; % Brightfield
    channels.Probe              = answer{2}; % TexRed
    channels.Cells              = answer{3}; % GFP        
    %% Defining Variables
    posnum          = 0; %#ok<NASGU>
    parts           = 0; %#ok<NASGU>
    znum            = 0; %#ok<NASGU>
    zPh2num         = 0; %#ok<NASGU>
    zTRnum          = 0; %#ok<NASGU>
    zGFPnum         = 0; %#ok<NASGU>
  %- Just in case - %
    Zname=1; %To ensure it still takes into account the name spacing
    PH2check   =  0;
    GFPcheck   =  0;
    TRcheck    =  0;%isempty(dTexRed); %If empty will be 1;   
    %% Determin BF Dir & #Positions,#Parts,#Zstacks, and #Frames
    dir_BF                      = strcat(dir_exportedBatchBranch,'\',channels.BF);
    cd                         (dir_BF);
    dBF                         = dir('*.tif'); % '*.tif' select only tif files in directory
    posnum                      = zirmiChVal('Part*_p'); %Position --user needs to be aware of this
    parts                       = zirmiChVal('*Part');
    zPh2num                     = zirmiChVal('*_z');
    cd                              (dir_exportedBatchBranch)       
        if posnum==2
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('  NOTE: Possibly More than 10 Positions  ')
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            pos_question=input('Are there really only 2 Positions?','s');
            if pos_question=='y';
            else
                    check                       =   0;
                    posnum                      =   -1;
                    while check == 0
                            posnum      =    posnum+1;                           
                           if posnum<10
                               ev          =    strcat('Part*_p0',num2str(posnum),'*'); %Problem: Preinjury and Parts have different Posnum
                           elseif posnum<100 
                               ev          =    strcat('Part*_p',num2str(posnum),'*'); %Problem: Preinjury and Parts have different Posnum
                           end;
                            f           =    dir(ev);
                            check       =    isempty(f);
                    end;   
            end;
        else
        end;
     %posnum   =  posnum-1;%-1; Don't minus 1 because p0 is a position.
    %% Probe Determining Color Filter Number of Parts, Positions, and Zstacks
    %-Probe 
    dir_TexRed                  = strcat(dir_exportedBatchBranch,'\',channels.Probe);
    cd                          (dir_TexRed)
    dTexRed                     = dir('*.tif'); % '*.tif' select only tif files in directory
    znum                        = zirmiChVal('*_z');
    zTRnum                      =  znum; %Incorporate in the Future?
    %% Cells Determining Color Filter Number of Parts, Positions, and Zstacks
    dir_GFP                     = strcat(dir_exportedBatchBranch,'\',channels.Cells);
    cd                           (dir_GFP);
    dGFP                        = dir('*.tif'); % '*.tif' select only tif files in directory
    zGFPnum                     = zirmiChVal('*_z');
    %% Separate into Arrays
    cds             ={dir_BF,dir_GFP,dir_TexRed};
    prev_folder     = fileparts(dir_exportedBatchBranch);
    cd              (prev_folder);               
    cd              (dir_exportedBatchBranch);
    Arr_Ph2         ={};
    Arr_GFP         ={};  %{Part}{Position}{All in Time001}
    Arr_TexRed      ={};
    for I=1:parts
        %% First For Loop Section
        cd(dir_BF)
        numfid=0;
        if I==parts
            str_shmidt='*njury';
            str_part=strcat(str_shmidt);
        else
            str_shmidt='Part';
            str_part=strcat(str_shmidt,num2str(I));
        end;
        allfids=length(dGFP);
        allfids(I)=allfids;
        Green={}; TeRe={};
        Green2={};
        GreenZ={};
        TeRe1={};
        TeReZ={};
        disp(str_part)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        %% Secondary For Loops
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:posnum 
            i=j-1;            
            switch j
                case{1}
                    if i<10
                        i = '0';
                    else
                        i = '00';
                        disp('Zirmi not capable of more than 99 fish/batch')
                    end;
                otherwise                            
                    if i<10
                        i=num2str(i);   
                    else
                        i=strcat('0',num2str(i));
                    end;
            end;
            num_middleImg   = round(znum/2);
            namePh2         = strcat(str_part,'*_z',num2str(num_middleImg),'*','p',i,'*'); 
            nameTR          = strcat(str_part,'*p',i,'*');
            nameGFP         = strcat(str_part,'*p',i,'*');
            arr_strPositionName{j}=strcat('P',i);
            %% Ph2 
            cd(dir_BF)
            ff          = dir(namePh2);
            numfids(I)  = length(ff);                
             DD         = [ff(:).datenum].'; % you may want to eliminate . and .. first.
            [DD,DD]     = sort(DD);
            PH2{j}      = {ff(DD).name}; % Cell array of names in order by datenum.
            %Note PH2 is arrayed {str} by Positions and then (str) frames.  No Z's
            %chose middle z~4 since we don't care about the stack at
            %this point.
            cd(dir_exportedBatchBranch);
            %% TexRed
            cd            (dir_TexRed);
            if I~=parts 
                gg      = dir(nameTR);
                EE      = [gg(:).datenum].'; % you may want to eliminate . and .. first.
                [EE,EE] = sort(EE);
                TeRe{j} = {gg(EE).name}; % Cell array of names in order by datenum.  TeRe{Position}{TimeInterval
                if znum>1
                    %% Z - TexRed
                    [TeRe1] = zirmiAdelim( numfids,I,i,str_part);
                    TeReZ{j}    = TeRe1;   % %TeReZ{Position}{Z}{frame}
                else    
                %% No TexRed Z
                display('No TexRed Z')
                end;
            else
                %% PRE INJURY%%%%%%%%%%%%%%%%%%
                %Preinjury has to be less than 10 frames or this won't work
                gg      = dir(nameTR);
                EE      = [gg(:).datenum].'; % you may want to eliminate . and .. first.
                [EE,EE] = sort(EE);
                TeRe{j} = {gg(EE).name}; % Cell array of names in order by datenum.  TeRe{Position}{TimeInterval
                if znum>1
                    for bi = 1:numfids(I)
                        ti=bi-1;
                       %--------Z---------
                       if bi==1
                           nameTR1  = strcat(str_part,'*t0','*p',num2str(i),'*');
                       elseif ti<10
                            nameTR1  = strcat(str_part,'*t',num2str(ti),'*p',num2str(i),'*');
                       else
%                             nameTR1  = strcat(str_part,'*t0',num2str(ti),'*P',num2str(i),'*');
                       end;
                        gg1       = dir(nameTR1);
                        EE1       = [gg1(:).datenum].'; % you may want to eliminate . and .. first.
                        [EE1,EE1] = sort(EE1);
                        TeRe1{bi} = {gg1(EE1).name}; % Cell array of names in order by datenum.  TeRe{Position}{TimeInterval               
                     end;
                    TeReZ{j}    = TeRe1;   % %TeReZ{Position}{Z}{frame}
                else
                end;
             end;
             cd(dir_exportedBatchBranch);
            %% GFP 
            cd(dir_GFP);
            if I~=parts
                hh=dir(nameGFP);
                HH = [hh(:).datenum].'; % you may want to eliminate . and .. first.
                [HH,HH] = sort(HH);
                Green{j} = {hh(HH).name}; % Cell array of names in order by datenum.  Green{Position}{TimeInterval
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if zGFPnum>1
                    %% Z - GFP
                    [Green2] = zirmiAdelim(numfids,I,i,str_part);
                    GreenZ{j}=Green2; %{Position}{all Z in time 1}       
                else
                    %% No GFP Z
                    display('No GFP Z')
                end;            
            else
                %% PRE INJURY%%%%%%%%%%%%%%%%%%
                hh=dir(nameGFP);
                HH = [hh(:).datenum].'; % you may want to eliminate . and .. first.
                [HH,HH] = sort(HH);
                Green{j} = {hh(HH).name}; % Cell array of names in order by datenum.  Green{Position}{TimeInterval
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if zGFPnum>1
                    for bi=1:numfids(I)
                       ti=bi-1;
                       %--------Z---------
                       if bi==1
                           nameGFP1  = strcat(str_part,'*t0','*p',num2str(i),'*');
                       elseif bi<10
                            nameGFP1  = strcat(str_part,'*t',num2str(ti),'*p',num2str(i),'*');
                       else
%                             nameTR1  = strcat(str_part,'*t0',num2str(ti),'*P',num2str(i),'*');
                       end;
                        mn=dir(nameGFP1);
                        MN = [mn(:).datenum].'; % you may want to eliminate . and .. first.
                        [MN,MN] = sort(MN);
                        Green2{bi} = {mn(MN).name}; % Cell array of names in order by datenum.  TeRe{Position}{TimeInterval                
                    end;
                    GreenZ{j}=Green2; %{Position}{all Z in time 1}       
                else
                end;
            end;
            cd(dir_exportedBatchBranch);
        end;
        Arr_Ph2{I}=PH2;
        Arr_GFP{I}=GreenZ;  %{Part}{Position}{All Z in Time001} therefore to call Arr_GFP{parts}{posnum}{numfid}(znum) **note parenthesis.
        Arr_TexRed{I}=TeReZ;
    end;
    PositionName=arr_strPositionName;
    disp('End: Sorting Names into Arrays')
end;
pause(2)
%% imreading images  
%numfidsPH = numel(dPh2); %number of tif files in folder
valsPhaseContrast       =   {};
valsPH2                 =   {};
valsTR                  =   {};
valsTexasRed            =   {};
valsGFP                 =   {};
valsGreen               =   {};
%%--Reads the grayscale/color image from the file specified by
%%the string FILENAME.  REteruns value in an array containing the image data
%% Create Folders and Save BF 
cd                    (dir_ZirmiFormattedData);
mkdir                 (name);%---> Move to processed data folder
cd                    (name)%-----> Move to experiment folder
ProcessedTemp_str   = cd; 
valsPh2             = {};
for J = 1:posnum 
    cd                   (ProcessedTemp_str) %Directory of Processed Experiment#
    Str_Posnum         = PositionName{J};
    mkdir                (Str_Posnum)
    cd                   (Str_Posnum)
    pos_folder         = cd; %Processed data Position directory
    [pa3, name3, ext3] = fileparts(cd);
    counter            = 1000;
    cnter              = 0;
    disp                (strcat('On Position',num2str(J),'!!!!!!!!!!!!!!'))
    valsPHZ            = {};
    for K=1:parts
        if K==parts % --> This is so that we can put the PreInjury In its own folder
                Str_Posnum         = PositionName{J};
                Str_PreInjuryPosnum= strcat(Str_Posnum,'-PreInjury');
                cd                   (ProcessedTemp_str)            
                mkdir                (Str_PreInjuryPosnum)
                cd                   (Str_PreInjuryPosnum)
                pos_folder         = cd; %Processed data Position directory
        end;
        disp                   (strcat('On Part',num2str(K),'/ out of',num2str(parts))) %displaying progress
        minfids            =   numfids(K);
        for j = 1:minfids
            %Count According to Part numfids cummulatively
            counter            =  counter+1; %IMPORTANT--COUNTING MECHANISM FOR FOLDERS accumulation
            cnter              =  cnter+1;
            cd                    (pos_folder) %Directory of Processed 'P#'
            disp                  (counter)
            %Change Folder Name to T000# format
            str_Tfoldername    =  strcat('T',num2str(counter));
            str_Tfoldername(2) =  '0';
            %Create Folder with Such Name and go to Directory
            mkdir                 (str_Tfoldername)                   
            cd                    (str_Tfoldername)               
            [pa4, name4, ext4] =  fileparts(cd);
            time_folder        =  cd;           
            valsBF             =  {0};
            %Grab Name and add T000# Format to such names
            str_Ph2name        =  Arr_Ph2{K}{J}{j};
            cd                   (dir_BF)
            valsPHZ{cnter}         =  imread(str_Ph2name); 
            cd                   (time_folder);
            str_SaveName       = strcat(str_Tfoldername,'-Ph2-',str_Ph2name);
            imwrite              (valsPHZ{cnter},str_SaveName,'Compression','none','Resolution',1);
            disp                 ([num2str(100*j/(minfids),'%2.1f') ' %BF' num2str(j) ' / ' num2str(minfids) ' ' num2str(J) ' out ' num2str(posnum) ' ' num2str(100*J/(posnum),'%2.1f')]);
        end;
    end;
    valsPh2{J}=valsPHZ;
%    valsTR{J}=valsTRnumfid;
%    valsGFP{J}=valsGnumfid;
end;
%% Save Zstacks GFP
cd                    (dir_ZirmiFormattedData);
mkdir                 (name);
cd                    (name)
ProcessedTemp_str   = cd; 
for J = 1:posnum 
    cd                   (ProcessedTemp_str) %Directory of Processed Experiment#
    Str_Posnum         = PositionName{J};
    mkdir                (Str_Posnum)
    cd                   (Str_Posnum)
    pos_folder         = cd; %Processed data Position directory
    [pa3, name3, ext3] = fileparts(cd);
    counter            = 1000;
    cnter              = 0;
    disp                (strcat('On Position',num2str(J),'!!!!!!!!!!!!!!'))
    for K=1:parts
        if K==parts % --> This is so that we can put the PreInjury In its own folder
                Str_Posnum         = PositionName{J};
                Str_PreInjuryPosnum= strcat(Str_Posnum,'-PreInjury');
                cd                   (ProcessedTemp_str)            
                mkdir                (Str_PreInjuryPosnum)
                cd                   (Str_PreInjuryPosnum)
                pos_folder         = cd; %Processed data Position directory
        end;
        disp                        (strcat('On Part',num2str(K),'/ out of',num2str(parts))) %displaying progress
        minfids                 =   numfids(K);
        for j = 1:minfids
            %Count According to Part numfids cummulatively
            counter             =  counter+1; %IMPORTANT--COUNTING MECHANISM FOR FOLDERS accumulation
            cnter               =  cnter+1;
            cd                    (pos_folder) %Directory of Processed 'P#'
            disp                  (counter)
            %Change Folder Name to T000# format
            str_Tfoldername     =  strcat('T',num2str(counter));
            str_Tfoldername(2)  =  '0';
            %Create Folder with Such Name and go to Directory
            mkdir                 (str_Tfoldername)                   
            cd                    (str_Tfoldername)               
            [pa4, name4, ext4]  =  fileparts(cd);
            time_folder         =  cd;           
            valsGFP             =  {0};
            %Grab Name and add T000# Format to such names
           if znum>1
                for z=1:znum
                    str_GFPname        =  Arr_GFP{K}{J}{j}{z};
                    cd                    (dir_GFP)
                    valsGFP{z}         =  imread(str_GFPname); 
                    cd                    (time_folder);
                    str_SaveName       =  strcat(str_Tfoldername,'-GFP-',str_GFPname);
                    imwrite               (valsGFP{z},str_SaveName,'Compression','none','Resolution',1);%,'BitDepth',16);
                end;
            else
                display('single Z for red')
                disp('no vals')
            end;
            disp                 ([num2str(100*j/(minfids),'%2.1f') ' % GFP' num2str(j) ' / ' num2str(minfids) ' ' num2str(J) ' out ' num2str(posnum) ' ' num2str(100*J/(posnum),'%2.1f')]);
        end;
    end;

end;
%% Save Zstacks TexRed Combine (Andre Transform)
%-Note Recently removed mkdir ---if it errors may be due to that.
cd                    (dir_ZirmiFormattedData);
mkdir                 (name);
cd                    (name)
ProcessedTemp_str   = cd; 
for J = 1:posnum 
    cd                   (ProcessedTemp_str) %Directory of Processed Experiment#
    Str_Posnum         = PositionName{J};
%     mkdir                (Str_Posnum)
    cd                   (Str_Posnum)
    pos_folder         = cd; %Processed data Position directory
    [pa3, name3, ext3] = fileparts(cd);
    counter            = 1000;
    cnter              = 0;
    disp                (strcat('On Position',num2str(J),'!!!!!!!!!!!!!!'))
    for K=1:parts
        if K==parts % --> This is so that we can put the PreInjury In its own folder
                Str_Posnum         = PositionName{J};
                Str_PreInjuryPosnum= strcat(Str_Posnum,'-PreInjury');
                cd                   (ProcessedTemp_str)            
%                 mkdir                (Str_PreInjuryPosnum)
                cd                   (Str_PreInjuryPosnum)
                pos_folder         = cd; %Processed data Position directory
        end;
        disp                        (strcat('On Part',num2str(K),'/ out of',num2str(parts))) %displaying progress
        minfids                 =   numfids(K);
        for j = 1:minfids
            %Count According to Part numfids cummulatively
            counter             =  counter+1; %IMPORTANT--COUNTING MECHANISM FOR FOLDERS accumulation
            cnter               =  cnter+1;
            cd                    (pos_folder) %Directory of Processed 'P#'
            disp                  (counter)
            %Change Folder Name to T000# format
            str_Tfoldername     =  strcat('T',num2str(counter));
            str_Tfoldername(2)  =  '0';
            %Create Folder with Such Name and go to Directory
%             mkdir                 (str_Tfoldername)   %GOING TO ASSUME AT THIS POINT FOLDER EXISTS                
            cd                    (str_Tfoldername)               
            [pa4, name4, ext4]  =  fileparts(cd);
            time_folder         =  cd;           
            valsTexRed             =  {0};
            %Grab Name and add T000# Format to such names
           if znum>1
               %PreFormat to determine dimensions
                str_TexRedname       =  Arr_TexRed{K}{J}{j}{1};
                cd                    (dir_TexRed)
                A                    =  imread(str_TexRedname); 
                [n,m]=size(A); % or [n,m]=size(tab2) because tab1 and tab2 have the same size
                E={};
                for z=1:znum
                    str_TexRedname        =  Arr_TexRed{K}{J}{j}{z};
                    cd                       (dir_TexRed)
                    valsTexRed{z}         =  imread(str_TexRedname); 
                    cd                       (time_folder);
                    str_SaveName          =  strcat(str_Tfoldername,'-TexRed-',str_TexRedname);
                    if z==1;
                        E= valsTexRed{z};
                    else
                        B= valsTexRed{z};
                        %find the maximum between tab1 and tab2 and then put the result into tab3
                        %row
                        for ii=1:n
                            %col
                            for jj=1:m
                                E(ii,jj)=max(E(ii,jj),B(ii,jj));
                            end
                        end
                    end;                             
                end;
            else
                disp('single Z for red')
                disp('no vals')
            end;
            imwrite              (E,str_SaveName,'Compression','none','Resolution',1);%,'BitDepth',16);
%             disp                 ([num2str(100*j/(minfids),'%2.1f') ' % TexRed' num2str(j) ' / ' num2str(minfids) ' ' num2str(J) ' out ' num2str(posnum) ' ' num2str(100*J/(posnum),'%2.1f')]);
        end;
    end;
end;
%% Save Iteration and Time Before
disp(strcat('Experiment:',name))
disp(strcat('Positions:',num2str(posnum)))
if posnum==8
    ti_d=65/60;
elseif posnum==9
    ti_d=74.3/60;
elseif posnum==10
    ti_d=85/60;
else
    ti_d=60/60;
end
disp            (strcat('Sampling Frequency :',...
                 num2str(ti_d),'(image per minute)'))
t_plateinput    =input('How many minutes post injury was first image taken? NOTE: Default is 30 MPI  ____:');
t_plate = t_plateinput;
disp            (strcat('Image #1 Taken  :',...
                 num2str(t_plate),'(minutes post injury)'))
%% Save
[pa2,n2,ext2]=fileparts(pa3);
cd(pa2);
% OriginalImg=valsPH2{1}{1};
str_savename=strcat(name,'.mat');
save(name,'name','ti_d','t_plate','numfids','PositionName','valsPh2','pa2','pa3','znum','posnum');
%disp for easy t
disp(name)
cd
dir
%% This is the end