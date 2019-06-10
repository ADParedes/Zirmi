function [ num_PreDf,num_PreBg,num_PreRw,num_PreMF,Img0,Tr0] = function_preInjuryCTF( dir_Xfolder,str_positionName,BPP,input40c )
%% function_preInjuryCTF Summary 
% This function is Used in Zirmi_D Scripts for CTF calculations
% If you want to include PreInjury Image data, you can use this. 
% NOTE: 40 circles is not used to determine preinjury ROS value
% This can be implemented 
%% Part 1: Perform Quantitation of ROS Baseline (PreInjury) at Region of Interest 
cd(dir_Xfolder)
dir_prename5            = strcat(str_positionName,'-PreInjury');
cd(dir_prename5)
dir_positionName               = cd;
dprename5           = dir('*T*');
boo_empty           = isempty(dprename5);
switch boo_empty
    case {1}
        warning ('There is no PreInjury Images')
        num_PreDf = NaN;
        num_PreBg = NaN;
        num_PreRw = NaN;
        num_PreMF = NaN;
        return
    otherwise
end;
    
len_d               = length(dprename5);
Arr_Pre_Bkgd        ={0};
Arr_Pre_Ifish       ={0};
Arr_Pre_Raw         ={0};
Arr_Pre_Difference  ={0};
Pre_Raw             =(0);
Pre_Bkgd            =(0);
Pre_Diff            =(0);
num_PreDf               = 0;
num_PreBg               = 0;
num_PreRw               = 0;
num_PreMF               = 0;
i_z                     = 1;
disp('Processing......')
disp(dir_prename5)
for i=1:1
    %% Set Parameters
    cd(dir_positionName)
    dir_time     = dprename5(i).name;
    cd(dir_time)
    str_TexRed  = strcat('*ex*');
    str_BF      = strcat('*Ph2*');
    str_GFP     = strcat('*GFP*');
    f_GFP       = dir(str_GFP);
    f_TexRed    = dir(str_TexRed);
    f_BF        = dir(str_BF);
    %% Select TexRed In Focus
        disp('Determining best z_stack');
        if length(f_TexRed)>1
            for ii=1:length(f_TexRed)
                figure(ii);imagesc(imread(f_TexRed(ii).name))
            end;
            i_z = input('Whats the best z-stack');
        else
            i_z = 1;
        end;
    everyother  = i;       
    str_OgImg = f_BF(1).name;
    seeMe       = f_TexRed(1).name; %We fused all zstacks together for TexRed, so should be the one and only
    %% Find Threshold for ROS
    OgImg=imread(str_OgImg);
    Img0=imadjust(OgImg);
    Tr0=imread(seeMe); 
        input_threshold=7000;
        input_happy='n';
        while input_happy ~= 'y'
            close all
            A=Tr0;
            A(A>=input_threshold)=1;
            A(A>1)=0;
            imagesc(A);
            save_Threshold=input_threshold;
            input_threshold=input(strcat('Set the Low Threshold, [',num2str(input_threshold),']'));
            if isempty(input_threshold)
                input_threshold=save_Threshold;
            else
            end;
            close all
            A(A>=input_threshold)=1;
            A(A<1)=0;
            imagesc(A)
            input_happy=input('Are you happy With Threshold','s');
        end;
        B=A;
        
        seeMe=f_TexRed(i_z).name;
        I=imread(seeMe);  %read image to beMe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%ROS Region
        if i==1 || i==3  %or statement
            disp('Select ROS Wound Region')
            close all; imagesc(Img0); title('ROS REGION');
            RosRegion=roipoly();
        else
        end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Background Region
        if i==1 
            disp('Select Background Region')
            close all; imagesc(Tr0); title('BACKGROUND REGION');
            BackgroundRegion = roipoly();
            close all;
        else
            disp(strcat('Keeping Previous Background Region on frame: ',num2str(i)))
        end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    arr_Pre_Bkgd    =(0);
    arr_Pre_Raw     =(0);
    arr_Pre_Diff    =(0);
    for j=1:length(f_TexRed)
        seeMe=f_TexRed(j).name;
        Tr0=imread(seeMe);  %read image to beMe
        BW=im2bw(Tr0,0.3);
        BW_clean1 = bwareaopen(BW, 2);
       %--Going to threshold for Notochord
        input_threshold=7000;
        input_happy='n';
        a=0;      
        while (isempty(input_happy)) || (input_happy ~= 'y')
            close all
            a=a+1;
            NC=Tr0;
            NC(NC>=input_threshold)=1;
            NC(NC>1)=0;   
            imagesc(NC)
            input_happy=input('Are you happy With Threshold','s');  
            if input_happy=='y'
                break
            elseif isempty(input_happy)
                break
            else                
                save_Threshold=input_threshold;
                input_threshold=input(strcat('Set the Low Threshold, [',num2str(input_threshold),']'));
                if isempty(input_threshold)
                    input_threshold=save_Threshold;
                else
                end;
            end;                               
        end;
     
        %--Threshold for Notochord
        se90 = strel('line', 3, 90); se0 = strel('line', 3, 0);            
        Bdil = imdilate(B, [se90 se0]); %figure, imshow(Adil), title('dilated gradient mask');
        BWdfill = imfill(Bdil, 'holes'); %figure, imshow(BWdfill); title('binary image with filled holes')             
        bw = bwareaopen(BWdfill, 5000);%Eliminating objects fewer than 5000 Pixels
        bw2=bwmorph(bw,'clean'); %clean not sure, can't hurt
        bw2=bwmorph(bw2,'thicken'); %thickens by just a little
        bw0=bw2;
        bwO=bwperim(bw0);
        bwO1=imdilate(bwO,[se90,se0]);
        bw_color_realsegout0=imoverlay(Tr0,bwO,[1 0 0]);%high is red
        %%%Converse is the Background
        se90bck = strel('line', 20, 90);
        se0bck = strel('line', 20, 0);
        bw_background = imdilate(BWdfill, [se90bck se0bck]);
        bw_background=~bw_background;
        %% USE 40 circles? 
        %input40c=input('Do 40 circles?','s');
        if input40c=='y'
%             disp('YES 40 circles')
            [MeanFluorescence0,Array_MF0]   = fortycircles(Tr0,bw0,10,1);  %10 stands for 10 pixels of diameter/1 stands for iterations;
        else
%             disp('NO 40 circles')
            MeanFluorescence0               = nan;
            Array_MF0                       = nan;
        end
        %Currently 40circles is not implemented for 
        %% To see the outline over the ROS Profile
        RR=RosRegion&bw;
        RR_outline=bwperim(RR);
        SegoutIntensity=Tr0;
        SegoutIntensity(RR_outline) = BPP;
        Fish = imadjust(Tr0, stretchlim(Tr0), []);
        I_fishRR=imoverlay(Fish,RR_outline,[1 0 0]);
        imshow(I_fishRR);  
        %% 
        RosRegionMask=RR;
        RawPixelIntensity=0;
        BackgroundPixelIntensity=0;
        CellPixelNum = sum(sum(RosRegionMask));
        BackgroundPixelNum = sum(sum(BackgroundRegion));%size(OriginalImgBF,1)*size(OriginalImgBF,2)-CellPixelNum; % or sum(sum(~CellMask))02/

        RawPixelIntensity = sum(sum(uint16(RosRegionMask).*Tr0))/CellPixelNum;
        BackgroundPixelIntensity = sum(sum(uint16(BackgroundRegion).*Tr0))/BackgroundPixelNum;
        DifferenceIntensity=[RawPixelIntensity-BackgroundPixelIntensity];
        arr_Pre_Bkgd(j)=BackgroundPixelIntensity;
        arr_Pre_Raw(j)=RawPixelIntensity;
        arr_Pre_Diff(j)=DifferenceIntensity;
        %(sum(sum(OriginalImgBF))-sum(sum(uint16(CellMask).*OriginalImgBF)))/BackgroundPixelNum;
    end;
        Arr_Pre_Bkgd{i}         = arr_Pre_Bkgd;
        Arr_Pre_Ifish{i}        = I_fishRR;
        Arr_Pre_Raw{i}          = arr_Pre_Raw;
        Arr_Pre_Difference{i}   = arr_Pre_Diff;
        Pre_Raw(i)              = arr_Pre_Raw(i_z);
        Pre_Bkgd(i)             = arr_Pre_Bkgd(i_z);
        Pre_Diff(i)             = arr_Pre_Diff(i_z);
        %     copyfile(combPh2img.name,combinePh2cd_s*tr)
end;   
%-Finalized Values 
num_PreDf = Pre_Diff;
num_PreBg = Pre_Bkgd;
num_PreRw = Pre_Raw;
num_PreMF = MeanFluorescence0;
Prex  =[num_PreDf,num_PreRw,num_PreBg];
%-Displaying
disp        (strcat('Background   : ',num2str(num_PreBg)))
disp        (strcat('Raw          : ',num2str(num_PreRw)))
disp        (strcat('Difference   : ',num2str(num_PreDf)))
disp        (strcat('40 Circles   : ',num2str(num_PreMF)))

end

