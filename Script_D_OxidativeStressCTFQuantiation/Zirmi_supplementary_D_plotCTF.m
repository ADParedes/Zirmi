%% Supplementary Script D - Plotting CTF Values
% 2018-01-03 
% STEP 1:  Loading matifiles
% STEP 2:  Run this script
% Version 1.3 2018-01-03
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% Note to User
clc
disp ('Best Method to Loading Data For this Plot')
disp ('Dragging .mat file from directory and placing it in command window')
%% Define Variables
if exist('Zirmi') 
    boo_CTF                  = isfield(Zirmi.CTF,'zeroPixelValues');             % 1 means it exists
    switch boo_CTF
        case {1}
            disp                    ('CTF fish values are registered')            
            wound_CTCF              = (Zirmi.intDen - ((Zirmi.areaROI).*Zirmi.CTF.zeroPixelValues));
        otherwise
            disp                    ('CTF fish values are NOT registered')
            warning                 ('Zirmi metaData not registered/loaded for CTF Plotting Script')
            disp                    ('Load CTF metaData - Script Discontinued')
            return
    end;    
else
    warning                 ('Zirmi metaData not registered/loaded for CTF Plotting Script')
    disp                    ('Load CTF metaData - Script Discontinued')
    return
end
%% Parameters Defined
% T_p, t_plate ti_d , Time, wound_CTCF, golditeration, Img1, ch_PH2** special (maybe always 8)
T_p                         = Zirmi.Plot.idx_frames;
Img1                        = Zirmi.Plot.I_bf_1;
golditeration               = Zirmi.Plot.sp_frames;
Time                        = Zirmi.Plot.axisX;         % Convereted to axisX %arr_time                 
frames                      = Zirmi.Plot.num_frames;    % Number of frames (lens of frame)
MPI_start                   = PARAMETERS.Parameter5;    % Minutes post injury (Start of Imaging) ; old: t_plate
SamplingFrequency           = PARAMETERS.Parameter4;    % Sampling Frequency; old: ti_d  

%% Parameters for a Table   %MAXIMUM LENGTH IS 500 or ERROR
warning             ('Limit of 500 values ... this can be changed')
lens_time           = length(Time);
lens_add            = 500-lens_time;
arr_add             = zeros(lens_add,1);
arr_add_1           = zeros(30-(length(T_p)*2),1);

S1.str_col          ={strcat('CTF-',POI.Parameter10d)}; % wound_CTF    - Standardize Measure
S2.str_col          = S1.str_col; 
S3.str_col          ={'Time(MPI)',...
        strcat(POI.Parameter10d,'-IntDen')};        % wound_IntDen -Integrated Density
S4.str_col          ={'Time(MPI)',...
        strcat(POI.Parameter10d,'Area(pix^2)')};    % Area of Wound 
S5.str_col          ={'Time(MPI)',...
        strcat(POI.Parameter10d,'Background')};     % mmpio/MeanFluorescence - Background or Baseline Fluorescence

S1.num              = num2cell(vertcat([vs(T_p);Time(T_p)'],[arr_add_1]));               % vert   
S2.num              = num2cell(vertcat(vs,arr_add));                           % vert
S3.num              = num2cell(vertcat([Time',wound_IntDen'],[arr_add,arr_add]));         % vert
S4.num              = num2cell(vertcat([Time',PixelNum'],[arr_add,arr_add]));         % vert
S5.num              = num2cell(vertcat([Time',mmpi0'],[arr_add,arr_add]));         % vert

S1.table            = [S1.str_col;S1.num];
S2.table            = [S2.str_col;S2.num];
S3.table            = [S3.str_col;S3.num];
S4.table            = [S4.str_col;S4.num];
S5.table            = [S5.str_col;S5.num];




%% Presetting Plot Parameters
% frames=115;
if exist('wound_CTCF')
    axisY= wound_CTCF';
elseif exist('nowound_CTCF')
    axisY= nowound_CTCF';
else
    disp('You done made an error')
end;

axisX                = round(((1:frames)*MPI_start+SamplingFrequency)-MPI_start);

ymax   =  9e+07 ;
ymin   =  1.0e+06 ;
xmax    = axisX(end)+5;
xmin    = SamplingFrequency;
%%---------------------
figure(1);plot((axisX(1:frames)),axisY); % Advanced box plot
    ylim([ymin ymax])
    xlim([xmin xmax])
    xlabel(' Minutes after Injury'); % Set the X-axis label
    ylabel(' Pixel Intensity abu'); % Set the X-axis label    
%% Moving Average 
FigHandle = figure;
set(FigHandle, 'Position', [400, 100, 800, 650]); %[x1, y1, x2, y2]
plot((axisX(1:frames)),axisY,'--'); % Advanced box plot
    axis([xmin xmax ymin ymax]) 
hold on 
    %-- moving average plot addition
window_size = 4;
simple=tsmovavg(axisY,'s',window_size,1); %note sample needs to be vertical**
semi_gaussian = [0.026 0.045 0.071 0.1 0.12 0.138];
semi_gaussian = [semi_gaussian fliplr(semi_gaussian)];
weighted      = tsmovavg(axisY,'w',semi_gaussian,1); 
holder        = weighted(6:end)';
wholder       = cat(2,holder,nan(1,5));

plot((axisX(1:frames)),simple,'color','r')

hold on
   %-- moving average plot addition
plot((axisX(1:frames)),wholder', 'color','b')
    legend('Raw','Simple Moving Average','Weighted Moving Average',...
    'Location','NorthWest')
    set(gca,'fontweight','bold','fontsize',18)
%     title('ROS Generation in Wound','fontsize',12)
    xlabel(' Minutes after Injury','fontsize',18,'fontweight','bold'); % Set the X-axis label
    ylabel('Corrected Total Fluorescence (units)','fontsize',18,'fontweight','bold'); % Set the X-axis label
    grid on
if exist('T_p')
    disp('T_p values text values');
    for i=1:length(T_p)
        txt1=strcat(num2str(round(axisY(T_p(i)))),'');
        txt2='\uparrow';
        text(axisX(T_p(i))-1,(axisY(T_p(i))-500000),txt2)
        text(axisX(T_p(i)),(axisY(T_p(i))-500000),txt1)
    end;
else
    disp('golditeration text values')
    for i=1:4%length(T_p)
        txt1=strcat(num2str(round(axisY(golditeration(i)))),'');
        txt2='\uparrow';
        text(axisX(golditeration(i))-1,(axisY(golditeration(i))-500000),txt2)
        text(axisX(golditeration(i)),(axisY(golditeration(i))-500000),txt1)
    end;
end;

%% Picture of The Smelly Fish
figure()    
imagesc(Img1); colormap('gray')   
currAxPos = get(gca,'position');
if currAxPos(end) ==1
    set(gca,'position',[  0.1300    0.1150    0.7800    0.8100 ]);axis on
else
    set(gca,'position',[0 0 1 1 ]);axis off
end
clear currAxPos
text(10,10,strcat('viaZirmi'),'Color','red','FontSize',20)
