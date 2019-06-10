%% Supplementary Script D - Plotting CTF Values
% 2018-01-03 
% STEP 1:  Loading matifiles
% STEP 2:  Run this script
% Version 1.3 2018-01-03
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
close all
%% Start new presentation
cd              (ADP.dir_metadat)
[f,n,e]         = fileparts(ADP.dir_metadat);
cd              (f)
mkdir           ('Zirmi_CTF')
cd              ('Zirmi_CTF')
mkdir           ('PowerPointviaZirmi')
cd              ('PowerPointviaZirmi')
exportToPPTX    ('close'); %just in case
pptfileexist    = strcat(POI.Parameter10d,'.pptx');
test_powerpoint = isempty(dir(pptfileexist));  
%---------------New Powerpoint Check
switch          test_powerpoint
    case        {1}  % It is empty
        isOpen      = exportToPPTX();
        if          ~isempty(isOpen),
            % If PowerPoint already started, then close first and then open a new one
            exportToPPTX('close');
        else
        end
        exportToPPTX('new','Dimensions',[12 6], ...
            'Title','Zirmi - notboxplot', ...
            'Author','Andre Daniel Paredes', ...
            'Subject','Automatically generated PPTX file', ...
            'Comments','This is for Zirmi users updated -2/13/2018');

        newFile     = exportToPPTX('saveandclose',POI.Parameter10d);
        exportToPPTX('open',POI.Parameter10d);
    otherwise        
        exportToPPTX('open',POI.Parameter10d);
end;
disp('END: Cleaned a New PPTX')         
%% Presetting Plot Parameters
% frames=115;
if exist('wound_CTCF')
    vs           = wound_CTCF';
elseif exist('nowound_CTCF')
    vs           = nowound_CTCF';
    wound_IntDen = nowound_IntDen;
    T_p          = golditeration(1:4);
    PixelNum     = ixelNum;
else
    disp('You done made an error')
end;

Time                = round(((1:frames)*ti_d+t_plate)-ti_d);

ymax   =  10e+10 ;
ymin   =  10e+2;
xmax    = Time(end)+5;
xmin    = t_plate;
%% Figure 1 
slideId     = exportToPPTX('addslide');

%% Moving Average 

FigHandle = figure;
grid on
set(gcf,'color','w');
set(FigHandle, 'Position', [400, 100, 800, 650]); %[x1, y1, x2, y2]
plot((Time(1:frames)),vs,'--'); % Advanced box plot
    axis([xmin xmax ymin ymax]) 
hold on 
    %-- moving average plot addition
window_size = 4;
simple=tsmovavg(vs,'s',window_size,1); %note sample needs to be vertical**
semi_gaussian = [0.026 0.045 0.071 0.1 0.12 0.138];
semi_gaussian = [semi_gaussian fliplr(semi_gaussian)];
weighted      = tsmovavg(vs,'w',semi_gaussian,1); 
holder        = weighted(6:end)';
wholder       = cat(2,holder,nan(1,5));

plot((Time(1:frames)),simple,'color','r')

hold on
   %-- moving average plot addition
plot((Time(1:frames)),wholder', 'color','b')
    legend('Raw','Simple Moving Average','Weighted Moving Average',...
    'Location','NorthEast')
    set(gca,'fontweight','bold','fontsize',18)
    set(gca, 'YScale', 'log')
%     title('ROS Generation in Wound','fontsize',12)
    xlabel(' Minutes Post Injury','fontsize',24,'fontweight','bold'); % Set the X-axis label
    ylabel('Corrected Total Fluorescence (units)','fontsize',24,'fontweight','bold'); % Set the X-axis label
    grid on
if exist('T_p')
    disp('T_p values text values');
    for i=1:length(T_p)
        txt1={'\downarrow',strcat(num2str(round(vs(T_p(i)))),'')};        
        text(Time(T_p(i))-1,(vs(T_p(i))-500000),txt1,'FontSize',16)
%         txt1=strcat(num2str(round(vs(T_p(i)))),'');
%         txt2='\uparrow';
%         text(Time(T_p(i))-1,(vs(T_p(i))-500000),txt2)
%         text(Time(T_p(i)),(vs(T_p(i))-500000),txt1)
    end;
else
    disp('golditeration text values')
     %golditeration instead of T_p
    for i=1:4%length(T_p)
        txt1={'\downarrow',strcat(num2str(round(vs(golditeration(i)))),'')};        
        text(Time(golditeration(i))-1,(vs(golditeration(i))-500000),txt1,'FontSize',16)

    end;
end;
exportToPPTX('addpicture',gcf,'Position',[1 0.5 8 6]);
%% Figure 2
figure()    
imagesc(Img1); 
colormap(gray.^(1.5))  
currAxPos         = get(gca,'position');
if currAxPos(end) ==1
    set                 (gca,'position',[  0.1300    0.1150    0.7800    0.8100 ]);axis on
else
    set                 (gca,'position',[0 0 1 1 ]);axis off
end
clear currAxPos
% text(10,10,strcat(name,name5),'Color','red','FontSize',20)
t=text             (200,480,'via Zirmi','Color','black','FontSize',32);
%%       
exportToPPTX     ('addpicture',gcf,'Position',[6 3.5 4 3])
 
 %% Clear Variables
%  clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 time
cd(ADP.dir_metadat)
cd('PowerPointviaZirmi')
newFile = exportToPPTX('saveandclose',POI.Parameter10d);
fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
disp('Zirmi Oxidative Stress Outcome Measures exported to Powerpoint')


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

%% EXCEL 
cd(ADP.dir_metadat)
[f,n,e] = fileparts(ADP.dir_metadat);
cd(f)
mkdir('Zirmi_CTF')
cd('Zirmi_CTF')
mkdir('ExcelviaZirmi')
cd('ExcelviaZirmi')
input_group         = input('What Group is This e.g. (Baseline,Control,3J,9J,18J...)','s');
file_excel         = strcat(input_group,'.xlsx');%DataType_Time

exist_excel        = exist(file_excel,'file');
            
switch exist_excel(1) 
        case{0} %file does not exist
            xlswrite(file_excel,S1.table,'S-CTF')
            disp('S-CTF Archived')
%             xlswrite(file_excel,S2.table,'CTF')
%             disp('CTF Archived')
%             xlswrite(file_excel,S3.table,'IntDen')
%             disp('IntDen Archived')
%             xlswrite(file_excel,S4.table,'Area')
%             disp('Area Archived')
%             xlswrite(file_excel,S5.table,'Background')
%             disp('Background Archived')
        otherwise %file does exist
            [NUM_1,TXT_1,RAW_1] = xlsread(file_excel,1);
            empty_raw           = length(RAW_1);
            switch empty_raw
                case {0,1} % EMPTY
                    disp('FYI Its Empty')
                    xlswrite(file_excel,S1.table,'S-CTF')
                    disp('S-CTF Updated')

                otherwise  % NOT EMPTY
                    new_S1  = [RAW_1,S1.table];
                    xlswrite(file_excel,new_S1,'S-CTF')
                    disp('S-CTF Updated')
                    return
            end;
end;
%%
winopen(file_excel)




