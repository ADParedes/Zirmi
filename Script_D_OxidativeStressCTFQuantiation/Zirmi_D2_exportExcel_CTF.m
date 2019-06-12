%% Supplementary Script D - Export To Excel
% STEP 1:  Loading matifiles
% STEP 2:  Run this script
% Version 1.4 2019-06-03
% Written By Andre Daniel Paredes | email @ andre.paredes@ymail.com
%% Note to User
clc
disp ('Best Method to Loading Data For this Plot')
disp ('Dragging .mat file from directory and placing it in command window')
%% Define Variables
if exist('Zirmi') 
    boo_CTF                  = isfield(Zirmi.Mock,'zeroPixelValues');             % 1 means it exists
    switch boo_CTF
        case {1}
            disp                    ('CTF fish values are registered')            
            Zirmi.CTF.CTF           = (Zirmi.CTF.intDen - ((Zirmi.CTF.areaROI).*Zirmi.Mock.zeroPixelValues));
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
idx                         = Zirmi.CTF.idx;
pixelValue                  = Zirmi.CTF.pixelValues - Zirmi.Mock.zeroPixelValues;
zeroPixelValue              = Zirmi.Mock.zeroPixelValues;
intDen                      = Zirmi.CTF.intDen;
woundArea                   = Zirmi.CTF.areaROI;
time                        = Zirmi.CTF.time;              
CTF_background              = Zirmi.CTF.CTF_background;
CTF_zirmi                   = Zirmi.CTF.CTF;

%% Parameters for a Table   %MAXIMUM LENGTH IS 500 or ERROR
warning             ('Limit of 500 values ... this can be changed')
lens_time           = length(time);
lens_add            = 500-lens_time;
arr_add             = zeros(lens_add,1);
arr_add_1           = zeros(30-(length(idx)*2),1);
str_col             ={POI.Parameter10d}; 
%% Selected Time Points (Golden Iterations)
S1.num              = num2cell(vertcat([pixelValue(idx)';time(idx)';arr_add_1]));               % vert   
S2.num              = num2cell(vertcat([zeroPixelValue(idx)';time(idx)';arr_add_1]));                             % vert
S3.num              = num2cell(vertcat([woundArea(idx)';time(idx)';arr_add_1]));  
S4.num              = num2cell(vertcat([intDen(idx)';time(idx)';arr_add_1])); 
S5.num              = num2cell(vertcat([CTF_background(idx)';time(idx)';arr_add_1])); 
S6.num              = num2cell(vertcat([CTF_zirmi(idx)';time(idx)';arr_add_1]));
str.sheets          = {'D-pixelValues','B-pixelValues','Area(pix^2)','Z-IntDen',...
                        'B-CTF','Z-CTF'};                       
C_head              = cell(1,length(str.sheets));
C_head(:)           = str_col;
C_body              =[S1.num,S2.num,S3.num,S4.num,S5.num,S6.num];
C_table             =[C_head;C_body];

%% EXCEL 
cd(ADP.dir_metadat)
[f,n,e]         = fileparts(ADP.dir_metadat);
cd              (f)
mkdir           ('Zirmi_CTF')
cd              ('Zirmi_CTF')
mkdir           ('viaZirmi_xlsx')
cd              ('viaZirmi_xlsx')
input_group      = input('What Group is This e.g. (Baseline,Control,3J,9J,18J...)','s');
file_excel       = strcat(input_group,'.xlsx');%DataType_Time
exist_excel      = exist(file_excel,'file');
            
switch exist_excel(1) 
        case{0} %file does not exist
            xlswrite(file_excel,[str_col;S1.num],'D-pixelValues')
            disp('pixelValues Archived')
            xlswrite(file_excel,[str_col;S2.num],'B-pixelValues')
            disp('zeroPixelValues Archived')
            xlswrite(file_excel,[str_col;S3.num],'Area(pix^2)')
            disp('Area Archived')
            xlswrite(file_excel,[str_col;S4.num],'Z-IntDen')
            disp('IntDen Archived')
            xlswrite(file_excel,[str_col;S5.num],'B-CTF')
            disp('CTF_background Archived')
            xlswrite(file_excel,[str_col;S6.num],'Z-CTF')
            disp('CTF_background Archived')
        otherwise %file does exist
            for i=1:length(str.sheets)
                [NUM_1,TXT_1,RAW_1] = xlsread(file_excel,str.sheets{i});
                empty_raw           = length(RAW_1);
                switch empty_raw
                    case {0,1} % EMPTY
                        disp('FYI Its Empty')
                        xlswrite(file_excel,C_table(:,1),str.sheets{i})
                        disp('S-CTF Updated')

                    otherwise  % NOT EMPTY
                        new_S1  = [RAW_1,C_table(:,1)];
                        xlswrite(file_excel,new_S1,str.sheets{i})
                        disp(str.sheets{i})

                end;
            end;
end;
winopen(file_excel)


