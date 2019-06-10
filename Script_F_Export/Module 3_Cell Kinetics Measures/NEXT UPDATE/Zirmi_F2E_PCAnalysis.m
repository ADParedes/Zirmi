%% Zirmi Package - Supplmentary PCA Analysis
% Excel Sheets contain a specific format
% This is a script read zirmi outputed excel sheets
% Script PCA ByTime version 01
% File      : User choses
% Update    : 2018-01-09
% exports   : SuperParameterFile.pptx
% By        : ADParedes 
% For Who   : OpenSource
% email     : andre.paredes@ymail.com
%% Note to User
% This Can do By Time point
%% Clean & Notify User
clc
clear                   all %#ok<CLALL>
close                   all
disp                    ('Start Principle Component Analysis') %Includes Pretest
disp                    ('***Ignores subjects with empty/NaN variable values***')
%% Define Array Parameters
%---Strings that will be needed for sorting
clear               str
str.Dose            = {'Baseline','0J','3J','9J','18J'};     
str.Variables       = {'Velocity(um/min)','StaticRatio','MeanderingRatio'}%,'Tortuosity(abu)',...
                       % 'ForwardRatio','ForwardtoBackwardRatio','WoundDistance(um)'};
str.timePoint       = {'60 min','90 min','120 min'};
str.Meta            = 'Unique ExpID';       
%---Numbers
num.Dose            = {-1,0,3,9,18};
%---Files to open (.xlsx files)
clear               filename
input_filesearch    = input('Type Choice as is: "cum","int","SM"','s');
filesearch.cum      ={'Cumulative*60*','Cumulative*90*','Cumulative*120*'};
filesearch.int      ={'Interval_GT*60*','Interval_GT*90*','Interval_GT*120*'};
filesearch.sm       ={'*SM*60*','*SM*90*','*SM*120*'};
switch input_filesearch
    case{'cum'}
        file_search  = filesearch.cum ;
    case{'int'}
        file_search  = filesearch.int;
    otherwise
        file_search  = filesearch.sm;
end;

filesearch.xls      =strcat('Zirmi_PCAnalysis_',input_filesearch,'.xlsx');
filesearch.ppt      =strcat('Zirmi_PCAnalysis_',input_filesearch);
%% Recognize User
hostname                    =   char( getHostName( java.net.InetAddress.getLocalHost ) );
ip                          =   char( getHostAddress( java.net.InetAddress.getLocalHost ) );
user                        =   getenv('UserName');    
clear               doi

if strcmp(hostname,'ADP-Chicago'); %this is the desktop
    disp('Did NOT select Labtop Directory')
    
    
elseif strcmp(hostname,'ADPAREDES') %this is the labtop
    disp('Selected Labtop Directory')
    boo_desktop             =   0;
    doi.main            = 'C:\Users\AndreDaniel\OneDrive - University of Illinois at Chicago\PhD Work\PhD Data\Data_Manuscript2\RAW';
    doi.save            = 'C:\Users\AndreDaniel\OneDrive - University of Illinois at Chicago\PhD Work\PhD Data\Data_Manuscript2\viaMatlab';
%---
    doi.functions       = 'C:\Users\AndreDaniel\Documents\Dropbox\Zirmi\Script_E_Functions';
else
    warning('Selected nothing');
    boo_destop              =   2;
    doi.main             = uigetdir('Select Directory where you keep excel files');
    doi.save            = uigetdir('Select Directory where you would like to export files');
    doi.functions       = uigetdir('Select Directory where you keep Zirmi Functions');
%---
end;
addpath(doi.functions);
disp('End: Selected Directory of Interest') 
%% Start new presentation
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[12 6], ...
    'Title','BARDA REPORT - Principle Component Analysis', ...
    'Author','Andre Daniel Paredes', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This is for Barda report 2017/2018 - Bartholomew Lab');
cd(doi.save)
newFile = exportToPPTX('saveandclose',filesearch.ppt);
exportToPPTX('open',filesearch.ppt);
disp('END: Cleaned a New PPTX')
%% Loop Determinants
lens_1                  = numel(str.Dose);       % 5 folders
lens_2                  = numel(file_search);    % 3 excel files
lens_3                  = numel(str.Variables);  % 7 Variables  - Lets customize this

Cell_timePoint               ={};
%---Loop Time Point -3- (e.g. 60,90,120)
for K                   = 1:lens_2
    %---Select Directory
    cd                      (doi.main)
    disp                    (str.timePoint{K})
    %---Define Variables
    arr_body                ={};
    arr_left                ={};
    super_body              ={};
    %% Loop Group (e.g. 0,3,9,18 J)
    for k                   = 1:lens_1
        %% Define Variables
        dirs                    = 0;
        filename                = 0;
        %---str
        str_header              ={};    
        str_observation_0       ={};
        str_dose_0              ={};
        str_variables           ={};
        %---num
        num_body_0              ={};
        num_dose_0              ={};
        %---concat arrays
        arr_left_0              ={};
        arr_body_0              ={};
        %---Final Array
        super_body_0            ={};        
        %% Select Directory 
        cd                      (doi.main)
        cd                      (str.Dose{k});  
        disp                    (str.Dose{k});
        %---Select File and Load Array
        dirs                    = dir(file_search{K});
        filename                = dirs.name;
        [enum,etxt,eraw]        = xlsread(filename,1);
        %% Deliniate Strings of interest
        str_header              = etxt(1,:);
            temp_observation    = strfind(str_header,str.Meta);
            idx_observation              = find(not(cellfun('isempty',temp_observation)));
        str_observations_0        = etxt(2:end,idx_observation);
            len_observations        = numel(str_observations_0);
        str_dose_0                = repmat(str.Dose(k),len_observations,1);
            temp_variables1     = strfind(str_header,str.Variables{1});
            temp_variables2     = strfind(str_header,str.Variables{end});
            idx_variable1       = find(not(cellfun('isempty',temp_variables1)));
            idx_variable2       = find(not(cellfun('isempty',temp_variables2)));
        str_variables           = str_header(idx_variable1:idx_variable2);
        arr_left_0                = [str_observations_0,str_dose_0];
        %% Deliniate Numbers of Interest
            temp_num_dose              = ones(len_observations,1)*num.Dose{k};
        num_dose_0              = num2cell(temp_num_dose);
        num_body_0              = enum(:,idx_variable1:idx_variable2);
        arr_body_0              = num2cell(num_body_0);
        %% Concat for Super Array
        super_body_0            = [str_observations_0,str_dose_0,arr_body_0];
        %% Super Concat
        switch k 
            case{1}
                num_body        = num_body_0;
                num_dose        = num_dose_0;
                str_dose        = str_dose_0;
                arr_body        = arr_body_0;
                arr_left        = arr_left_0;
                super_body      = super_body_0;
            otherwise
                num_body        = vertcat(num_body,num_body_0);
                num_dose        = vertcat(num_dose,num_dose_0);
                str_dose        = vertcat(str_dose,str_dose_0);
                arr_body        = vertcat(arr_body,arr_body_0);
                arr_left        = vertcat(arr_left,arr_left_0);
                super_body      = vertcat(super_body,super_body_0);
        end;
%         disp                    (filename)
    end;
    cd(doi.save)
    %% Save .xlsx file (Excel)
%     Cell_timePoint{K}    = super_body;
%     xlswrite             (filesearch.xls, Cell_timePoint{K},str.timePoint{K},'A1')
%     fullfilename         = fullfile(doi.save,filesearch.xls);
%     fprintf             ('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',...
%                             fullfilename,fullfilename);
    %% Save .pptx file (Powerpoint)
    input_als               ='n'; %NO ALS  
    str.Groups              = str.Dose;
    arr.Group               = cell2mat(num_dose);
    num.Groups              = num.Dose;
    str.var1                = str_variables;
    function_zirmi_pca     (num_body,str.timePoint{K},arr,num,str,input_als);

            

end;

disp('END: Loop Madness')
%% Finished
cd(doi.save)
newFile = exportToPPTX('saveandclose',filesearch.ppt);
fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);

disp('FIN')
