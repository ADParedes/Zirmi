%% Heatmap Figure 3 Preperation
%  clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 time
close all
clc
%% Hashtag Parameters
str_type        ={'interval','intervalsm','cum'};

str_metric      ={'velocity','staticratio', 'meandering',... %1-3
                                'tortuosity','forward','FBratio'...%4-6
                                'WoundScoreUm','WoundScore1234',...%7-9
                                'WoundStartUm'};                   %9   
m               = 1;
type            = 2;

Temporal         ={'T1','T2','T3'};
Time             = [1,2,3];
str_y           ={'Time Domain'};
Spatial        ={'S1','S2','S3','S4'};
str_x           ={'Space Domain'};
%% Determine Distance from Wound Start Position (same as score)     
d               =[SM.(str_metric{9}).(str_type{type}){1}',SM.(str_metric{9}).(str_type{type}){2}',...
                                SM.(str_metric{9}).(str_type{type}){3}'];                
[row,col]       = size(d);
clear sorted_d_1 sorted_indices sorted_id sorted_score_1

for i=1:col
    A                   = d(:,i);
    [sorted, indices ]  = sort(A);
    sorted_d_1{i}         = sorted;
    sorted_indices{i}   = indices;
    sorted_id{i}        = Worthy.ID.(str_type{type})(indices)';
    sorted_score_1{i}     = SM.WoundScore1234.(str_type{type}){i}(indices)'; 

end; 
d_id            = horzcat(d,Worthy.ID.(str_type{type})');
sorted_d        = horzcat(sorted_d_1{:});
sorted_indices  = horzcat(sorted_indices{:}); %Is what we use to get the proper order of the 

sorted_id       = horzcat(sorted_id{:});

sorted_d_id     = horzcat(sorted_d,sorted_id);

sorted_score_1    = horzcat(sorted_score_1{:});
%% Loop to do all metrics

    metric          =[SM.(str_metric{m}).(str_type{type}){1}',SM.(str_metric{m}).(str_type{type}){2}',...
                        SM.(str_metric{m}).(str_type{type}){3}'];
[row,col]       = size(d);                   
    for i=1:col
        B                      = metric(:,i);
        sorted_metric_1{i}       = B(sorted_indices(:,i));
        sorted_id_checkB{i}    = Worthy.ID.(str_type{type})(sorted_indices(:,i));
        C                      = cell(1,row);
        C(:)                   = Temporal(i);
        D{i}                   = ones(1,row)'*Time(i);
        sorted_temporal_1{i}     = C';
        sorted_score{i}     = SM.WoundScore1234.(str_type{type}){i}(sorted_indices(:,i))'; 
    end;

Distance            = vertcat(sorted_d_1{:});
Metric              = vertcat(sorted_metric_1{:});
SpatialDomain       = vertcat(sorted_score{:});
TemporalDomain      = vertcat(sorted_temporal_1{:});
TimeDomain          = vertcat(D{:});
T                   = table(SpatialDomain,TimeDomain,Distance,Metric);   
% tb1                 = readtable(T);
% 
%  h = HeatMap(T,'TimeDomain','SpatialDomain','ColorVariable','Metric');

%% heatmap
minValue        = min(min(Metric));
maxValue        = max(max(Metric));
clims           = [minValue maxValue];
Metric_1            = horzcat(sorted_metric_1{:});
SpatialDomain_1     = horzcat(sorted_score{:});
TimeDomain          = horzcat(D{:});
heatmap         = Metric_1';
hmo             = HeatMap(heatmap,'Colormap','jet','DisplayRange',round(maxValue),... %round(maxValue)
                              'Annotate','on', 'AnnotColor','k',...
                              'LabelsWithMarkers',true);
         tits   = addTitle(hmo,strcat(str_metric{m},' - ',str_type{type}));
         Hx     = addXLabel(hmo,'Spatial');%   - Add HeatMap x-axis (column) label.
         Hy     = addYLabel(hmo,'Temporal');%   - Add HeatMap y-axis (row) label.
 close();
%% plot hetamap
h               = plot(hmo);         
         tits   = addTitle(hmo,strcat(str_metric{m},' -  ',str_type{type}));
         Hx     = addXLabel(hmo,'Spatial');%   - Add HeatMap x-axis (column) label.
         Hy     = addYLabel(hmo,'Temporal');%   - Add HeatMap y-axis (row) label.
%%
