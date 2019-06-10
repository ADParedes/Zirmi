% 5/06/17 Evaluation of 5 metrics of interest 
% meandering index, velocity, static ratio, directionality and persistence
close all;
warning('Must Add path of Zirmi Script_E_Functions to Run this Script')
dir_desktop='C:\Users\AndreDaniel\Dropbox\Zirmi\Script_E_Functions';
addpath(dir_desktop)
%% Define/Select Parameters
clear s input_metric
micronsPerPixel             = PARAMETERS.Parameter2;   % LateralPixelResolution  
coe                         = micronsPerPixel/PARAMETERS.Parameter4;
str_all_timepoints          =['T1: 30- 59';'T2: 60- 89';'T3: 90-120';'T4: 30- 90';'T5: 30-120'];
str_int_timepoints          =['T1: 30- 59';'T2: 60- 89';'T3: 90-120'];
str_cum_timepoints          =['T4: 30- 90';'T5: 30-120';' * ALL *  '];
[N_ALL,~]                           = size(str_all_timepoints);
[N,~]                           = size(str_int_timepoints);
str_Metricname              ={'AbsVelocity','StaticRatio', 'MeanderingRatio',...
                                'Tortuosity','ForwardRatio','ForwardtoBackwardRatio'...
                                'WoundDistance','WoundScore'};
names                       = str_Metricname;
for ind                     = 1:length(str_Metricname)
    s.(names{ind})              = 1;
end;
 
screen                      = POI.Parameter12;
                            
rawMetrics_MetricType       ={SM.velocity,SM.staticratio,SM.meandering,...
                                SM.tortuosity,SM.forward,SM.FBratio,...
                                SM.WoundScoreUm,SM.WoundScore1234};   
%% Compute Metrics of interest 
disp('Please note Velocity units are micrometeres per minute')
input_metric                =input(strcat('Please Write Selected Metric e.g "velocity" or \n', ...
                                '"staticratio" or "meandering" or "tortuosity" or \n',...
                                '"forward","FtoB","FBratio","WoundScoreUm"'),'s');
    switch input_metric
        case{'velocity'}
            %%%--for velocity5
            t1                  = SM.(input_metric).interval{1}*coe;
            t2                  = SM.(input_metric).interval{2}*coe;
            t3                  = SM.(input_metric).interval{3}*coe;
            yy                  = SM.(input_metric).cum{1}*coe;
            t4                  = SM.(input_metric).cum{2}*coe;
            t5                  = SM.(input_metric).cum{3}*coe;
            yyy                 = SM.(input_metric).intervalsm{1}*coe;
            yyy1                = SM.(input_metric).intervalsm{2}*coe;
            yyy2                = SM.(input_metric).intervalsm{3}*coe;
            z                   = SM.(input_metric).phagosight*coe;
            PS                  = handles.distanceNetwork.absVelocity*coe;
            ymax    = max(t2)+1;
            ymin    = min(t2)-1;       
            input_metric ='Absolute Velocity (\mum/min)';
        case{'WoundScoreUm'}
            t1                  =SM.(input_metric).interval{1}*-1;
            t2                  =SM.(input_metric).interval{2}*-1;
            t3                  =SM.(input_metric).interval{3}*-1;
            yy                  =SM.(input_metric).cum{1}*-1; 
            t4                  =SM.(input_metric).cum{2}*-1;
            t5                  =SM.(input_metric).cum{3}*-1;
            yyy                 =SM.(input_metric).intervalsm{1}*-1;
            yyy1                =SM.(input_metric).intervalsm{2}*-1;
            yyy2                =SM.(input_metric).intervalsm{3}*-1;
            PS                  = handles.distanceNetwork.oriVelocity*88;
            z                   =handles.distanceNetwork.oriVelocity(Worthy.ID.cum)*88;    
            input_metric        ='Oriented Net Distance (\mum)';
            ymax                = max(t5)+0.2;
            ymin                = min(t1)-0.2;
        otherwise
            t1                   =SM.(input_metric).interval{1};
            t2                  =SM.(input_metric).interval{2};
            t3                  =SM.(input_metric).interval{3};
            yy                  =SM.(input_metric).cum{1};
            t4                 =SM.(input_metric).cum{2};
            t5                 =SM.(input_metric).cum{3};
            yyy                 =SM.(input_metric).intervalsm{1};
            yyy1                =SM.(input_metric).intervalsm{2};
            yyy2                =SM.(input_metric).intervalsm{3};
            
            switch input_metric
                case{'staticratio'}
                    PS                  = handles.distanceNetwork.staticRatio;
                    z                   =handles.distanceNetwork.staticRatio(Worthy.ID.cum);  
                case{'meandering'}
                    PS                  = handles.distanceNetwork.meanderRatio;
                    z                   =handles.distanceNetwork.meanderRatio(Worthy.ID.cum);   
                case{'tortuosity'}
                    PS                  = handles.distanceNetwork.tortuosity;
                    z                   =handles.distanceNetwork.tortuosity(Worthy.ID.cum);
                case{'forward'}
                    PS                  = handles.distanceNetwork.forwardRatioTot;
                    z                   =handles.distanceNetwork.forwardRatioTot(Worthy.ID.cum);
                otherwise
                    PS                  = NaN;
                    z                   = NaN;
            end;
            ymax    = max(t2)+0.2;
            ymin    = min(t2)-0.2;
    end;
%% Start new presentation
cd(ADP.dir_metadat)
mkdir('PowerPointviaZirmi')
cd('PowerPointviaZirmi')
 exportToPPTX('close'); %just in case
pptfileexist = strcat(POI.Parameter10d,'.pptx');
test_powerpoint = isempty(dir(pptfileexist));  

switch test_powerpoint
    case{1}  % It is empty
        isOpen  = exportToPPTX();
        if ~isempty(isOpen),
            % If PowerPoint already started, then close first and then open a new one
            exportToPPTX('close');
        end
        exportToPPTX('new','Dimensions',[12 6], ...
            'Title','Zirmi - notboxplot', ...
            'Author','Andre Daniel Paredes', ...
            'Subject','Automatically generated PPTX file', ...
            'Comments','This is for Zirmi users updated -2/13/2018');

        newFile = exportToPPTX('saveandclose',POI.Parameter10d);
        exportToPPTX('open',POI.Parameter10d);
    otherwise        
        exportToPPTX('open',POI.Parameter10d);
end;
disp('END: Cleaned a New PPTX')    
    
%% Statistics
npstat=(0);
pstat=(0);
%NOTE: h is only a logic signal to whether it hit the %5 statistical
%significance, it has nothing to do with parametric vs non parametric
%distribution
%     y                   =input_metric.interval{1};
%     y1                  =input_metric.interval{2};
%     y2                  =input_metric.interval{3};
%     yy                  =input_metric.cum{1};
%     yy1                 =input_metric.cum{2};
%     yy2                 =input_metric.cum{3};
%     z                   =input_metric.phagosight;
%Wilcoxon signed rank test.  - two sided nonparametric paired -p value
% [p,h,stats]    =signtest(yy1,z); %needs to be the same number of elements?!! what
% disp(p)
[p,h,stats]    =signrank(yy,t4); %needs to be the same number of elements?!! what    
disp([h p]); npstat(1)=p;
[p,h,stats]    =signrank(yy,t5); %needs to be the same number of elements?!! what    
disp([h p]); npstat(2)=p;
[p,h,stats]    =signrank(t4,t5); %needs to be the same number of elements?!! what 
disp([h p]); npstat(3)=p;

[p,h,stats]    =signrank(yyy,yyy1); %needs to be the same number of elements?!! what    
disp([h p]); npstat(4)=p;
[p,h,stats]    =signrank(yyy,yyy2); %needs to be the same number of elements?!! what    
disp([h p]); npstat(5)=p;
[p,h,stats]    =signrank(yyy1,yyy2); %needs to be the same number of elements?!! what 
disp([h p]); npstat(6)=p;



% tttest is for two sided parametric paired p value
%cum
[h,p,stats]    =ttest(yy,t4); %needs to be the same number of elements?!! what    
disp([h p]); pstat(1)=p;
[h,p,stats]    =ttest(yy,t5); %needs to be the same number of elements?!! what    
disp([h p]); pstat(2)=p;
[h,p,stats]    =ttest(t4,t5); %needs to be the same number of elements?!! what    
disp([h p]); pstat(3)=p;
%intervalSM
[h,p,stats]    =ttest(yyy,yyy1); %needs to be the same number of elements?!! what    
disp([h p]); pstat(4)=p;
[h,p,stats]    =ttest(yyy,yyy2); %needs to be the same number of elements?!! what    
disp([h p]); pstat(5)=p;
[h,p,stats]    =ttest(yyy1,yyy2); %needs to be the same number of elements?!! what    
disp([h p]); pstat(6)=p;
%phago

% ttest2 is for two sided parametric unpaired p value
% [h,p,stats]    =ttest2(y,y2); %needs to be the same number of elements?!! what    
% disp([h p])
% [h,p,stats]    =ttest2(yy1,z); %needs to be the same number of elements?!! what    
% disp([h p])
% [h,p,stats]    =ttest2(yy2,z); %needs to be the same number of elements?!! what    
% disp([h p])
% [h,p,stats]    =ttest2(yyy,z); %needs to be the same number of elements?!! what    
% disp([h p])
% [h,p,stats]    =ttest2(yyy1,z); %needs to be the same number of elements?!! what    
% disp([h p])
% [h,p,stats]    =ttest2(yyy2,z); %needs to be the same number of elements?!! what    
% disp([h p])
zstat       =[npstat;pstat];
yyys        =[yyy;yyy1;yyy2];     %intervalsm -- 70% per total track as defined by that interval (therefore needs to be in ALL intervals)
yys         =[yy;t4;t5];          %cum
% ys          =[y;y1;y2];             %interval  -- 90% per that interval
% We don't do this because it doesn't always equal the same stuff
zzz         = vertcat(yys,yyys);
% ys=[y;y1;y2];    
%% -ALL TIME DOMAINS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Slide 1
slideId     = exportToPPTX('addslide');

figure('Color',[1 1 1]);
set(gcf,'Name','Time Domains T1-T5','Position',[1 1 1100  420]);
ylim([ymin ymax]);
 notBoxPlot(t1,1,'style','sdline')
  hold on
 notBoxPlot(t2,2,'style','sdline')
  hold on
 notBoxPlot(t3,3,'style','sdline')
  hold on
 notBoxPlot(t4,4,'style','sdline')
   hold on
 notBoxPlot(t5,5,'style','sdline')
    set(gca,'XLim',[0 N_ALL+1],'XTick',1:N_ALL,'XTickLabel',str_all_timepoints,'fontweight','bold','fontsize',14);
     title(strcat('Tracks detectable in',{' '},num2str(PARAMETERS.ParameterA*100),'% of all frames'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
%   title('Meandering Index')
%   set(gca,'fontweight','bold','fontsize',18)
exportToPPTX('addpicture',gcf,'Position',[1 0.5 8 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Interval A','Position',[1 0 3 0.5],'Vert','bottom');
%% -Cumulative%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Color',[1 1 1]);ylim([ymin ymax])
%  plot(x,y,'x') 
%   hold on
%  plot(x+1,y1,'o')
%  plot(x+2,y2,'d')
%  plot(x+3,y3,'v')
 notBoxPlot(t4,1,'style','sdline')
  hold on

 notBoxPlot(t5,2,'style','sdline')
  hold on
 notBoxPlot(z,3,'style','sdline')
    set(gca,'XLim',[0 N+1],'XTick',1:N,'XTickLabel',str_cum_timepoints,'fontweight','bold','fontsize',14);
     title(strcat('Tracks detectable in',{' '},num2str(PARAMETERS.ParameterA*100),'% of all frames'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
% Top Right Corner
exportToPPTX('addpicture',gcf,'Position',[6 0.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Cumulative A','Position',[6 0 3 0.5],'Vert','bottom');
%% -Selected Metric Interval (intervalsm)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Color',[1 1 1]);ylim([ymin ymax])
%  plot(x,y,'x') 
%   hold on
%  plot(x+1,y1,'o')
%  plot(x+2,y2,'d')
%  plot(x+3,y3,'v')
 notBoxPlot(yyy,1,'style','sdline')
  hold on
 notBoxPlot(yyy1,2,'style','sdline')
  hold on
 notBoxPlot(yyy2,3,'style','sdline')
    set(gca,'XLim',[0 N+1],'XTick',1:N,'XTickLabel',str_int_timepoints,'fontweight','bold','fontsize',14);
     title(strcat('Tracks detectable in',{' '},num2str(POI.Parameter_gtA*100),'% of each time domain'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
% Bottom Left Corner
exportToPPTX('addpicture',gcf,'Position',[1 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Interval B','Position',[1 3 3 0.5],'Vert','bottom');

%% -Phagosight%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(z)
else
    %% Statistics
    %Wilcoxon signed rank test.  - two sided nonparametric paired -p value
    [p,h,stats]    =signrank(yy,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(7)=p;
    [p,h,stats]    =signrank(t4,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(8)=p;
    [p,h,stats]    =signrank(t5,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(9)=p;
    [p,h,stats]    =signrank(yyy,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(10)=p;
    [p,h,stats]    =signrank(yyy1,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(11)=p;
    [p,h,stats]    =signrank(yyy2,z); %needs to be the same number of elements?!! what    
    disp([h p]); npstat(12)=p;
    % ttest is for two sided parametric paired p value
    [h,p,stats]    =ttest(yy,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(7)=p;
    [h,p,stats]    =ttest(t4,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(8)=p;
    [h,p,stats]    =ttest(t5,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(9)=p;
    [h,p,stats]    =ttest(yyy,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(10)=p;
    [h,p,stats]    =ttest(yyy1,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(11)=p;
    [h,p,stats]    =ttest(yyy2,z); %needs to be the same number of elements?!! what    
    disp([h p]); pstat(12)=p;
    
    yys         =[yy;t4;t5;z];          %cum
    %% Plot
        figure('Color',[1 1 1]);ylim([ymin ymax]) 
            notBoxPlot(PS,1,'style','sdline')
        set(gca,'XLim',[0 2],'XTick',1,'XTickLabel','All frames','fontweight','bold','fontsize',14);
         title('\itPhagoSight\it all tracks')
        ylabel(input_metric,'fontsize',18,'fontweight','bold')
    % Bottom Right Corner
    exportToPPTX('addpicture',gcf,'Position',[6 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
%     exportToPPTX('addtext','Cumulative B','Position',[6 3 3 0.5],'Vert','bottom');        
end;

%% Clear Variables
%  clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2 time
cd(ADP.dir_metadat)
cd('PowerPointviaZirmi')
newFile = exportToPPTX('saveandclose',POI.Parameter10d);
fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
disp('Zirmi Cell Kinetics Outcome Measures exported to Powerpoint')
