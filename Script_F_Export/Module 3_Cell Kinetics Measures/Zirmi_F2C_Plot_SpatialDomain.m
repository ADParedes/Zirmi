% 5/06/17 Evaluation of 5 metrics of interest 
% meandering index, velocity, static ratio, directionality and persistence
close all;
warning('Must Add path of Zirmi Script_E_Functions to Run this Script')
warning('Must have ran S1A TimeDomain and Excel Export first')
dir_desktop='C:\Users\AndreDaniel\Dropbox\Zirmi\Script_E_Functions';
addpath(dir_desktop)
%% Define/Select Parameters
clear s input_metric
micronsPerPixel             = PARAMETERS.Parameter2;   % LateralPixelResolution  
coe                         = micronsPerPixel/PARAMETERS.Parameter4;
str_int_timepoints          =['S1';'S2';'S3';'S4'];
str_cum_timepoints          =['S1';'S2';'S3';'S4'];
[N,~]                           = size(str_int_timepoints);
str_Metricname              ={'AbsVelocity','StaticRatio', 'MeanderingRatio',...
                                'Tortuosity','ForwardRatio','ForwardtoBackwardRatio'...
                                'WoundDistance','WoundScore'};
names                       = str_Metricname;
for ind                     = 1:length(str_Metricname)
    s.(names{ind})              = 1;
end;
                            
                            
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
            t1_s1                =INT{1}.xls.Score1.vals(:,2)*coe; 
            t1_s2                 =INT{1}.xls.Score2.vals(:,2)*coe;
            t1_s3                  =INT{1}.xls.Score3.vals(:,2)*coe;
            t1_s4                   =INT{1}.xls.Score4.vals(:,2)*coe; 
            t3_s1                 =INT{3}.xls.Score1.vals(:,2)*coe; 
            t3_s2                 =INT{3}.xls.Score2.vals(:,2)*coe; 
            t3_s3                 =INT{3}.xls.Score3.vals(:,2)*coe; 
            t3_s4                =INT{3}.xls.Score4.vals(:,2)*coe; 
            t5_s1                =CUM{3}.xls.Score1.vals(:,2)*coe; 
            t5_s2                =CUM{3}.xls.Score2.vals(:,2)*coe; 
            t5_s3                =CUM{3}.xls.Score3.vals(:,2)*coe;
            t5_s4                =CUM{3}.xls.Score4.vals(:,2)*coe;
            z                   = SM.(input_metric).phagosight*coe;
            PS                  = handles.distanceNetwork.absVelocity*coe;
            ymax    = max(t1_s2)+1;
            ymin    = min(t1_s2)-1;       
            input_metric ='Absolute Velocity (\mum/min)';
        case{'WoundScoreUm'}

            t1_s1                =INT{1}.xls.Score1.vals(:,8)*-1; 
            t1_s2                 =INT{1}.xls.Score2.vals(:,8)*-1;
            t1_s3                  =INT{1}.xls.Score3.vals(:,8)*-1;
            t1_s4                   =INT{1}.xls.Score4.vals(:,8)*-1; 
            t3_s1                 =INT{3}.xls.Score1.vals(:,8)*-1; 
            t3_s2                 =INT{3}.xls.Score2.vals(:,8)*-1; 
            t3_s3                 =INT{3}.xls.Score3.vals(:,8)*-1; 
            t3_s4                =INT{3}.xls.Score4.vals(:,8)*-1; 
            t5_s1                =CUM{3}.xls.Score1.vals(:,8)*-1; 
            t5_s2                =CUM{3}.xls.Score2.vals(:,8)*-1; 
            t5_s3                =CUM{3}.xls.Score3.vals(:,8)*-1;
            t5_s4                =CUM{3}.xls.Score4.vals(:,8)*-1;
            PS                  = WoundScoreUm.cum{1}*-1;
            z                   =WoundScoreUm.cum{1}(Worthy.ID.cum)*-1;    
            input_metric        ='Oriented Net Distance (\mum)';
            allinclusive        = [t1_s1;t1_s2;t1_s3;t1_s4;...
                                    t3_s1;t3_s2;t3_s3;t3_s4;...
                                    t5_s1;t5_s2;t5_s3;t5_s4;...
                                    PS'];
            ymax                = max(allinclusive)+5;
            ymin                = min(allinclusive)-5;
       otherwise            
            switch input_metric
                case{'staticratio'}
                    PS                  = handles.distanceNetwork.staticRatio;
                    z                   = handles.distanceNetwork.staticRatio(Worthy.ID.cum);  
                    cc                  = 3;
                case{'meandering'}
                    PS                  = handles.distanceNetwork.meanderRatio;
                    z                   =handles.distanceNetwork.meanderRatio(Worthy.ID.cum); 
                    cc                  = 4;
                case{'tortuosity'}
                    PS                  = handles.distanceNetwork.tortuosity;
                    z                   = handles.distanceNetwork.tortuosity(Worthy.ID.cum);
                    cc                  = 5;
                case{'forward'}
                    PS                  = handles.distanceNetwork.forwardRatioTot;
                    z                   = handles.distanceNetwork.forwardRatioTot(Worthy.ID.cum);
                    cc                  = 6;
                case{'FtoB'}
                    PS                  = NaN;
                    z                   = NaN;
                    cc                  = 7;
                otherwise
                    PS                  = NaN;
                    z                   = NaN;
                    cc                  = 1;
                    warning('You did not type a measure correctly, TRY AGAIN')
            end;
            t1_s1                =INT{1}.xls.Score1.vals(:,cc); 
            t1_s2                 =INT{1}.xls.Score2.vals(:,cc);
            t1_s3                  =INT{1}.xls.Score3.vals(:,cc);
            t1_s4                   =INT{1}.xls.Score4.vals(:,cc); 
            t3_s1                 =INT{3}.xls.Score1.vals(:,cc); 
            t3_s2                 =INT{3}.xls.Score2.vals(:,cc); 
            t3_s3                 =INT{3}.xls.Score3.vals(:,cc); 
            t3_s4                =INT{4}.xls.Score4.vals(:,cc); 
            t5_s1                =CUM{3}.xls.Score1.vals(:,cc); 
            t5_s2                =CUM{3}.xls.Score2.vals(:,cc); 
            t5_s3                =CUM{3}.xls.Score3.vals(:,cc);
            t5_s4                =CUM{3}.xls.Score4.vals(:,cc);
            ymax    = max(t1_s2)+0.2;
            ymin    = min(t1_s2)-0.2;
    end;
%% Start new presentation
cd(ADP.dir_metadat)
[filepath,name,ext] = fileparts(ADP.dir_metadat)
mkdir('PowerPointviaZirmi')
cd('PowerPointviaZirmi')
 exportToPPTX('close'); %just in case
% isOpen  = exportToPPTX();
% if ~isempty(isOpen),
%     % If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
% exportToPPTX('new','Dimensions',[12 6], ...
%     'Title','Zirmi - notboxplot', ...
%     'Author','Andre Daniel Paredes', ...
%     'Subject','Automatically generated PPTX file', ...
%     'Comments','This is for Zirmi users updated -2/13/2018');
% newFile = exportToPPTX('saveandclose',POI.Parameter10d);

exportToPPTX('open',POI.Parameter10d);
disp('END: Cleaned a New PPTX')    
    
    
%% -T1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Slide 1
slideId     = exportToPPTX('addslide');
slideId     = exportToPPTX('addslide');
figure('Color',[1 1 1]);ylim([ymin ymax])

 notBoxPlot(t1_s1,1,'style','sdline')
    hold on
 notBoxPlot(t1_s2,2,'style','sdline')
    hold on
 notBoxPlot(t1_s3,3,'style','sdline')
    hold on
 notBoxPlot(t1_s4,4,'style','sdline')
    set(gca,'XLim',[0 N+1],'XTick',1:N,'XTickLabel',str_int_timepoints,'fontweight','bold','fontsize',14);
     title(strcat('Tracks detectable in',{' '},num2str(POI.ParameterA*100),'% of T1'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
%   title('Meandering Index')
%   set(gca,'fontweight','bold','fontsize',18)
exportToPPTX('addpicture',gcf,'Position',[1 0.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Interval A','Position',[1 0 3 0.5],'Vert','bottom');
%% -T3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color',[1 1 1]);ylim([ymin ymax])
 notBoxPlot(t3_s1,1,'style','sdline')
    hold on
 notBoxPlot(t3_s2,2,'style','sdline')
    hold on
 notBoxPlot(t3_s3,3,'style','sdline')
    hold on
 notBoxPlot(t3_s4,4,'style','sdline')
    set(gca,'XLim',[0 N+1],'XTick',1:N,'XTickLabel',str_cum_timepoints,'fontweight','bold','fontsize',14);
      title(strcat('Tracks detectable in',{' '},num2str(POI.ParameterA*100),'% of T3'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
% Top Right Corner
exportToPPTX('addpicture',gcf,'Position',[6 0.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Cumulative A','Position',[6 0 3 0.5],'Vert','bottom');
%% -T5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color',[1 1 1]);ylim([ymin ymax])
 notBoxPlot(t5_s1,1,'style','sdline')
    hold on
 notBoxPlot(t5_s2,2,'style','sdline')
    hold on
 notBoxPlot(t5_s3,3,'style','sdline')
    hold on
 notBoxPlot(t5_s4,4,'style','sdline')
    set(gca,'XLim',[0 N+1],'XTick',1:N,'XTickLabel',str_cum_timepoints,'fontweight','bold','fontsize',14);
      title(strcat('Tracks detectable in',{' '},num2str(POI.ParameterA*100),'% of T5'))
    ylabel(input_metric,'fontsize',18,'fontweight','bold')
% Bottom Left Corner
exportToPPTX('addpicture',gcf,'Position',[1 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
% exportToPPTX('addtext','Interval B','Position',[1 3 3 0.5],'Vert','bottom');
%% Statistics
npstat=(0);
pstat=(0);
%NOTE: h is only a logic signal to whether it hit the %5 statistical
%significance, it has nothing to do with parametric vs non parametric
%distribution

%% -Phagosight%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(z)
else
    %% Statistics
    %Wilcoxon signed rank test.  - two sided nonparametric paired -p value

    %% Plot
        figure('Color',[1 1 1]);ylim([ymin ymax]) 
            notBoxPlot(PS,1,'style','sdline')
        set(gca,'XLim',[0 2],'XTick',1,'XTickLabel','All frames','fontweight','bold','fontsize',14);
         title('All tracks')
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
