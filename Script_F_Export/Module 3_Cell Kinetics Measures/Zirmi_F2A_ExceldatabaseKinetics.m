%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%ZEBRAFISH EXCEL DOCTOR IT
%Move Excel  GTs 60, 90, 120, 150, 180 and max to excel.
%Updated 2017/11/24
clearvars -except Worthy SM Frame POI PARAMETERS ADP PhagoSight handles dataIn dataL dataR ch_GFP ch_Ph2
disp('START: Zirmi F: Create Excel Tables')
[c1] = clock;
c1_fix = fix(c1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recognized GUI-ed Parameters
clear Str Num Dirs Cell Array % New Structure Array
Str.name                = POI.Parameter10d;
Str.expnum              = POI.Parameter10a(3:end);
Str.MetricType          ={'Velocity(um/min)','StaticRatio', 'MeanderingRatio',...  Velocity is still pixel per frame
                            'Tortuosity(abu)','ForwardRatio','ForwardtoBackwardRatio'...
                            'WoundDistance(um)','WoundScore(1-4)'};
Str.MetricName          ={'Velocity','StaticRatio', 'MeanderingRatio',...
                                'Tortuosity','ForwardRatio','ForwardtoBackwardRatio'...
                                'WoundDistance','WoundScore','ID'};
Str.MetricSM            ={'velocity','staticratio', 'meandering',...
                                'tortuosity','forward','FBratio'...
                                'WoundScoreUm','WoundScore1234'};
Str.Sheets              ={'SumAll','MpAll','Score1','Score2',...
                            'Score3','Score4'};
Str.DataType            ={'interval','intervalsm','cum',...
                            'phagosight'};
Str.DataTypeFile        ={'Interval','IntervalSM','Cumulative'};
Str.GT                  ={'GT_60','GT_90','GT_120','max'}; 
Str.Sheet1              ={'min/frame','MPI','Unique ExpID',...
                            'Velocity(um/min)','StaticRatio', 'MeanderingRatio',...
                            'Tortuosity(abu)','ForwardRatio','ForwardtoBackwardRatio'...
                            'WoundDistance(um)','WoundScore(1-4)',...   %11 Points
                            'STD1','STD2','STD3','STD4','STD5','STD6','STD7','STD8'}; %18 Points
Str.Sheet23456          ={'Unique ExpID','Unique MpID',...
                            'Velocity(um/min)','StaticRatio', 'MeanderingRatio',...
                            'Tortuosity(abu)','ForwardRatio','ForwardtoBackwardRatio'...
                            'WoundDistance(um)','WoundScore(1-4)'}; %10 Points    
                                                         
Num.exp                 = str2num(Str.expnum); %#ok<ST2NM>
Num.mp                  = length(handles.finalLabel); %never changes this is the numofTracks
Num.micronsPerPixel     = PARAMETERS.Parameter2;
Num.SamplingFrequency   = PARAMETERS.Parameter4;
Num.MPI                 = PARAMETERS.Parameter5;
Num.TimeIntervals       = 3;  %only dooing GT 60 90 and 120
Num.Scores              = 4;
Num.Sheets              = 6;
Num.DataType            = length(Str.DataType);

Dirs.Main               = ADP.dir_metadat;

Arrays.MetricSM         ={SM.velocity,SM.staticratio,SM.meandering,...
                            SM.tortuosity,SM.forward,SM.FBratio,...
                            SM.WoundScoreUm,SM.WoundScore1234,Worthy.ID};
display('End: Recognize Parameters')
%% Directory Branching/Sorting For Saving .xls files
switch ADP.boo2 
    case {1,0}
        disp('Hello ADP')
    %-Experiment Type Descriptor
%         Str.ExpSet             ='test';
         Str.ExpSet           = input('User-defined Identifier e.g. Baseline,0J,3J,9J,DoubleWound,Manuscript1 : ','s');
        if Num.exp<46
            Str.Fish           = 'mpegdendra2';
        elseif Num.exp>46
            Str.Fish           = 'mpxdendra2';
        else
            disp('error in Directory Sorting')
        end;

    case {2}
        disp('Hello User')
        Str.ExpSet           = input(strcat('Please Input Experiment Set Descriptor',...
                                    '(E.g. Set#1 or 9J'),'s');
        Str.Fish           = input(stcat('Please Input Tg(Fish) Descriptor',...
                                   '(E.g. mpeg or mpx dendra2 '),'s'); 
    otherwise 
        warning('Zirmi Cannot Recognize Your Computer - will error')   
end;
display('End: Directory Labeling')
%% THE FOR LOOP  
for J=1:Num.TimeIntervals
    %% Main Directory
    clear FishMean RawMetrics  %Structure Array for the entire fish mean
    clear Interval IntervalSM Cumulative
    cd(Dirs.Main)
    %% Sort into Time Specific Structure Arrays
    for ind             = 1:length(Str.MetricName);
        switch Str.MetricName{ind}
            case {'ID'}
                Interval.(Str.MetricName{ind})              = Arrays.MetricSM{ind}.interval{J};
                IntervalSM.(Str.MetricName{ind})            = Arrays.MetricSM{ind}.intervalsm;   
                Cumulative.(Str.MetricName{ind})            = Arrays.MetricSM{ind}.cum;
            otherwise  
                Interval.(Str.MetricName{ind}).all          = Arrays.MetricSM{ind}.interval{J};
                IntervalSM.(Str.MetricName{ind}).all        = Arrays.MetricSM{ind}.intervalsm{J};   
                Cumulative.(Str.MetricName{ind}).all        = Arrays.MetricSM{ind}.cum{J};
                
                Interval.(Str.MetricName{ind}).mean          = nanmean(Arrays.MetricSM{ind}.interval{J});
                IntervalSM.(Str.MetricName{ind}).mean        = nanmean(Arrays.MetricSM{ind}.intervalsm{J});
                Cumulative.(Str.MetricName{ind}).mean        = nanmean(Arrays.MetricSM{ind}.cum{J});
                
                Interval.(Str.MetricName{ind}).std          = nanstd(Arrays.MetricSM{ind}.interval{J});
                IntervalSM.(Str.MetricName{ind}).std        = nanstd(Arrays.MetricSM{ind}.intervalsm{J});
                Cumulative.(Str.MetricName{ind}).std        = nanstd(Arrays.MetricSM{ind}.cum{J});                
        end;
    end;

%-Metric#7: Forward to Backward (normalized)  
%-Metric#8: WoundScoreUm  - is the distance to the wound gap from initial
%point to the final point | negative(-) means closer ; positive(+) means moved farther away.
%(ADP Labnotebook pg 92)
    disp('End: Sorted into Time Structure Arrays');
    %% Sort into Wound Score Structure Arrays
    i   = 0;
    for I = 1:Num.Sheets; 
        temp1           = Str.Sheets{I};
        switch temp1
            case{'Score1','Score2','Score3','Score4'}
                i                       = i+1;
                Interval.(temp1)        = find(Interval.WoundScore.all == i);
                IntervalSM.(temp1)      = find(IntervalSM.WoundScore.all == i);
                Cumulative.(temp1)      = find(Cumulative.WoundScore.all == i);
                for ind                 = 1:length(Str.MetricName);                   
                   switch Str.MetricName{ind}
                       case {'ID'}
                            Interval.xls.(Str.Sheets{I}).vals(:,1)          = Interval.(Str.MetricName{ind})(Interval.(temp1))';
                            IntervalSM.xls.(Str.Sheets{I}).vals(:,1)        = IntervalSM.(Str.MetricName{ind})(IntervalSM.(temp1))';   
                            Cumulative.xls.(Str.Sheets{I}).vals(:,1)        = Cumulative.(Str.MetricName{ind})(Cumulative.(temp1))';  
                       otherwise                          
                            Interval.xls.(Str.Sheets{I}).vals(:,ind+1)          = Interval.(Str.MetricName{ind}).all(Interval.(temp1))';
                            IntervalSM.xls.(Str.Sheets{I}).vals(:,ind+1)        = IntervalSM.(Str.MetricName{ind}).all(IntervalSM.(temp1))';   
                            Cumulative.xls.(Str.Sheets{I}).vals(:,ind+1)        = Cumulative.(Str.MetricName{ind}).all(Cumulative.(temp1))'; 
                   end;
                end;
            case{'MpAll'}
                for ind                 = 1:length(Str.MetricName);                   
                   switch Str.MetricName{ind}
                       case {'ID'}
                            Interval.xls.(Str.Sheets{I}).vals(:,1)          = Interval.(Str.MetricName{ind})';
                            IntervalSM.xls.(Str.Sheets{I}).vals(:,1)        = IntervalSM.(Str.MetricName{ind})';   
                            Cumulative.xls.(Str.Sheets{I}).vals(:,1)        = Cumulative.(Str.MetricName{ind})';  
                       otherwise                          
                            Interval.xls.(Str.Sheets{I}).vals(:,ind+1)          = Interval.(Str.MetricName{ind}).all';
                            IntervalSM.xls.(Str.Sheets{I}).vals(:,ind+1)        = IntervalSM.(Str.MetricName{ind}).all';   
                            Cumulative.xls.(Str.Sheets{I}).vals(:,ind+1)        = Cumulative.(Str.MetricName{ind}).all'; 
                   end;
                end;
            otherwise
                for ind2                 = 1:length(Str.MetricName);
                    switch Str.MetricName{ind2}
                        case {'ID'}
%                             Interval.xls.(Str.Sheets{I}).vals(:,ind2)         = Interval.(Str.MetricName{ind2})';
%                             IntervalSM.xls.(Str.Sheets{I}).vals(:,ind2)        = IntervalSM.(Str.MetricName{ind2})';   
%                             Cumulative.xls.(Str.Sheets{I}).vals(:,ind2)        = Cumulative.(Str.MetricName{ind2})';  
                        otherwise
                            Interval.xls.(Str.Sheets{I}).vals(:,ind2)          = Interval.(Str.MetricName{ind2}).mean';
                            IntervalSM.xls.(Str.Sheets{I}).vals(:,ind2)        = IntervalSM.(Str.MetricName{ind2}).mean';   
                            Cumulative.xls.(Str.Sheets{I}).vals(:,ind2)        = Cumulative.(Str.MetricName{ind2}).mean';     
                            
                            Interval.xls.(Str.Sheets{I}).std(:,ind2)          = Interval.(Str.MetricName{ind2}).std';
                            IntervalSM.xls.(Str.Sheets{I}).std(:,ind2)        = IntervalSM.(Str.MetricName{ind2}).std';   
                            Cumulative.xls.(Str.Sheets{I}).std(:,ind2)        = Cumulative.(Str.MetricName{ind2}).std'; 
                    end;
                end
        end;
    end;
      disp('End: Sorted into Wound Score Structure Arrays');
    %% Sheet 1: Mean
    clear Zirmi Sheets
    cd(Dirs.Main)
    mkdir(Dirs.Main,Str.Fish);cd(Str.Fish)
    Dirs.Fish          = cd;
    mkdir(Dirs.Fish,Str.ExpSet);cd(Str.ExpSet)
    Dirs.ExpSet          = cd;
    Zirmi.Interval     = Interval;
    Zirmi.IntervalSM   = IntervalSM;
    Zirmi.Cumulative   = Cumulative;

    %Note: This is For GT_'J'
    for ind4           = 1:length(Str.DataTypeFile);
        cd(Dirs.ExpSet) % Relocate To ExperimentSet (E.g. 0J,3J,18J...)
        
        for ind5            = 1:length(Str.Sheets);
            temp1               = Str.Sheets{ind5};
            temp2               = Zirmi.(Str.DataTypeFile{ind4}).xls.(Str.Sheets{ind5}).vals;
            [row,col]           = size(temp2);      
            switch temp1
                case{'MpAll','Score1','Score2','Score3','Score4'}
                            placeholder         = num2cell(temp2);
                            tempA               = cell(row,1);
                            tempA(:)            = {Str.name};
                            Sheets.num{ind5}        =placeholder;
                            Sheets.raw{ind5}        =[tempA,placeholder];
                            Sheets.title{ind5}  = [Str.Sheet23456;Sheets.raw{ind5}];
                    
                otherwise 
                            temp3               = Zirmi.(Str.DataTypeFile{ind4}).xls.(Str.Sheets{ind5}).std;
                            placeholder2        = num2cell(temp2);
                            placeholder3        = num2cell(temp3);
                            Sheets.num{ind5}        =[placeholder2,placeholder3];
                            Sheets.raw{ind5}        =[[Num.SamplingFrequency],...
                                                        Num.MPI,...
                                                        Str.name,[Sheets.num{ind5}]];
                            Sheets.title{ind5}  = [Str.Sheet1;Sheets.raw{ind5}];
            end;
            file_excel         = strcat(Str.DataTypeFile{ind4},...
                                            '_',Str.GT{J},'.xlsx');%DataType_Time
            exist_excel(ind5)        = exist(file_excel,'file');
            
            switch exist_excel(1) 
                case{0} %file does not exist
%                     xlswrite(file_excel,Sheets.title{ind5},ind5)
                otherwise %file does exist
                    [NUM_1,TXT_1,RAW_1] = xlsread(file_excel,ind5);
                    empty_raw           = isempty(RAW_1);
                    switch empty_raw
                        case{0} %is not empty
                            switch temp1 
                                case{'MpAll','Score1','Score2','Score3','Score4'}
                                    Sheets.combo{ind5}=[RAW_1;Sheets.raw{ind5}];
                                    xlswrite(file_excel,Sheets.combo{ind5},ind5)
                                otherwise
                                    Sheets.combo{ind5}=[RAW_1;Sheets.raw{ind5}];
                                    xlswrite(file_excel,Sheets.combo{ind5},ind5)
                            end;
                                                        
                        otherwise %is empty
%                             xlswrite(file_excel,Sheets.title{ind5},ind5) 
                    end;
            end;
                    
           
        end;
        disp(strcat('CheckPoint DataType-',Str.DataTypeFile{J}))
    end;
      

    disp(strcat('Excel End_',Str.GT{J}))
end;
disp(Str.name)
[c2]    = clock;
c2_fix  = fix(c2);
str_clock = strcat('Start |',num2str(c1_fix(4)),':',num2str(c1_fix(5)),':',num2str(c1_fix(6))',...
                   'END: |',num2str(c2_fix(4)),':',num2str(c2_fix(5)),':',num2str(c2_fix(6)));
disp(str_clock)               
disp('FINISHED: Zirmi_F Transfer to .xls Files')