   %% Separate into Arrays
    cds             ={dir_BF,dir_GFP,dir_TexRed};
    prev_folder     = fileparts(dir_exportedBatchBranch);
    cd              (prev_folder);               
    cd              (dir_exportedBatchBranch);
    Arr_Ph2         ={};
    Arr_GFP         ={};  %{Part}{Position}{All in Time001}
    Arr_TexRed      ={};
    for I=1%:parts
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
            %% Determining Fish Position Maximum - Cannot exceed 99 fish per analysis
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
                display('No TexRed Z or only one image file - not coded for this')
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