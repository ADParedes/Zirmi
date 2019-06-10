function [ TeRe1 ] = zirmiAdelim( numfids,I,i,str_part)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
TeRe1 = cell(numfids(I),1);
for bi = 1:numfids(I)
   ti=bi-1;
   if numfids(I)>99
       %% Over 100 time Points
       if bi==1
           nameTR1  = strcat(str_part,'*t000','*p',num2str(i),'*');
       elseif ti<10
            nameTR1  = strcat(str_part,'*t00',num2str(ti),'*p',num2str(i),'*');
       elseif ti<100
            nameTR1  = strcat(str_part,'*t0',num2str(ti),'*P',num2str(i),'*');
       else
            nameTR1  = strcat(str_part,'*t',num2str(ti),'*P',num2str(i),'*');
       end;
   elseif numfids(I)<10
       %% Only 10 time Points
       nameTR1  = strcat(str_part,'*t',num2str(ti),'*p',num2str(i),'*');
   elseif numfids(I)==10
       %% Less than 10 Time Points 
       if bi==1
           nameTR1  = strcat(str_part,'*t00','*p',num2str(i),'*');
       elseif ti<10
            nameTR1  = strcat(str_part,'*t0',num2str(ti),'*p',num2str(i),'*');
       else
            nameTR1  = strcat(str_part,'*t10',num2str(ti),'*P',num2str(i),'*');
       end;                  
   else
       %% More than 10 Less than 100 Time Points 
       if ti==0
           nameTR1  = strcat(str_part,'*t00','*p',num2str(i),'*');
       elseif ti<10
            nameTR1  = strcat(str_part,'*t0',num2str(ti),'*p',num2str(i),'*');
       elseif ti<100
            nameTR1  = strcat(str_part,'*t',num2str(ti),'*p',num2str(i),'*');
       else
%                             nameTR1  = strcat(str_part,'*t0',num2str(ti),'*P',num2str(i),'*');
       end;
   end;
    gg1       = dir(nameTR1);
    EE1       = [gg1(:).datenum].'; % you may want to eliminate . and .. first.
    [EE1,EE1] = sort(EE1);
    TeRe1{bi} = {gg1(EE1).name}; % Cell array of names in order by datenum.  TeRe{Position}{TimeInterval               
end;


