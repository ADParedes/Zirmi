function [ output ] = zirmiChVal( input )
% input is a string
% Assumption is that you have more than 1
% Assumption that you have no more than 1000
%   Returns the Limits of parts, positions, and Z stack automatically
    check                       =   0;
    output                      =   0;
    %--Positions-------------------------------------------------
    %Stores Position Numbers into variable
    %If variable Destination ev has values for Position # it increase
    %Once Position # is high enough it stops and knows Max Position Length
        while check == 0  %%Note we are assuming there are Z/if not will LOOP CONTINOUSLY
                output      =    output+1;
                ev          =    strcat(input,num2str(output),'*'); %Problem: Preinjury and Parts have different Posnum
                f           =    dir(ev);
                check       =    isempty(f);
                 if output>1000
                     output    =   1;
                     break
                 end;
        end;
     %posnum   =  posnum-1;%-1; Don't minus 1 because p0 is a position.

end

