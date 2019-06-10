%% Select folder with handles
input_16bit     = ADP.boo1;
name            = POI.Parameter10a;
name5           = POI.Parameter11c;
str_Presumptious= POI.Parameter10c;
%%
cd(str_Presumptious)
disp(name5)
disp(name)
SelectedFolder_str=uigetdir
[pa5, name5, ext5] = fileparts(SelectedFolder_str); %break up file name
cd(SelectedFolder_str)
disp(name5)
disp(name)
% SelectedFile_str=uigetfile

% if Folder_newHandles ~= 0 %Unnecessary if handles = handlesNew
%       [oo,pp]=fileparts(Ha_new)
%        plotNeutrophilMovie(oo,1,find(handles.distanceNetwork.numHops>10));  
% else

%       [oo,pp]=fileparts(Ha_new)
      [F1]= plotNeutrophilMovie('handles.mat',3);
%               plotNeutrophilMovie(oo,6,find(handles.distanceNetwork.numHops>2)); %labels--so sick
    % F = plotNeutrophilMovie(Ha(1:10));
% end;
% display('pause for another video?')
% pause()
% [F4]=plotNeutrophilMovie(cd,4);