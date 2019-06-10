%% Plots ROS
%% Load Prevoius wound_CTCF
 cd('C:\Users\AndreDaniel\OneDrive\PhD Data\m_ROS\CTCF');
 str_inputCTCFname=input('What was the Experiment and Position & Wound/noWound E.g. AB001P0Wound?','s');
 str_CTCFname=strcat(str_inputCTCFname,'.mat');
 load(str_CTCFname);
 

 %%
 wound_IntDen  = (PixelNum).* mfi; %Integrated Density
wound_CTCF            = wound_IntDen - ((PixelNum).*MeanFluorescence);%Corrected Total Cell Fluorescence 
%%
pause();
namename        =strcat(name,name5);
disp('Load the MeanFlourescence of Wound')
cd('C:\Users\AndreDaniel\OneDrive\PhD Data\m_ROS\CTCF');
woundname    =strcat(namename,'Wound');

save(woundname,'mmpi0','mfi','PixelNum','wound_IntDen','wound_CTCF','Time','frames','arr_T','T_p','name5','name','t_plate','ti_d','golditeration','Img1');
disp('END: Saved Wound CTCF');