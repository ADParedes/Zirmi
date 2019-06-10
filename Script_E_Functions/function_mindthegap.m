% function [mmp_I bwfinal ] = mindthegap(I_TR,bw_color_realsegout,bw,diameter)
% Shortcut summary goes here
%- presets
close all;
clc;
 I           = I_TR;%Tr0;
 canvas      = I_TR;
 i           = 1;
  diameter    = 60; 
 [ix,iy]     = size(I);

 x10=0;y10=0;
%%%%Draw Circles 
if i==1 || i==3  %or statement
%     f=figure(1);imagesc(I); title('ROS REGION');    
%     [image_x_len, image_y_len, bw_woundgap, x10, y10] = roipoly();
    %--
    disp('Select right above notochord line to being CTCF')
    f=figure(1);imagesc(bw_color_realsegout); title('Outline Wound Gap');
    [x, y]          = getline(f);
    bw_woundgap     = poly2mask(x,y,ix,iy);
    exist_woundgap  = sum(sum(bw_woundgap));
    if exist_woundgap>0
        exist_woundgap          = 1;
        bw_woundgap_outline    = bwperim(bw_woundgap);
        figure(2);
        Segout_woundgap        = imoverlay(I,bw_woundgap,[1 0 0]);
        imshow(Segout_woundgap);
    else
    end;

else
end;
%% Calculate the Metrics 
%--find how many points have been established 

num_points=length(x);
d       =(0);
n       =(0);
x1      = 0;
y1      = 0;
%-circle centers
Arr_xy  ={0};
Arr_x   ={0};
Arr_y   ={0};
%-masks
CMASK={0}; %  = combo;
RCMASK={0};%  = raw_combo;
for i=1:(num_points-1)
    xy          = 0; % saves two vertical arrays of centers of circles
    x1          = x(i);
    y1          = y(i);
    % Considering for loop condition this is impossible 
    if i == num_points        
        x2 = x(1);
        y2 = y(1);
    else
        x2 = x(i+1);
        y2 = y(i+1);
    end;
    f=figure(2+i);imagesc(I); 
    hold on
    plot(x,y,'LineWidth',2,'Color','r')
    plot(x1,y1,'c*','LineWidth',2)
    plot(x2,y2,'y*','LineWidth',2)
    hold on
    % calculate metrics
    d(i)        = pdist([x1,y1;x2,y2]);
    radius      = diameter/2;
    n           = d(i);
    num         = round(n); %number of circle points per line
    %--figure out the line equation for the one drawn
    coefficients= polyfit([x1, x2], [y1, y2], 1);
    a           = coefficients (1);
    b           = coefficients (2);
     xmin       = min(x(i:i+1));
     xmax       = max(x(i:i+1));
     ymin       = min(y(i:i+1));
     ymax       = max(y(i:i+1));
     funct      = linspace(xmin,xmax, num); % Adapt n for resolution of graph
    fun_y       = a*funct+b;  
    xy          = [funct',fun_y']; 
     Arr_xy{i}  = xy;
     Arr_x{i}   = funct;
     Arr_y{i}   = fun_y;
 %--Presets for next for loop -------------------------------------%
     lens       = length(funct);    
     Img_cmask  = {0};
     x_j        = 0;
     y_j        = 0;
 %--For Loop for Establishing Circles ------------------------------%
     for j=1:lens
        %---showing purposes
        x_j                 = funct(j);
        y_j                 = fun_y(j);
        [xoc yoc]           = circle([x_j y_j], radius, 1000); % outer circle
        plot(xoc, yoc, '-','linewidth',2,'color','g')%'color',0.5.*rand(1,3));
        %---masking purposes
        cx=x_j  ;   cy=y_j;     r=radius;
        [xm,ym]             = meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
        c_mask              =((xm.^2+ym.^2)<=r^2); 
        cc_mask             = c_mask & bw; %corrected circle mask to be inside fish only
        MaskPixelNum        = sum(sum(cc_mask));        
        MaskPixelIntensity  = sum(sum(uint16(cc_mask).*Tr0))/MaskPixelNum;
        mmpi(j)             = MaskPixelIntensity;
        Img_cmask{j}        = cc_mask;
        if j==1
            raw_combo   = c_mask;
            combo       = cc_mask;
        else
            raw_combo   = raw_combo + c_mask;
            raw_combo   = im2bw(raw_combo);
            combo       = combo + cc_mask;
            combo       = im2bw(combo);
        end;

        Img_combomask{j}    = combo;

        hold on
    end;
    CMASK{i}    = combo;
    RCMASK{i}   = raw_combo;
    
end;
%%  Final Image
%After you found a, You can get b from your equation y=a*x+b 
%Creating Circle Mask BW
%-Raw BW before filtering

raw_blank   = I_TR;
blank       = I_TR;
 for kk=1:i      
     if kk==1
         %--ultimate combo
         ultimate_mask  = CMASK{kk};
         %--ultimate raw
         ultimateraw_mask= RCMASK{kk};
     else
         %--ultimate combo
         ultimate_mask  = ultimate_mask + CMASK{kk};
         ultimate_mask  = im2bw(ultimate_mask);
         %--ultimate raw
         ultimateraw_mask   = ultimateraw_mask + RCMASK{kk};
         ultimateraw_mask   = im2bw(ultimateraw_mask);
     end;
     blank      = imoverlay(blank,bwperim(CMASK{kk}),0.3.*rand(1,3));
     raw_blank  = imoverlay(raw_blank,bwperim(RCMASK{kk}),0.6.*rand(1,3));
 end;
RUSegout    = imoverlay(I_TR,bwperim(ultimateraw_mask),[1 1 1]);
USegout     = imoverlay(I_TR,bwperim(ultimate_mask),[1 1 0]);
 
figure(); imshow(RUSegout); title('Raw Mask'); hold on  
plot(x,y,'LineWidth',2,'Color','r')   
hold off
figure(); imshow(USegout); title('Fish Mask'); hold on  
plot(x,y,'LineWidth',2,'Color','r')   
hold off

% imshow(bw_woundgap)
%% Mind the Wound Gap
disp('This is where we need to mind the gap')
disp('Edit here, its wrong')
    if exist_woundgap>0
        exist_woundgap          = 1;
        bwfinal = ~ultimate_mask | bw_woundgap;
        bwfinal = ~bwfinal;
%         figure();imshow(bwfinal);
    else
    end;    
FinalSegout     = imoverlay(I_TR,bwperim(bwfinal),[1 1 0]);
figure(); imshow(FinalSegout); title('Fish Mask'); hold on  
plot(x,y,'LineWidth',2,'Color','r')   
hold off
% close all;
            MaskPixelNum        = sum(sum(bwfinal));        
            MaskPixelIntensity  = sum(sum(uint16(bwfinal).*I_TR))/MaskPixelNum;
            mmp_I             = MaskPixelIntensity;
