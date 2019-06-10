function [  ] = function_zirmi_pca( data,str_test,arr,num,str,input_als)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




%% Find the Prinicipal Compoents of the ingredents data
    % latent = eignevalues of covariance matrix
    % filtered.Status         = temp.Survival(idx_toi);  %Group S vs NS
    % arr.Group               = temp.Group(idx_toi); %Group Radiation Group
switch input_als
    case ('n')
        [coeff,score,~ ,~,explained, ~]     = pca(data);
    otherwise
         [coeff,score,~ ,~,explained, ~]    = pca(data,'algorithm','als');
end;
%---Sort the Variables and Groups
idx_baseline                            = find(arr.Group ==num.Groups{1});
idx_0J                                  = find(arr.Group ==num.Groups{2});
idx_3J                                  = find(arr.Group ==num.Groups{3});
idx_9J                                  = find(arr.Group ==num.Groups{4});
idx_18J                                 = find(arr.Group ==num.Groups{5});

disp('End: Derive Coefficent or eigen values or Influence')
%% Slide 1

slideId     = exportToPPTX('addslide');

%--Figure 1 SCREE PLOT 
figure('Renderer','zbuffer'); pareto(explained);
    set(gca,'Fontsize',14)
    tx      = xlabel('Principal Component','FontSize',18,...
                    'FontWeight','bold','Color','k');
    ty      = ylabel('Variance Explained (%)','FontSize',18,...
                    'FontWeight','bold','Color','k');%'Position',[2.5 50 -1],...
% Upper Left Corner
exportToPPTX('addpicture',gcf,'Position',[1 0.5 3 2],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','SCREE PLOT','Position',[1 0 3 0.5],'Vert','bottom');
close(gcf);
disp('End: Scree Plot')
%--Figure 2 Dose(J/cm^2) Score Plot
A_score=score((idx_baseline),:);
B_score=score((idx_0J),:);
C_score=score((idx_3J),:);
D_score=score((idx_9J),:);
E_score=score((idx_18J),:);
figure('Name','Score Plot  2-D','Renderer','zbuffer')
    plot(A_score(:,1),A_score(:,2),'r+','LineWidth',4); hold on
    plot(B_score(:,1),B_score(:,2),'b+','LineWidth',4); hold on
    plot(C_score(:,1),C_score(:,2),'g+','LineWidth',4); hold on
    plot(D_score(:,1),D_score(:,2),'k+','LineWidth',4); hold on
    plot(E_score(:,1),E_score(:,2),'m+','LineWidth',4); hold on
    legend(str.Groups(:),'Location','best')
    grid('on')
%     xlim([-100,100]) %This fixes axis'
%     ylim([-100,100]) %This fixes axis'
    set(gca,'Fontsize',14);
    xlabel('1st Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
    ylabel('2nd Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
    ax = gca;
    ax.LineWidth = 3;
% Bottom Left Corner
exportToPPTX('addpicture',gcf,'Position',[1 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Dose(J/cm^2) Score PLOT','Position',[1 3 3 0.5],'Vert','bottom');
close(gcf);
disp('End: Radiation Score Plot');   
%--Figure 3 Surival Status Score Plot 
% S_score=score((idx_S),:);
% NS_score=score((idx_NS),:);
figure('Name','Survival Status Score Plot','Renderer','zbuffer')
%     plot(S_score(:,1),S_score(:,2),'r+','LineWidth',4); hold on
%     plot(NS_score(:,1),NS_score(:,2),'b+','LineWidth',4); hold on
%     legend(str.Status(:),'Location','best')
    grid('on')
%     xlim([-100,100])
%     ylim([-100,100])
    set(gca,'Fontsize',14);
    xlabel('1st Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
    ylabel('2nd Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
    ax = gca;
    ax.LineWidth = 3;
% Bottom Right Corner
exportToPPTX('addpicture',gcf,'Position',[6 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Surivival Status Score PLOT','Position',[6 3 3 0.5],'Vert','bottom');
close(gcf);
disp('Score End: Plot')
%--Figure 4 BiPlot
figure('Name','Bi Plot 2-D','Renderer','zbuffer')
    biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',str.var1);
    set(gca,'Fontsize',14)
    xlabel('1st Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
    ylabel('2nd Principal Component','FontSize',18,...
           'FontWeight','bold','Color','k')
%     mins    = min(min(coeff(:,1:2)));
%     maxs    = max(max(coeff(:,1:2)));
%     axis([mins maxs mins maxs]);
    view(2)
% Top Right Corner
exportToPPTX('addpicture',gcf,'Position',[6 0.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Bi PLOT 2-D','Position',[6 0 3 0.5],'Vert','bottom');
close(gcf);
disp('Score End: Plot')
% Add note to the slide
exportToPPTX('addnote',str_test,'FontWeight','bold');
% Add text to cell
exportToPPTX('addtext',str_test, ...
    'Position',[8 0 4 2], ...  %currently right [4 0 4 2] Unused is center
    'HorizontalAlignment','center', ...
    'Color','k');
%% Slide 2
slideId     = exportToPPTX('addslide');
%--- 3 Dimensional Score plot
figure('Renderer','zbuffer')
    biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',str.var1);
%     mins    = min(min(coeff(:,1:3)));
%     maxs    = max(max(coeff(:,1:3)));
%     axis([mins maxs mins maxs]);
    view([30 40]);
exportToPPTX('addpicture',gcf,'Position',[6 0.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Bi PLOT 3-D','Position',[6 0 3 0.5],'Vert','bottom');
close(gcf);

figure('Name','Score Plot  3-D','Renderer','zbuffer')
    scatter3(A_score(:,1),A_score(:,2),A_score(:,3),'r+','LineWidth',10); hold on
    scatter3(B_score(:,1),B_score(:,2),B_score(:,3),'b+','LineWidth',10); hold on
    scatter3(C_score(:,1),C_score(:,2),C_score(:,3),'g+','LineWidth',10); hold on
    scatter3(D_score(:,1),D_score(:,2),D_score(:,3),'k+','LineWidth',10); hold on
    scatter3(E_score(:,1),E_score(:,2),E_score(:,3),'m+','LineWidth',10); hold on
    axis equal
    legend(str.Groups(:),'Location','best')
    % xlim([-100,100])
    % ylim([-100,100])
    % zlim([-100,100])
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    % Upper Left Corner
exportToPPTX('addpicture',gcf,'Position',[1 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Radiation 3D Score PLOT','Position',[1 3 3 0.5],'Vert','bottom');
close(gcf);

%---
figure('Name','Score Plot  3-D','Renderer','zbuffer')
%     scatter3(S_score(:,1),S_score(:,2),S_score(:,3),'r+','LineWidth',10); hold on
%     scatter3(NS_score(:,1),NS_score(:,2),NS_score(:,3),'b+','LineWidth',10); hold on
    axis equal
%     legend(str.Status(:),'Location','best')
    % xlim([-100,100])
    % ylim([-100,100])
    % zlim([-100,100])
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    % Upper Left Corner
exportToPPTX('addpicture',gcf,'Position',[6 3.5 4 3],'EdgeColor',[0 0 0.8],'LineWidth',3);
exportToPPTX('addtext','Surivival 3D Score PLOT','Position',[6 3 3 0.5],'Vert','bottom');
close(gcf);

exportToPPTX('addnote',str_test,'FontWeight','bold');
% Add text to cell
exportToPPTX('addtext',str_test, ...
    'Position',[8 0 4 2], ...  %currently right [4 0 4 2] Unused is center
    'HorizontalAlignment','center', ...
    'Color','k');
end

