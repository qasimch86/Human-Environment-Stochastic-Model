function SocLearn_Figure(X1, Y1, Y2, Y3, Y4)
%CREATEFIGURE(X1,Y1,Y2,Y3,Y4)
%  X1:  vector of x data
%  Y1:  vector of y data
%  Y2:  vector of y data
%  Y3:  vector of y data
%  Y4:  vector of y data

%  Auto-generated by MATLAB on 12-May-2014 21:29:41

% Create figure
figure1 = figure('PaperType','<custom>','PaperSize',[53.25 75.33],...
    'PaperUnits','centimeters');

% Create axes
axes('Parent',figure1,'ZColor',[1 1 1],'YTickLabel','','YTick',zeros(1,0),...
    'YColor',[1 1 1],...
    'XTickLabel','',...
    'XTick',zeros(1,0),...
    'XColor',[1 1 1],...
    'Position',[0.05405 0.05368 0.9346 0.8975],...
    'CLim',[0 1]);

% Create title
title('Variation in time-to-crosspatch-infestation due to social norms',...
    'FontSize',16);

% Create xlabel
xlabel('Social norm "n"','FontSize',16,'Color',[0 0 0]);

% Create ylabel
ylabel('Time-to-cross-patch-infestation t_c_r_o_s_s (in years)','FontSize',16,'Color',[0 0 0]);

% Create axes
axes1 = axes('Parent',figure1,...%'YTickLabel',{'0','50','100','150','200','250','300'},...%     'YTick',[0 50 100 150 200 250 300],...
    'Position',[0.1003 0.5549 0.4026 0.3536],...
    'FontSize',12);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 300]);
grid(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Parent',axes1,'LineWidth',2,'DisplayName','Patch 1',...
    'Color',[0 0 0]);

% Create title
title('Patch 1','FontSize',14);

% Create axes
axes2 = axes('Parent',figure1,'Position',[0.5811 0.5586 0.3968 0.3516],...
    'FontSize',12);
grid(axes2,'on');
hold(axes2,'all');

% Create plot
plot(X1,Y2,'Parent',axes2,'LineWidth',2,'DisplayName','Patch 1',...
    'Color',[0 0 0]);

% Create title
title('Patch 2','FontSize',14);

% Create axes
axes3 = axes('Parent',figure1,'Position',[0.09882 0.1044 0.407 0.3571],...
    'FontSize',12);
grid(axes3,'on');
hold(axes3,'all');

% Create plot
plot(X1,Y3,'Parent',axes3,'LineWidth',2,'DisplayName','Patch 1',...
    'Color',[0 0 0]);

% Create title
title('Patch 3','FontSize',14);

% Create axes
axes4 = axes('Parent',figure1,'Position',[0.5811 0.1026 0.3968 0.3608],...
    'FontSize',12);
grid(axes4,'on');
hold(axes4,'all');

% Create plot
plot(X1,Y4,'Parent',axes4,'LineWidth',2,'DisplayName','Patch 1',...
    'Color',[0 0 0]);

% Create title
title('Patch 4','FontSize',14);
