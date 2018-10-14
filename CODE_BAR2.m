function CODE_BAR2(y1, y2, y3, y4)
%CREATEFIGURE(Y1,Y2,Y3,Y4)
%  Y1:  bar y
%  Y2:  bar y
%  Y3:  bar y
%  Y4:  bar y
 
%  Auto-generated by MATLAB on 11-Aug-2014 21:44:41
 
%% Create figure
figure1 = figure('FileName','C:\Users\Ali\Desktop\Chris Bauch\Programs\stochastic\Working6\Patches arrangements\Compeletly connected\stat_cross.fig');
 
%% Create axes
axes1 = axes(...
  'FontSize',12,...
  'Position',[0.1031 0.06725 0.8822 0.8631],...
  'XTick',[1 2 3 4],...
  'YGrid','on',...
  'Parent',figure1);
ylim(axes1,[0 100]);
title(axes1,'Statistics of cross patch infestation');
xlabel(axes1,'Number of Patches');
ylabel(axes1,'Cross Patch infestation');
hold(axes1,'all');
 
%% Create bar
bar1 = bar(y1,...
  'Parent',axes1,...
  'DisplayName','Number of Infestation',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.8314 0.8157 0.7843],...
  'ShowBaseLine','off');
 
%% Create bar
bar2 = bar(y2,...
  'Parent',axes1,...
  'DisplayName','Average Time',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.502 0.502 0.502],...
  'ShowBaseLine','off');
 
%% Create bar
bar3 = bar(y3,...
  'Parent',axes1,...
  'DisplayName','Variance',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.2471 0.2471 0.2471],...
  'ShowBaseLine','off');
 
%% Create bar
bar4 = bar(y4,...
  'Parent',axes1,...
  'DisplayName','Standard Deviation',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0 0 0],...
  'ShowBaseLine','off');
 
%% Create legend
legend1 = legend(...
  axes1,{'Number of Infestation','Average Time','Variance','Standard Deviation'},...
  'FontSize',12,...
  'Position',[0.1101 0.7446 0.3339 0.1727]);
 
%% Create axes
axes2 = axes(...
  'FontSize',12,...
  'Position',[0.5106 0.3621 0.4534 0.535],...
  'XTick',[1 2 3 4],...
  'YGrid','on',...
  'Parent',figure1);
ylim(axes2,[0 10]);
hold(axes2,'all');
 
%% Create bar
bar5 = bar(y1,...
  'Parent',axes2,...
  'DisplayName','Number of Infestation',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.8314 0.8157 0.7843],...
  'ShowBaseLine','off');
 
%% Create bar
bar6 = bar(y2,...
  'Parent',axes2,...
  'DisplayName','Average Time',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.502 0.502 0.502],...
  'ShowBaseLine','off');
 
%% Create bar
bar7 = bar(y3,...
  'Parent',axes2,...
  'DisplayName','Variance',...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0.2471 0.2471 0.2471],...
  'ShowBaseLine','off');
 
%% Create bar
bar8 = bar(y4,...
  'Parent',axes2,...
  'BarWidth',0.2,...
  'LineWidth',2,...
  'EdgeColor','none',...
  'FaceColor',[0 0 0],...
  'ShowBaseLine','off');
 
