function h=bargraph(series, error)

h = bar(series);
set(h,'BarWidth',1); % The bars will now touch each other

% set(gca,'YGrid','on')
% set(gca,'GridLineStyle','-')
% set(gca,'XTicklabel','Prestim| Stim period')
% set(get(gca,'YLabel'),'String','Mean FR (spikes/sec)','Fontsize',12)
% set(h(1),'FaceColor','k')
% set(h(2),'FaceColor','b')
% if sum(series(2,:))         % this is a hack; need better solution
%     set(h,'BarWidth',1); % The bars will now touch each other
%     set(gca,'XTicklabel','Prestim| Stimulus period| Blank trials')
%     numgroups = size(series,1);
% else
%     numgroups = 1;
%     series = series(1,:);
% end

lh = legend('No light','Light');
% set(lh,'Location','BestOutside','Orientation','horizontal')

hold on;

numgroups = size(series, 1);
numbars = size(series, 2);

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, series(:,i), error(:,i), 'k', 'linestyle', 'none');
end
hold off
return