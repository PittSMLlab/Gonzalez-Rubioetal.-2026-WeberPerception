function Boxplot(data,color,plotIndivFlag,plotMean,position)

% Width and scaling factor
if position>100
    w=100;
elseif position<10
    w=20;
end

% Make the boxplot

if plotIndivFlag
    % Individual data points (jittered)
    jitter=(rand(size(data))-0.5)*w/2;
    % jitter=-0.5*w/2;
    scatter(position+jitter+w,data,20,color,'filled','MarkerFaceAlpha',0.4); % ,'MarkerEdgeColor','k'
end

% Boxplot manually
q=quantile(data,[0.25 0.5 0.75]);
IQR=q(3)-q(1);
% whisker_low=q(1)-1.5*IQR;
% whisker_high=q(3)+1.5*IQR;
whisker_low=max(min(data),q(1)-1.5*IQR);
whisker_high=min(max(data),q(3)+1.5*IQR);

% Draw box
rectangle('Position',[position-w/2,q(1),w,IQR],'EdgeColor',color,'LineWidth',1.5);

% Draw median line
line([position-w/2,position+w/2],[q(2),q(2)],'Color',color,'LineWidth',2);

% Whiskers
line([position,position],[whisker_low,q(1)],'Color',color,'LineWidth',1);
line([position, position],[q(3),whisker_high],'Color',color,'LineWidth',1);

% Mean dot
if plotMean
    mean_val=mean(data,'omitnan');
    scatter(position,mean_val,40,'d','MarkerEdgeColor','k','MarkerFaceColor',color); % change from d to 0 if you want to change the marker shape
end

end