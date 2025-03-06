 %Criteria for activation can be determined by 1)must reach a fold change that is some percent max of fold all CRs (FC threshold)  
  %2) must reach a PR that is some percent of the max PR for each CR or 3) all
  %the CRs
function [percent_activated,timeDelayOn,timeDelayOff]=TimeDelay_FC(mCherry,PR,timeVector,threshold)

%fig1=figure;
timeDelayOn=NaN(1,size(mCherry,1));
for i=1:size(mCherry,1)
  if nansum(isnan(PR(i,1:14)))>=6
      timeDelayOn(i)=NaN;
  else

ind=find(mCherry(i,:)>=threshold);
if isempty(ind)==1
    timeDelayOn(i)=inf;
else
timeDelayOn(i)=timeVector(min(ind));
end
  end
%   figure(fig1)
%   %plot(timeVector-timeDelayOn(i),PR(i,:))
%   hold on
%   patch(timeVector-timeDelayOn(i),i.*ones(length(PR(i,:)),1),PR(i,:),PR(i,:),'FaceColor','none','EdgeColor','interp','Linewidth',2)
% c=colorbar;
% view(3)
end
numberOfCells=sum(~isnan(timeDelayOn));
percent_activated=sum(timeDelayOn>0 & timeDelayOn~=inf)/numberOfCells;
figure
subplot(2,1,1)
histogram(timeDelayOn,-100:20:800)
timeDelayOff=NaN;
xlabel('Time delay on (min)')
title(strcat('N=', string(numberOfCells),' Percent activated', string(percent_activated*100),'%'))
subplot(2,1,2)
plot(timeVector,sum(~isnan(mCherry)))

hold on
plot(timeVector,sum(mCherry>=threshold))
xlim([-150,850])
legend('Number of cells segmented','Number of cells activated','Location','northoutside')
xlabel('Time (min)')
ylabel('Number of cells')
end