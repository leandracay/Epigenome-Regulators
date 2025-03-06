    %Criteria for activation can be determined by 1)must reach a fold change that is some percent max of fold all CRs (FC threshold)  
  %2) must reach a PR that is some percent of the max PR for each CR or 3) all
  %the CRs
function [timeDelayOn,timeDelayOff]=TimeDelay(mCherry,PR,timeVector)
fig1=figure;
for i=1:size(mCherry,1)
  if nansum(isnan(PR(i,1:14)))>=3
      timeDelayOn(i)=NaN;
  else
[PRsorted,I]=sort(PR(i,1:40),'MissingPlacement','first');
maxPR(i)=nanmedian(PRsorted(end-2:end));
ind=find(PR(i,:)>=0.5.*maxPR(i));
if isnan(maxPR(i))==1
    timeDelayOn(i)=NaN;
elseif maxPR(i)<=0
    timeDelayOn(i)=NaN;
else
timeDelayOn(i)=timeVector(min(ind));
end
  end
  figure(fig1)
  %plot(timeVector-timeDelayOn(i),PR(i,:))
  hold on
  patch(timeVector-timeDelayOn(i),i.*ones(length(PR(i,:)),1),PR(i,:),PR(i,:),'FaceColor','none','EdgeColor','interp','Linewidth',2)
c=colorbar;
view(3)
end
figure
histogram(timeDelayOn,-100:20:800)
timeDelayOff=NaN;
xlabel('Time delay on (min)')


end
