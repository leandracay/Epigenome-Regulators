 %Criteria for activation can be determined by 1)must reach a fold change that is some percent max of fold all CRs (FC threshold)  
  %2) must reach a PR that is some percent of the max PR for each CR or 3) all
  %the CRs
function [percent_activated,timeDelayOn,timeDelayOff,numberOfCells]=TimeDelay_PR2(mCherry,PR,timeVector,threshold)
timeDelayOn=NaN(1,size(mCherry,1));
timeDelayOff=NaN(1,size(mCherry,1));
activation=PR>=threshold;
deactivated=PR<threshold;
numberOfCells_Activated=sum(activation);
numberOfCells=sum(~isnan(PR));
percent_activated=numberOfCells_Activated./numberOfCells;
%percent_activated(isnan(PR))=NaN;
for i=1:length(timeDelayOn)
    ind=min(find(activation(i,1:28)==1));
    j=0;

    if isempty(ind)==1 || sum(activation(i,:))==1
        timeDelayOn(i)=inf;
    else
    while j==0
        if activation(i,ind+1)==0
            ind=ind+1;
        else
            j=1;
            ind=ind;
        end
    end
        timeDelayOn(i)=timeVector(ind);
    end
    
    ind2=min(find(deactivated(i,ind:end)==1));
    j=0;
    while j==0
        if ind2+ind+1>=52
            ind2=[];
        else
        if deactivated(i,ind2+ind+1)==0
            ind2=ind2+1;
        else
            ind2=ind2;
            j=1;
        end
        end
    end
    
    if isempty(ind2)==1
        timeDelayOff(i)=inf;
    elseif timeDelayOn(i)>=360
        timeDelayOff(i)=NaN;
    elseif timeDelayOn(i)==inf
        timeDelayOff(i)=NaN;
    elseif activation(i,ind2+1)==0
        
    else
        timeDelayOff(i)=timeVector(ind2+ind);
    end
end
figure
%subplot(2,1,1)
histogram(timeDelayOn,-100:20:800)
% hold on
% histogram(timeDelayOff,-100:20:800)
xlabel('Time delay (min)')
xlim([-150,900])
title(strcat('N=', string(max(numberOfCells)),' Percent activated', string(max(percent_activated*100)),'%'))
%
figure
% subplot(2,1,2)
fraction_activated = numberOfCells_Activated./numberOfCells;
plot(timeVector,fraction_activated);
%plot(timeVector,numberOfCells)
percent_activated=nanmax(percent_activated);
%hold on
%plot(timeVector,numberOfCells_Activated)
%plot(timeVector,sum(deactivated))
xlim([-150,900])
ylim([0,1])
legend('Number of cells segmented','Number of cells activated','Number of cells deactivated','Location','northoutside')
xlabel('Time (min)')
ylabel('Number of cells')
end
