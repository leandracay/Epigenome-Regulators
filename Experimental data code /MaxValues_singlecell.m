
function [maxFC,maxPR,stdFC,stdPR,maxPR_activated]=MaxValues_singlecell(mCherry,PR)

%Select only cells with more than 10 values in between time 100 and 700
%(frames 11-41)
  numberOfValues_mCherry=~isnan(mCherry(:,11:41));
  numberOfValues_PR=~isnan(PR(:,11:41));
  mCherry(sum(numberOfValues_mCherry,2)<10,:)=[];
  PR(sum(numberOfValues_PR,2)<10,:)=[];
  maxFC=NaN(1,size(mCherry,1));
  maxPR=NaN(1,size(PR,1));
  for i=1:size(mCherry)
      
  [FCsorted,I]=sort(mCherry(i,:),'MissingPlacement','first');
maxFC(i)=nanmedian(FCsorted(end-2:end));
  end
  for i=1:size(PR,1)
[PRsorted,I]=sort(PR(i,:),'MissingPlacement','first');
maxPR(i)=nanmedian(PRsorted(end-2:end));

  end
  
  figure
  histogram(maxFC,10)
  xlim([1,2.1])
  ylim([0,150])
  ylabel('number of cells')
  xlabel('maximum fold change')
  stdFC=nanstd(maxFC);
  maxFC=nanmean(maxFC);


  figure
  xlim([1,2.1])
  histogram(maxPR)
  ylabel('maximum PR')
  stdFC=nanstd(maxFC);
  stdPR=nanstd(maxPR);
  maxPR_activated=nanmean(maxPR(maxPR>=3.5e-4));
  maxFC=nanmean(maxFC);
  maxPR=nanmean(maxPR);
  subplot(2,1,1)
  histogram(maxFC,10)
  xlim([1,2.1])
  ylabel('maximum fold change')
  subplot(2,1,2)
  histogram(maxPR)
  ylabel('maximum PR')
  stdFC=nanstd(maxFC);
  stdPR=nanstd(maxPR);
  maxPR_activated=nanmean(maxPR(maxPR>=3.5e-4));
  maxFC=nanmean(maxFC);
  maxPR=nanmean(maxPR);
end
