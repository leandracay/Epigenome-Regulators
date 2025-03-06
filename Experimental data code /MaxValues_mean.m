function [maxFC,maxPR]=MaxValues_mean(meanmCherry,meanPR)


  

[FCsorted,I]=sort(meanmCherry,'MissingPlacement','first');
maxFC=nanmedian(FCsorted(end-2:end));

  
[PRsorted,I]=sort(meanPR,'MissingPlacement','first');
maxPR=nanmedian(PRsorted(end-2:end));

  

end
