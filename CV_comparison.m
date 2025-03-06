function [mCherry_CV,PR_CV,maxmCherry_CV,maxPR_CV]=CV_comparison(mCherry,PR,timeVector)
mCherry_CV=NaN(1,length(timeVector));
PR_CV=NaN(1,length(timeVector));
for i=1:size(mCherry,2)
mCherry_CV(i)=.5*(prctile(mCherry(:,i),84.13)-prctile(mCherry(:,i),15.87))/nanmedian(mCherry(:,i));%Robust CV
%mCherry_CV(i)=nanstd(mCherry(:,i))/nanmean(mCherry(:,i));
PR_CV(i)=.5*(prctile(PR(:,i),84.13)-prctile(PR(:,i),15.87))/nanmedian(PR(:,i));%Robust CV
%PR_CV(i)=nanstd(PR(:,i))/nanmean(PR(:,i));
end
figure
subplot(2,1,1)
plot(timeVector,mCherry_CV)
ylabel('CV for fold change')
subplot(2,1,2)
plot(timeVector,PR_CV)
ylabel('CV for PR')
xlabel('Time (min)')
maxmCherry_CV=nanmax(mCherry_CV);
maxPR_CV=nanmax(PR_CV);

end