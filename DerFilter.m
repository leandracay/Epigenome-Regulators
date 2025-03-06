function [mCherry_filtered,PR_filtered] = DerFilter(mCherry,PR)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
meanPR=nanmean(PR);
upperbound=meanPR+3.*nanstd(PR);


lowerbound=meanPR-3.*nanstd(PR);
  

mCherry_filtered=mCherry;
ind=find(PR<lowerbound | PR>upperbound);
mCherry_filtered(ind)=NaN;
PR_filtered=PR;
PR_filtered(ind)=NaN;
rangeOfCells=1:10;
% figure
% subplot(2,2,1)
% plot(mCherry(rangeOfCells,:)')
% subplot(2,2,2)
% plot(mCherry_filtered(rangeOfCells,:)')
% subplot(2,2,3)
% plot(PR(rangeOfCells,:)')
% subplot(2,2,4)
% plot(PR_filtered(rangeOfCells,:)')
end

