function [ slope,tslope ] = CalcSlope(fluorescenceMatrixN,time )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sizeMatrix=size(fluorescenceMatrixN);
%smooth curve averaging 15 points
n=1;%cell

for n=1:1:sizeMatrix(1)
%       figure
%   hold on
%   plot(time,fluorescenceMatrixN(n,:)','-k')
%   title(sprintf('Fluorescence for cell %d', n)), xlabel('Time (min)'), ylabel('Fluorescence (au)'), axis([0 450 0 .1]) 
%     %make t and fluorescentMatrix without NaN or 0 values
    fluorescence=fluorescenceMatrixN(n,:);
    t=time(~isnan(fluorescence));
    fluorescence=fluorescence(~isnan(fluorescence));
    Sizefluor=size(fluorescence);
    i=1;%beginning point
k=14;%ending point
if Sizefluor<=14
    sprintf('Cell %d has less than 14 points',n)
    slope(n,1)=0;
    
else
   while k<Sizefluor(2)
   Y=fluorescence(i:k);
   X=t(i:k);
   A=fit(X',Y','poly1');
   B=coeffvalues(A);
   y=B(1)*X+B(2);
   %plot(X,y);
   %1st derivative
   slope(n,i)=B(1);
   i=i+1;
   k=k+1;
   end 
end
end
tslope=t(15:Sizefluor(2));
end
