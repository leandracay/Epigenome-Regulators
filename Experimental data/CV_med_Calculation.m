%%CV_med_Caluclation requires flow cytometry data
%Written by Jessica Lee 20200718
%Updated 20200919 to change input to condition number and use a single
%spreadsheet.
%Modified to analyze data from 20210128 and 20210129. Add file to
%directory.Output is mCherry median',Robust CV of
%mCherry as percent,FSC-A,cell/event count, robust CV of FSC-A as percent,
%CV of mCherry as percent

function [G]=CV_med_Calculation

[allData,names]=xlsread('Summary.xlsx');
Strain=names(3:end,1);
WellID=names(3:end,12);
Round=allData(:,9);
Plate=allData(:,10);

for i=1:length(Round)
   i
   if Round(i)==1
       date=20210128;
       date2='2021-01-28';
   elseif Round(i)==2
           date=20210129;
           date2='2021-01-29';
   end
   if Plate(i)==2
       plateNumber='2.0';
   elseif Plate(i)==3
       plateNumber='3.0';
   elseif Plate(i)==3.1
       plateNumber='3.1';
   end
   S=strcat(string(date),'_',string(plateNumber),'_');
   for j=1:2
       
   if j==1
   Cond_file_name=strcat('Exported flow data/',S,'light/','export_JL',date2,'__',WellID(i),'.csv');
   else
   Cond_file_name=strcat('Exported flow data/',S,'dark/','export_JL',date2,'__',WellID(i),'.csv');
   end
   
if exist(Cond_file_name)==0
    CV(i,j)=NaN;
    Fluor(i,j)=NaN;
else
Y=csvread(Cond_file_name,1,0);
mCherry=Y(:,6);
FSC_A=Y(:,1);
FSC_H=Y(:,2);
SSC_A=Y(:,3);

%mCherry=mCherry./FSC_A; %size correction used for Figures 3, 5,6
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
radius=.7;
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate);
figure
subplot(2,2,1)
scatter(FSC_A,SSC_A,1,'filled')
xlim([0,200])
ylim([0,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,2)
scatter(FSC_A(Gate),SSC_A(Gate),1,'filled')
xlim([0,200])
ylim([0,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,3)
scatter(FSC_A(Gate),FSC_H(Gate),1,'filled')
xlim([0,200])
ylim([0,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,4)
scatter(FSC_A(Both_Gates),FSC_H(Both_Gates),1,'filled')
xlim([0,200])
ylim([0,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
mCherry=mCherry./FSC_A;
ln_mCherry=log(mCherry);
std_ln=std(ln_mCherry);
%CV(i)=sqrt(exp(std_ln.^2)-1);% logarithmic approx
%CV(i)=std_ln/mean(ln_mCherry);
CV(i,j)=std(mCherry(Both_Gates))./mean(mCherry(Both_Gates));
CV_robust(i,j)=.5*(prctile(mCherry(Both_Gates),84.13)-prctile(mCherry(Both_Gates),15.87))/median(mCherry(Both_Gates));%Robust CV
Fluor(i,j)=median(mCherry(Both_Gates));
hypot(i,j)=chi2gof(mCherry(Both_Gates));
count(i,j)=sum(Both_Gates);
CV_FSC(i,j)=.5*(prctile(Y(Both_Gates,1),84.13)-prctile(Y(Both_Gates,1),15.87))/median(Y(Both_Gates,1));%Robust CV
FSC(i,j)=median(FSC_A);
end
   end
end
G_light=table(Fluor(:,1),CV_robust(:,1).*100,FSC(:,1),count(:,1),CV_FSC(:,1).*100,CV(:,1).*100);
writetable(G_light,'light.xlsx')

G_dark=table(Fluor(:,2),CV_robust(:,2).*100,FSC(:,2),count(:,2),CV_FSC(:,2).*100,CV(:,2).*100);
writetable(G_dark,'dark.xlsx')





end