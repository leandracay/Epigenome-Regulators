%Analysis for VP16 matrix experiment. Written on 20200121 by Jessica Lee

%Import data. The Excel file should have the setup of columns as follows:
%mCherry median, robust mCherry CV (in percent), FSC-A median, count, FSC-A CV,
%mCherry CV, FSC-A robust CV. 

%Outputs 
%names: names of CRs 
%Table1=[foldChange1_average',fluorescence_dark_average',fluorescence_light_average',foldChange3',averageDarkCV',averageLightCV'];
%foldChange1:fold change by dividing each strain individually and normalizing to the
%strain incubated in the dark. foldChange1=lightMedianNew./(darkMedianNew);
%fluorescence: normalized to reporter only (JY28)
%foldChange3: fold change by dividing then averaging to the strain incubated in
%the dark. foldChange3_notaveraged=(lightMedianNew)./(averageDarkMedian);

function [names,Table1,foldChange1,fluorescence_light,foldChange3_notaveraged,lightCVNew,darkCVNew,fluorescence_dark,fileNew]=CRAnalysisFunction
close all
colors=[115/255,184/255,191/255;166/255,144/255,33/255;242/255,216/255,87/255;64/255,58/255,47/255;242/255,82/255,46/255];
%excelFileName=input('What is data file name? (make sure to include .xlsx) ','s');
[allData,tex]=xlsread('Summary.xlsx');
Strain=tex(3:end,1);
Round=allData(:,9);
plateNumber=allData(:,10);
OD600=allData(:,8);

[allData,tex]=xlsread('light.xlsx');
mCherryMedian_light=allData(:,1);
mCherryCV_light=allData(:,2)./100;%CV is 6. Robust CV is column 2
amplitude=6e10;%au

[allData,tex]=xlsread('dark.xlsx');
mCherryMedian_dark=allData(:,1);
mCherryCV_dark=allData(:,2)./100;

%Remove samples with less than 250 cells
count=allData(:,4);
mCherryMedian_light(count<=500)=NaN;
mCherryCV_light(count<=500)=NaN;


%Only analyze 1 plate
%mCherryMedian(Date~=20200318 | plateNumber~=3.1)=NaN;
%mCherryCV(Date~=20200117)=NaN;
% mCherryMedian(Date~=20200117)=NaN;r
% Strain(Date~=20200117)={'NaN'};
% mCherryMedian(plateNumber~=2)=NaN;
% mCherryCV(plateNumber~=2)=NaN;
% % Y11_average=1.1425;


listOfStrains=unique(Strain);


% [x,names]=xlsread('Excel Files/CRAnalysis','List of CRs');
% TransName=names(:,1);
% CRName=names(:,2);

%Normalize plates to JY145

listOfRounds=[1, 2];
JY145_ind=strcmp(Strain,'JY145');
Y11_ind=strcmp(Strain,'Y11');
%    
    JY145_light_overall_average=mCherryMedian_light(JY145_ind==1);
    outliers=isoutlier(JY145_light_overall_average);%remove outliers
    JY145_light_overall_average(outliers)=NaN;
    JY145_light_overall_average=nanmean(JY145_light_overall_average);
    JY145_dark_overall_average=mCherryMedian_dark(JY145_ind==1);
    outliers=isoutlier(JY145_dark_overall_average);
    JY145_dark_overall_average(outliers)=NaN;
    JY145_dark_overall_average=nanmean(JY145_dark_overall_average);

for k=1:2
    ind_Date_Y11=find(Round==k & Y11_ind==1);
    Y11_Date_average=nanmean(mCherryMedian_dark(ind_Date_Y11)); %Doesn't matter
    
    ind_Date_p2_light=find(Round==k & plateNumber==2);
    mCherryMedian_light(ind_Date_p2_light)=mCherryMedian_light(ind_Date_p2_light)-Y11_Date_average;
    JY145_Date_p2_light=mCherryMedian_light(Round==k & plateNumber==2 & JY145_ind==1);
    JY145_Date_p2_light(isoutlier(JY145_Date_p2_light))=NaN;
    JY145_Date_p2_light=nanmean(JY145_Date_p2_light);
    mCherryMedian_light(ind_Date_p2_light)=mCherryMedian_light(ind_Date_p2_light)./JY145_Date_p2_light.*JY145_light_overall_average;
    
    ind_Date_p2_dark=find(Round==k & plateNumber==2);
    mCherryMedian_dark(ind_Date_p2_dark)=mCherryMedian_dark(ind_Date_p2_dark)-Y11_Date_average;
    JY145_Date_p2_dark=mCherryMedian_dark(Round==k & plateNumber==2 & JY145_ind==1);
    JY145_Date_p2_dark(isoutlier(JY145_Date_p2_dark))=NaN;
    JY145_Date_p2_dark=nanmean(JY145_Date_p2_dark);
    mCherryMedian_dark(ind_Date_p2_dark)=(mCherryMedian_dark(ind_Date_p2_dark))./JY145_Date_p2_dark.*JY145_dark_overall_average;
    
    ind_Date_p3_light=find(Round==k & plateNumber==3);
    mCherryMedian_light(ind_Date_p3_light)=(mCherryMedian_light(ind_Date_p3_light)-Y11_Date_average);
    JY145_Date_p3_light=mCherryMedian_light(Round==k & plateNumber==3 &  JY145_ind==1);
    JY145_Date_p3_light(isoutlier(JY145_Date_p3_light))=NaN;
    JY145_Date_p3_light=nanmean(JY145_Date_p3_light);
    mCherryMedian_light(ind_Date_p3_light)=(mCherryMedian_light(ind_Date_p3_light))./JY145_Date_p3_light.*JY145_light_overall_average;
    
    ind_Date_p3_dark=find(Round==k & plateNumber==3);
    mCherryMedian_dark(ind_Date_p3_dark)=(mCherryMedian_dark(ind_Date_p3_dark)-Y11_Date_average);
    JY145_Date_p3_dark=mCherryMedian_dark(Round==k & plateNumber==3 & JY145_ind==1);
    JY145_Date_p3_dark(isoutlier(JY145_Date_p3_dark))=NaN;
    JY145_Date_p3_dark=nanmean(JY145_Date_p3_dark);
    mCherryMedian_dark(ind_Date_p3_dark)=(mCherryMedian_dark(ind_Date_p3_dark))./JY145_Date_p3_dark.*JY145_dark_overall_average; 
    
    ind_Date_p31_light=find(Round==k & plateNumber==3.1);
    mCherryMedian_light(ind_Date_p31_light)=(mCherryMedian_light(ind_Date_p31_light)-Y11_Date_average);
    JY145_Date_p31_light=mCherryMedian_light(Round==k & plateNumber==3.1 & JY145_ind==1);
    JY145_Date_p31_light(isoutlier(JY145_Date_p31_light))=NaN;
    JY145_Date_p31_light=nanmean(JY145_Date_p31_light);
    mCherryMedian_light(ind_Date_p31_light)=(mCherryMedian_light(ind_Date_p31_light))./JY145_Date_p31_light.*JY145_light_overall_average;
    
    ind_Date_p31_dark=find(Round==k & plateNumber==3.1);
    mCherryMedian_dark(ind_Date_p31_dark)=(mCherryMedian_dark(ind_Date_p31_dark)-Y11_Date_average);
    JY145_Date_p31_dark=mCherryMedian_dark(Round==k & plateNumber==3.1 & JY145_ind==1);
    JY145_Date_p31_dark(isoutlier(JY145_Date_p31_dark))=NaN;
    JY145_Date_p31_dark=nanmean(JY145_Date_p31_dark);
    mCherryMedian_dark(ind_Date_p31_dark)=(mCherryMedian_dark(ind_Date_p31_dark))./JY145_Date_p31_dark.*JY145_dark_overall_average;

end

lightMedian=mCherryMedian_light;
darkMedian=mCherryMedian_dark;
lightCV=mCherryCV_light;
darkCV=mCherryCV_dark;      
% lightFiles=file(lightIndices);
% darkFiles=file(darkIndices);

%find Y11 average
Y11_average=nanmean([lightMedian(strcmp(Strain,'Y11')); darkMedian(strcmp(Strain,'Y11'))]);




%get CRs
[~,tex]=xlsread('Summary.xlsx','ORF');
%[~,tex2]=xlsread('Summary.xlsx','GO_function');

for i=1:length(Strain)
    CRs(i)=tex(strcmp(tex(:,1),Strain(i,1)),2);
    
end

%fold change 3
foldChange3=(lightMedian-Y11_average)./(darkMedian-Y11_average);
t=table(CRs', foldChange3, lightCV, darkCV, Strain,'VariableNames',{'CR','FC3','CV_light','CV_dark', 'Strain'});

Stats=grpstats(t,{'CR','Strain'},{'mean','sem','meanci'});




figure
[Stats,ind]=sortrows(Stats,7);
CRs=CRs(ind);
bar(1:height(Stats),Stats.mean_CV_light)
ylim([0,3])
xlabel('Strain')
ylabel('Robust CV')
%set(gca,'YScale','log')
set(gca,'XTick',1:height(Stats))
set(gca,'xTickLabel',Stats.CR)
set(gca,'XTickLabelRotation',90)
set(gca,'Yscale','log')
hold on
errorbar(1:height(Stats),Stats.mean_CV_light,Stats.sem_CV_light,'linestyle','none','color','k','linewidth',1,'capsize',0);
set(gcf,'position',[440,225,1000,450]);
hold off

figure
[Stats,ind]=sortrows(Stats,4);

bar(1:height(Stats),Stats.mean_FC3)
ylim([0,150])
set(gca,'XTick',1:height(Stats))
set(gca,'xTickLabel',Stats.CR)
set(gca,'XTickLabelRotation',90)
set(gca,'Yscale','log')
hold on
errorbar(1:height(Stats),Stats.mean_FC3,Stats.sem_FC3,'linestyle','none','color','k','linewidth',1,'capsize',0);
xlabel('Strain')
ylabel('Fold change')
set(gcf,'position',[440,225,1000,450]);
hold off

figure
log_FC3=log10(Stats.mean_FC3(2:88));
log_CV=(Stats.mean_CV_light(2:88));
CRs=Stats.CR(2:88);
FC_1=.75*(log_FC3(end)-log_FC3(1))+log_FC3(1);
FC_2=.5*(log_FC3(end)-log_FC3(1))+log_FC3(1);
FC_3=.25*(log_FC3(end)-log_FC3(1))+log_FC3(1);
FC_4=0;

CVcutoff=range(Stats.mean_CV_light)/2;

hold on
% area(linspace(10^FC_1,170,10),3.*ones(1,10))
% area(linspace(10^FC_2,10^FC_1,10),3.*ones(1,10))
% area(linspace(10^FC_3,10^FC_2,10),3.*ones(1,10))
% area(linspace(10^FC_4,10^FC_3,10),3.*ones(1,10))
% area(linspace(10^FC_1,170,10),CVcutoff.*ones(1,10))
% area(linspace(10^FC_2,10^FC_1,10),CVcutoff.*ones(1,10))
% area(linspace(10^FC_3,10^FC_2,10),CVcutoff.*ones(1,10))
% area(linspace(10^FC_4,10^FC_3,10),CVcutoff.*ones(1,10))
text(Stats.mean_FC3,Stats.mean_CV_light,Stats.CR)
%set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([1,450])
ylim([0,3])
xlabel('Fluorescence')
ylabel('CV')

Group(log_FC3>FC_1 & log_CV>CVcutoff)=1;
Group(log_FC3<=FC_1 & log_FC3>FC_2 & log_CV>CVcutoff)=2;
Group(log_FC3<=FC_2 & log_FC3>FC_3 & log_CV>CVcutoff)=3;
Group(log_FC3<=FC_3 & log_CV>CVcutoff)=4;

Group(log_FC3>FC_1 & log_CV<=CVcutoff)=5;
Group(log_FC3<=FC_1 & log_FC3>FC_2 & log_CV<=CVcutoff)=6;
Group(log_FC3<=FC_2 & log_FC3>FC_3 & log_CV<=CVcutoff)=7;
Group(log_FC3<=FC_3 & log_CV<=CVcutoff)=8;
output=[Stats(2:88,:), table(Group','VariableNames',{'group'})];
output=sortrows(output,13);
list={'hpa2','rad6','nut1','tod6','med8','cpr1','ies1','rtt102','med6','sus1','btt1','chz1','hir3','caf4','bur2'};
for i=1:length(list)
   index=strcmp(Stats.CR,list(i));
   text(Stats.mean_FC3(index),Stats.mean_CV_light(index),list(i),'color',colors(5,:))
end
set(gca,'linewidth',1.5)
% set(gca,'font',8)
set(gcf,'Position',[400,378,500,400])

end
