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
%close all
%excelFileName=input('What is data file name? (make sure to include .xlsx) ','s');
[allData,tex]=xlsread('Summary.xlsx');
Strain=tex(3:end,1);
Round=allData(:,9);
plateNumber=allData(:,10);

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

%Outlier analysis. Determine which samples are outliers for their
%conditions

i=1;
lightMedianNew=zeros(80,length(listOfStrains));
darkMedianNew=zeros(80,length(listOfStrains));
fileNew=num2cell(NaN(80,length(listOfStrains)));
lightOutliers=0;
names=num2cell(NaN(length(listOfStrains),1));
while i<=length(listOfStrains)

  
  tempLightMedian=zeros(80,1);
  tempDarkMedian=zeros(80,1);
  tempLightCV=zeros(80,1);
  tempDarkCV=zeros(80,1);
  %tempLightFiles=num2cell(zeros(80,1));
indices=strcmp(Strain,listOfStrains(i));
  tempLightMedian(1:sum(indices))=lightMedian(indices);
  tempDarkMedian(1:sum(indices))=darkMedian(indices);
  %tempLightFiles(1:sum(Light_indices1))=lightFiles(Light_indices1);
  tempLightCV(1:sum(indices))=lightCV(indices);
  tempDarkCV(1:sum(indices))=darkCV(indices);
  
   lightOutliers=isoutlier(lightMedian(indices));
   lightOutliersCV=isoutlier(lightCV(indices));
   darkOutliers=isoutlier(darkMedian(indices));
   darkOutliersCV=isoutlier(darkCV(indices));
   
   tempLightMedian(lightOutliers)=NaN;
   tempDarkMedian(darkOutliers)=NaN;
   
   tempLightCV(lightOutliersCV)=NaN;
   tempDarkCV(darkOutliersCV)=NaN;
   %tempLightFiles(lightOutliers)={NaN};
   tempLightMedian(tempLightMedian==0)=NaN;
   tempDarkMedian(tempDarkMedian==0)=NaN;
   tempLightCV(tempLightCV==0)=NaN;
   tempDarkCV(tempDarkCV==0)=NaN;
   
   lightMedianNew(:,i)=tempLightMedian;
   darkMedianNew(:,i)=tempDarkMedian;
   %fileNew(:,i)=tempLightFiles;
   lightCVNew(:,i)=tempLightCV;
   darkCVNew(:,i)=tempDarkCV;
   
   %Find the mean for each set of replicates
   averageLightMedian(i)=nanmean(tempLightMedian);
   averageDarkMedian(i)=nanmean(tempDarkMedian);
   averageLightCV(i)=nanmean(tempLightCV);
   averageDarkCV(i)=nanmean(tempDarkCV);
    i=i+1;
end

%find Y11 average
Y11_average=(averageDarkMedian(strcmp(listOfStrains,'Y11'))+averageLightMedian(strcmp(listOfStrains,'Y11')))/2;


LightCVstd=nanstd(lightCVNew);



%foldChange1: fold change by dividing each strain individually and normalizing to the
%strain incubated in the dark
foldChange1=lightMedianNew./(darkMedianNew);
X=categorical(listOfStrains);
foldChange1_average=nanmean(foldChange1);
[foldChange1_sorted,ind]=sort(foldChange1_average);
names_sorted=names(ind);
foldChange1Std=nanstd(foldChange1);

%fluorescence: fold change by dividing each strain individually and normalizing to the
%reporter (JY28)
% Reporter_average=averageDarkMedian(strcmp(listOfStrains,'JY28'));
% fluorescence_light=(lightMedianNew)./(Reporter_average);
% fluorescence_dark=(darkMedianNew)./(Reporter_average);
% fluorescence_light_average=nanmean(fluorescence_light);
% fluorescence_dark_average=nanmean(fluorescence_dark);
% fluorescence_Std=[nanstd(fluorescence_light); nanstd(fluorescence_dark)];
% fluorescence_Std=fluorescence_Std';

%foldChange3: fold change by dividing then averaging to the strain incubated in
%the dark
foldChange3_notaveraged=(lightMedianNew)./(averageDarkMedian);
foldChange3=nanmean(foldChange3_notaveraged);
foldChange3_stdev=nanstd(foldChange3_notaveraged);




Table1=[foldChange1_average',foldChange3',averageDarkCV',averageLightCV'];
% %Table1.Properties.VariableNames={'CR','Fold_Change_1','Dark_FC_2', 'Light_FC_2','Fold_Change_3', 'Dark_CV','Light_CV'};

[b,i]=sort(foldChange3);
bar(1:length(foldChange3),foldChange3(i))
ylim([0,100])
set(gca,'XTick',1:length(foldChange3))
set(gca,'xTickLabel',listOfStrains(i)')
set(gca,'XTickLabelRotation',90)
hold on
errorbar(1:length(foldChange3),foldChange3(i),foldChange3_stdev(i),'linestyle','none','color','k','linewidth',1,'capsize',0);
xlabel('Strain')
ylabel('Fold change')
hold off
figure
lightCVNew_average=nanmean(lightCVNew);
[b,i]=sort(lightCVNew_average);

bar(1:length(lightCVNew_average),lightCVNew_average(i))
ylim([0,3])
xlabel('Strain')
ylabel('Robust CV')
%set(gca,'YScale','log')
set(gca,'XTick',1:length(lightCVNew))
set(gca,'xTickLabel',listOfStrains(i)')
set(gca,'XTickLabelRotation',90)
hold on
errorbar(1:length(lightCVNew_average),lightCVNew_average(i),LightCVstd(i),'linestyle','none','color','k','linewidth',1,'capsize',0);
hold off
figure
text(foldChange3,lightCVNew_average,listOfStrains)
xlim([0,100])
ylim([0,3])
end
