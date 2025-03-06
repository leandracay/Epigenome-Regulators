close all
clear all
activation_basis=input('Activation basis (1-for Fold change, 2-for PR)');
max_basis=input('Basis for calculating max FC and PR (1-for mean, 2-for single cell)');
colors=[242/255,234/255,47/255;97/255,191/255,122/255;48/255,159/255,217/255;71/255,65/255,154/255];
%close all
[data]=readtable('Summary.xlsx');
Strain=table2array(data(:,1));
nameOfFolder=table2array(data(:,4));
date=table2array(data(:,2));
intensity=table2array(data(:,5));
segmentationStatus=table2array(data(:,8));
backgroundCell=table2array(data(:,12));
%backgroundCell(isnan(backgroundCell))=1; % for now, set the background cell as row 1. Will check this individually later
CR=table2array(data(:,3));
meanmCherry=NaN(length(CR),52);
meanPR=NaN(length(CR),52);
meanmCurvature=NaN(length(CR),52);
meanTimeDelayOn=NaN(length(CR),1);
meanTimeDelayOff=NaN(length(CR),1);
percent_activated=NaN(length(CR),1);
maxFC=NaN(length(CR),1);
maxPR=NaN(length(CR),1);
maxmCherry_CV=NaN(length(CR),1);
maxPR_CV=NaN(length(CR),1);
stdFC=NaN(length(CR),1);
stdPR=NaN(length(CR),1);
numberOfCells=NaN(length(CR),52);

timeVector=20.*(0:1:(52-1))-100;%in min

mCherry_CV=NaN(length(CR),52);
PR_CV=NaN(length(CR),length(timeVector));
maxPR_activated=NaN(length(CR),1);
FC_480=NaN(length(CR),1);
PR_480=NaN(length(CR),1);
CV_480=NaN(length(CR),1);

for i=1:length(Strain)
    nameOfFileOpen=strcat(string(date(i)),'/',string(nameOfFolder(i)),'/mCherry.mat');
    if exist(nameOfFileOpen)
        mCherry=open(nameOfFileOpen);
%mCherry.fluorescenceMatrix=nameOfFileOpen;
%backgroundCell=1;
backgroundmCherry=mCherry.fluorescenceMatrix(backgroundCell(i),:);
mCherry.fluorescenceMatrix(backgroundCell(i),:)=NaN;%also removed background (row 1 if not specified)
mCherrySubBackground=mCherry.fluorescenceMatrix./backgroundmCherry; 
%Remove cells with fewer than 5 frames
numberOfFrames=sum(~isnan(mCherrySubBackground),2);
mCherrySubBackground(numberOfFrames<10,:)=NaN;

%Fold change based on fluorescence before blue light (normalized)
darkmCherry=nanmean(mCherrySubBackground(:,1:5),2);
normalizedmCherry=(mCherrySubBackground)./nanmedian(darkmCherry);
avgdarkmCherry=nanmedian(darkmCherry)
mCherrySubBackground=normalizedmCherry;

     [mCherry,PR,curvature]=Smoothing(mCherrySubBackground);
     mCherry=mCherry(:,1:52);
     PR=PR(:,1:52);
     curvature=curvature(:,1:52);
minFrames=10;
numberOfFrames=sum(~isnan(mCherry),2);
mCherry(numberOfFrames<minFrames,:)=[];
PR(numberOfFrames<minFrames,:)=[];
curvature(numberOfFrames<minFrames,:)=[];
PR(PR==0)=NaN;
% 
meanmCherry(i,1:length(timeVector))=nanmean(mCherry);%add per cell basis
meanPR(i,1:length(timeVector))=nanmean(PR);%add per cell basis
meanCurvature(i,1:length(timeVector))=nanmean(curvature);
upperlimit=80;
if activation_basis==1
    [percent_activated(i),timeDelayOn,timeDelayOff,]=TimeDelay_FC(mCherry,PR,timeVector,1.2);%Based off of fold change
    sgtitle(strcat(CR(i), ' int=', string(intensity(i)),' FC threshold=1.2'))%based off of fold change
elseif activation_basis==2

[percent_activated(i),timeDelayOn,timeDelayOff,numberOfCells(i,:)]=TimeDelay_PR2(mCherry,PR,timeVector,3.5e-4);%Based off of PR
sgtitle(strcat(CR(i), ' int=', string(intensity(i)),' PR threshold=3.5e-4'))%based off of PR change
end
%Removed unactivated ones
stdTimeDelayOn(i)=std(timeDelayOn(timeDelayOn~=inf));
meanTimeDelayOn(i)=nanmean(timeDelayOn(timeDelayOn~=inf));
meanTimeDelayOff(i)=nanmean(timeDelayOff(timeDelayOff~=inf));



% saveas(gcf,strcat('TimePlots/Fold change',string(Strain(i)),'_',string(intensity(i)),'_',string(date(i)),'.tif'))
[~,I]=sort(max(mCherry',[],'omitnan'));
if max_basis==1
[maxFC(i),maxPR(i)]=MaxValues_mean(meanmCherry(i,:),meanPR(i,:));% based on mean values

elseif max_basis==2
  [maxFC(i),maxPR(i),stdFC(i),stdPR(i),maxPR_activated(i)]=MaxValues_singlecell(mCherry,PR);% based on mean values
end
sgtitle(strcat(CR(i), ' int=', string(intensity(i))))
[mCherry_CV(i,:),PR_CV(i,:),maxmCherry_CV(i),maxPR_CV(i)]=CV_comparison(mCherry,PR,timeVector);
sgtitle(strcat(CR(i), ' int=', string(intensity(i))))
maxFC
FC_480(i)=nanmean(mCherry(:,30));
PR_480(i)=nanmean(PR(:,30));
CV_480(i)=nanstd(mCherry(:,30))./nanmean(mCherry(:,30));
Fluor=figure;
set(gcf,'Position',[441,310,913,489])
x=[timeVector,NaN];
z=[mCherry, NaN(size(mCherry,1),1)];

 %subplot(3,2,1)
 %hold on
 %for j=i:size(mCherry,1)
 % colormap(parula)   
 %patch(x,j.*ones(size(x)),z(I(j),:),z(I(j),:),'FaceColor','none','EdgeColor','interp','Linewidth',2)
 %c=colorbar;
 % view(3)
 % end
 % %patch(500,size(mCherry,1)+10,upperlimit,'FaceColor','none','EdgeColor','interp','Linewidth',2)
 % ylim([1,size(mCherry,1)])
 % %set(c,'Limits',[0,upperlimit])
 % xlabel('time (min)')
 % ylabel('cell')
 % zlabel('fluorescence')
 % %zlim([0,upperlimit])
 % hold off
 % 
 % subplot(1,2,1)
 % plot(timeVector,mCherry(1:10,:))
 % hold on
 % %plot(timeVector,nanmean(mCherry),'r','linewidth',2)
 % xlabel('time (min)')
 % ylabel('fluorescence')
 % xlim([-80,960])
 % ylim([1,1.8])
 % zlim([0,upperlimit])
 % hold off
 % 
 % z=[PR,NaN(size(mCherry,1),1)];
 % subplot(3,2,3)
 % for j=i:size(mCherry,1)
 %  colormap(parula)   
 % patch(x,j.*ones(size(x)),z(I(j),:),z(I(j),:),'FaceColor','none','EdgeColor','interp','Linewidth',2)
 % c=colorbar;
 % view(3)
 % end
 % xlabel('time (min)')
 % ylabel('cell')
 % zlabel('PR')
 % %zlim([-.2,.2])
 % %set(c,'Limits',[-.2,.2])
 % hold off
 %  subplot(1,2,2)
 %  plot(timeVector,PR(1:10,:))
 % hold on
 % %plot(timeVector,nanmean(PR),'r','linewidth',2)
 % xlabel('time (min)')
 % ylabel('PR')
 % xlim([-80,960])
 % ylim([-1.5e-3,2e-3])
 % zlim([-.2,.2])
 % hold off

% z=[curvature,NaN(size(mCherry,1),1)];
% subplot(3,2,5)
%  for j=i:size(mCherry,1)
%  colormap(parula)   
%  patch(x,j.*ones(size(x)),z(I(j),:),z(I(j),:),'FaceColor','none','EdgeColor','interp','Linewidth',2)
%  c=colorbar;
%  view(3)
%  end
%  xlabel('time (min)')
%  ylabel('cell')
%  zlabel('curvature')
%  %zlim([-.002,.002])
%  %set(c,'Limits',[-.002,.002])
%  hold off
%  subplot(3,2,6)
%  plot(timeVector,curvature,'color',[.5,.5,.5])
%  hold on
%  plot(timeVector,nanmean(curvature),'r')
%  xlabel('time (min)')
%  ylabel('curvature')
%  zlim([-.002,.002])
%  hold off
 
%  sgtitle(strcat(string(CR(i)),'/',string(intensity(i)),'/',string(date(i))))
%  saveas(gcf,strcat('TimePlots/',string(CR(i)),'_',string(intensity(i)),'_',string(date(i)),'.tif'))
%   else
%      mCherry=NaN; 
  end

end

CR(18)=[];
meanmCherry(18,:)=[];
mCherry_CV(18,:)=[];
intensity(18)=[];
meanPR(18,:)=[];
meanTimeDelayOn(18)=[];%Remove extra VP16
stdTimeDelayOn(18)=[];
meanTimeDelayOff(18)=[];
percent_activated(18)=[];
maxPR(18)=[];
stdPR(18)=[];
maxFC(18)=[];
stdFC(18)=[];
maxmCherry_CV(18)=[];
maxPR_CV(18)=[];
maxPR_activated(18)=[];
FC_480(18)=[];
CV_480(18)=[];
numberOfCells(18,:)=[];

CR(34)=[];
meanmCherry(34,:)=[];
mCherry_CV(34,:)=[];
intensity(34)=[];
meanPR(34,:)=[];
meanTimeDelayOn(34)=[];%Remove extra VP16
stdTimeDelayOn(34)=[];
meanTimeDelayOff(34)=[];
percent_activated(34)=[];
maxPR(34)=[];
stdPR(34)=[];
maxFC(34)=[];
stdFC(34)=[];
maxmCherry_CV(34)=[];
maxPR_CV(34)=[];
maxPR_activated(34)=[];
FC_480(34)=[];
CV_480(34)=[];
numberOfCells(34,:)=[];

CR(52)=[];
meanmCherry(52,:)=[];
mCherry_CV(52,:)=[];
intensity(52)=[];
meanPR(52,:)=[];
meanTimeDelayOn(52)=[];%Remove extra VP16
stdTimeDelayOn(52)=[];
meanTimeDelayOff(52)=[];
percent_activated(52)=[];
maxPR(52)=[];
stdPR(52)=[];
maxFC(52)=[];
stdFC(52)=[];
maxmCherry_CV(52)=[];
maxPR_CV(52)=[];
maxPR_activated(52)=[];
FC_480(52)=[];
CV_480(52)=[];
numberOfCells(52,:)=[];
%% Make heatmap of average
xlabels=num2cell(timeVector);
xlabels(1:2:end)={""};
xlabels(4:4:end)={""};
indices=sum(isnan(meanmCherry),2);
meanmCherry(indices==54,:)=[];
intensity(indices==54)=[];
CR(indices==54)=[];
meanPR(indices==54,:)=[];
meanTimeDelayOn(indices==54)=[];
meanTimeDelayOff(indices==54)=[];
percent_activated(indices==54)=[];
maxPR(indices==54)=[];
stdPR(indices==54)=[];
maxFC(indices==54)=[];
stdFC(indices==54)=[];
maxmCherry_CV(indices==54)=[];
maxPR_CV(indices==54)=[];
maxPR_activated(indices==54)=[];

customColormap=[1,1,1;0.984803921568628,0.961601307189543,0.966013071895425;0.969607843137255,0.923202614379085,0.932026143790850;0.954411764705882,0.884803921568627,0.898039215686275;0.939215686274510,0.846405228758170,0.864052287581699;0.924019607843137,0.808006535947712,0.830065359477124;0.908823529411765,0.769607843137255,0.796078431372549;0.893627450980392,0.731209150326797,0.762091503267974;0.878431372549020,0.692810457516340,0.728104575163399;0.863235294117647,0.654411764705882,0.694117647058824;0.848039215686275,0.616013071895425,0.660130718954248;0.832843137254902,0.577614379084967,0.626143790849673;0.817647058823529,0.539215686274510,0.592156862745098;0.802450980392157,0.500816993464052,0.558169934640523;0.787254901960784,0.462418300653595,0.524183006535948;0.772058823529412,0.424019607843137,0.490196078431373;0.756862745098039,0.385620915032680,0.456209150326797;0.741666666666667,0.347222222222222,0.422222222222222;0.726470588235294,0.308823529411765,0.388235294117647;0.711274509803922,0.270424836601307,0.354248366013072;0.696078431372549,0.232026143790850,0.320261437908497;0.680882352941177,0.193627450980392,0.286274509803922;0.665686274509804,0.155228758169935,0.252287581699346;0.650490196078431,0.116830065359477,0.218300653594771;0.635294117647059,0.0784313725490196,0.184313725490196];
customColormap1=[1,1,1;0.974166666666667,0.989166666666667,0.978333333333333;0.948333333333333,0.978333333333333,0.956666666666667;0.922500000000000,0.967500000000000,0.935000000000000;0.896666666666667,0.956666666666667,0.913333333333333;0.870833333333333,0.945833333333333,0.891666666666667;0.845000000000000,0.935000000000000,0.870000000000000;0.819166666666667,0.924166666666667,0.848333333333333;0.793333333333333,0.913333333333333,0.826666666666667;0.767500000000000,0.902500000000000,0.805000000000000;0.741666666666667,0.891666666666667,0.783333333333333;0.715833333333333,0.880833333333333,0.761666666666667;0.690000000000000,0.870000000000000,0.740000000000000;0.664166666666667,0.859166666666667,0.718333333333333;0.638333333333333,0.848333333333333,0.696666666666667;0.612500000000000,0.837500000000000,0.675000000000000;0.586666666666667,0.826666666666667,0.653333333333333;0.560833333333333,0.815833333333333,0.631666666666667;0.535000000000000,0.805000000000000,0.610000000000000;0.509166666666667,0.794166666666667,0.588333333333333;0.483333333333333,0.783333333333333,0.566666666666667;0.457500000000000,0.772500000000000,0.545000000000000;0.431666666666667,0.761666666666667,0.523333333333333;0.405833333333333,0.750833333333333,0.501666666666667;0.380000000000000,0.740000000000000,0.480000000000000];
 customColormap2=[1,1,1;0.997104073145889,0.997216327687926,0.998801472748493;0.994208146291779,0.994432655375852,0.997602945496986;0.991312219437668,0.991648983063778,0.996404418245479;0.988416292583558,0.988865310751704,0.995205890993973;0.985520365729447,0.986081638439630,0.994007363742466;0.982624438875337,0.983297966127556,0.992808836490959;0.979728512021226,0.980514293815482,0.991610309239452;0.976832585167116,0.977730621503408,0.990411781987945;0.973936658313005,0.974946949191334,0.989213254736438;0.971040731458894,0.972163276879259,0.988014727484932;0.968144804604784,0.969379604567185,0.986816200233425;0.965248877750673,0.966595932255111,0.985617672981918;0.962352950896563,0.963812259943037,0.984419145730411;0.959457024042452,0.961028587630963,0.983220618478904;0.956561097188341,0.958244915318889,0.982022091227397;0.953665170334231,0.955461243006815,0.980823563975890;0.950769243480120,0.952677570694741,0.979625036724384;0.947873316626010,0.949893898382667,0.978426509472877;0.944977389771899,0.947110226070593,0.977227982221370;0.942081462917789,0.944326553758519,0.976029454969863;0.939185536063678,0.941542881446445,0.974830927718356;0.936289609209567,0.938759209134371,0.973632400466849;0.933393682355457,0.935975536822297,0.972433873215343;0.930497755501346,0.933191864510223,0.971235345963836;0.927601828647236,0.930408192198149,0.970036818712329;0.924705901793125,0.927624519886075,0.968838291460822;0.921809974939015,0.924840847574001,0.967639764209315;0.918914048084904,0.922057175261926,0.966441236957808;0.916018121230793,0.919273502949852,0.965242709706301;0.913122194376683,0.916489830637778,0.964044182454795;0.910226267522572,0.913706158325704,0.962845655203288;0.907330340668462,0.910922486013630,0.961647127951781;0.904434413814351,0.908138813701556,0.960448600700274;0.901538486960241,0.905355141389482,0.959250073448767;0.898642560106130,0.902571469077408,0.958051546197260;0.895746633252019,0.899787796765334,0.956853018945753;0.892850706397909,0.897004124453260,0.955654491694247;0.889954779543798,0.894220452141186,0.954455964442740;0.887058852689688,0.891436779829112,0.953257437191233;0.884162925835577,0.888653107517038,0.952058909939726;0.881266998981467,0.885869435204964,0.950860382688219;0.878371072127356,0.883085762892890,0.949661855436712;0.875475145273245,0.880302090580816,0.948463328185206;0.872579218419135,0.877518418268742,0.947264800933699;0.869683291565024,0.874734745956668,0.946066273682192;0.866787364710914,0.871951073644594,0.944867746430685;0.863891437856803,0.869167401332520,0.943669219179178;0.860995511002693,0.866383729020445,0.942470691927671;0.858099584148582,0.863600056708371,0.941272164676164;0.855203657294471,0.860816384396297,0.940073637424658;0.852307730440361,0.858032712084223,0.938875110173151;0.849411803586250,0.855249039772149,0.937676582921644;0.846515876732140,0.852465367460075,0.936478055670137;0.843619949878029,0.849681695148001,0.935279528418630;0.840724023023919,0.846898022835927,0.934081001167123;0.837828096169808,0.844114350523853,0.932882473915616;0.834932169315697,0.841330678211779,0.931683946664110;0.832036242461587,0.838547005899705,0.930485419412603;0.829140315607476,0.835763333587631,0.929286892161096;0.826244388753366,0.832979661275557,0.928088364909589;0.823348461899255,0.830195988963483,0.926889837658082;0.820452535045145,0.827412316651409,0.925691310406575;0.817556608191034,0.824628644339335,0.924492783155068;0.814660681336923,0.821844972027261,0.923294255903562;0.811764754482813,0.819061299715187,0.922095728652055;0.808868827628702,0.816277627403112,0.920897201400548;0.805972900774592,0.813493955091038,0.919698674149041;0.803076973920481,0.810710282778964,0.918500146897534;0.800181047066371,0.807926610466890,0.917301619646027;0.797285120212260,0.805142938154816,0.916103092394521;0.794389193358149,0.802359265842742,0.914904565143014;0.791493266504039,0.799575593530668,0.913706037891507;0.788597339649928,0.796791921218594,0.912507510640000;0.785701412795818,0.794008248906520,0.911308983388493;0.782805485941707,0.791224576594446,0.910110456136986;0.779909559087597,0.788440904282372,0.908911928885479;0.777013632233486,0.785657231970298,0.907713401633973;0.774117705379375,0.782873559658224,0.906514874382466;0.771221778525265,0.780089887346150,0.905316347130959;0.768325851671154,0.777306215034076,0.904117819879452;0.765429924817044,0.774522542722002,0.902919292627945;0.762533997962933,0.771738870409928,0.901720765376438;0.759638071108822,0.768955198097854,0.900522238124932;0.756742144254712,0.766171525785779,0.899323710873425;0.753846217400601,0.763387853473705,0.898125183621918;0.750950290546491,0.760604181161631,0.896926656370411;0.748054363692380,0.757820508849557,0.895728129118904;0.745158436838270,0.755036836537483,0.894529601867397;0.742262509984159,0.752253164225409,0.893331074615890;0.739366583130048,0.749469491913335,0.892132547364384;0.736470656275938,0.746685819601261,0.890934020112877;0.733574729421827,0.743902147289187,0.889735492861370;0.730678802567717,0.741118474977113,0.888536965609863;0.727782875713606,0.738334802665039,0.887338438358356;0.724886948859496,0.735551130352965,0.886139911106849;0.721991022005385,0.732767458040891,0.884941383855343;0.719095095151275,0.729983785728817,0.883742856603836;0.716199168297164,0.727200113416743,0.882544329352329;0.713303241443053,0.724416441104669,0.881345802100822;0.710407314588943,0.721632768792595,0.880147274849315;0.707511387734832,0.718849096480521,0.878948747597808;0.704615460880722,0.716065424168446,0.877750220346301;0.701719534026611,0.713281751856372,0.876551693094795;0.698823607172501,0.710498079544298,0.875353165843288;0.695927680318390,0.707714407232224,0.874154638591781;0.693031753464279,0.704930734920150,0.872956111340274;0.690135826610169,0.702147062608076,0.871757584088767;0.687239899756058,0.699363390296002,0.870559056837260;0.684343972901948,0.696579717983928,0.869360529585753;0.681448046047837,0.693796045671854,0.868162002334247;0.678552119193727,0.691012373359780,0.866963475082740;0.675656192339616,0.688228701047706,0.865764947831233;0.672760265485505,0.685445028735632,0.864566420579726;0.669864338631395,0.682661356423558,0.863367893328219;0.666968411777284,0.679877684111484,0.862169366076712;0.664072484923174,0.677094011799410,0.860970838825205;0.661176558069063,0.674310339487336,0.859772311573699;0.658280631214953,0.671526667175262,0.858573784322192;0.655384704360842,0.668742994863188,0.857375257070685;0.652488777506731,0.665959322551114,0.856176729819178;0.649592850652621,0.663175650239039,0.854978202567671;0.646696923798510,0.660391977926965,0.853779675316164;0.643800996944400,0.657608305614891,0.852581148064658;0.640905070090289,0.654824633302817,0.851382620813151;0.638009143236178,0.652040960990743,0.850184093561644;0.635113216382068,0.649257288678669,0.848985566310137;0.632217289527957,0.646473616366595,0.847787039058630;0.629321362673847,0.643689944054521,0.846588511807123;0.626425435819736,0.640906271742447,0.845389984555616;0.623529508965626,0.638122599430373,0.844191457304110;0.620633582111515,0.635338927118299,0.842992930052603;0.617737655257404,0.632555254806225,0.841794402801096;0.614841728403294,0.629771582494151,0.840595875549589;0.611945801549183,0.626987910182077,0.839397348298082;0.609049874695073,0.624204237870003,0.838198821046575;0.606153947840962,0.621420565557929,0.837000293795069;0.603258020986851,0.618636893245855,0.835801766543562;0.600362094132741,0.615853220933781,0.834603239292055;0.597466167278630,0.613069548621706,0.833404712040548;0.594570240424520,0.610285876309632,0.832206184789041;0.591674313570409,0.607502203997558,0.831007657537534;0.588778386716299,0.604718531685484,0.829809130286027;0.585882459862188,0.601934859373410,0.828610603034521;0.582986533008078,0.599151187061336,0.827412075783014;0.580090606153967,0.596367514749262,0.826213548531507;0.577194679299856,0.593583842437188,0.825015021280000;0.574298752445746,0.590800170125114,0.823816494028493;0.571402825591635,0.588016497813040,0.822617966776986;0.568506898737525,0.585232825500966,0.821419439525480;0.565610971883414,0.582449153188892,0.820220912273973;0.562715045029304,0.579665480876818,0.819022385022466;0.559819118175193,0.576881808564744,0.817823857770959;0.556923191321082,0.574098136252670,0.816625330519452;0.554027264466972,0.571314463940596,0.815426803267945;0.551131337612861,0.568530791628522,0.814228276016438;0.548235410758751,0.565747119316447,0.813029748764931;0.545339483904640,0.562963447004373,0.811831221513425;0.542443557050530,0.560179774692299,0.810632694261918;0.539547630196419,0.557396102380225,0.809434167010411;0.536651703342308,0.554612430068151,0.808235639758904;0.533755776488198,0.551828757756077,0.807037112507397;0.530859849634087,0.549045085444003,0.805838585255890;0.527963922779977,0.546261413131929,0.804640058004384;0.525067995925866,0.543477740819855,0.803441530752877;0.522172069071756,0.540694068507781,0.802243003501370;0.519276142217645,0.537910396195707,0.801044476249863;0.516380215363534,0.535126723883633,0.799845948998356;0.513484288509424,0.532343051571559,0.798647421746849;0.510588361655313,0.529559379259485,0.797448894495342;0.507692434801203,0.526775706947411,0.796250367243836;0.504796507947092,0.523992034635337,0.795051839992329;0.501900581092982,0.521208362323263,0.793853312740822;0.499004654238871,0.518424690011189,0.792654785489315;0.496108727384760,0.515641017699114,0.791456258237808;0.493212800530650,0.512857345387040,0.790257730986301;0.490316873676539,0.510073673074966,0.789059203734794;0.487420946822429,0.507290000762892,0.787860676483288;0.484525019968318,0.504506328450818,0.786662149231781;0.481629093114208,0.501722656138744,0.785463621980274;0.478733166260097,0.498938983826670,0.784265094728767;0.475837239405986,0.496155311514596,0.783066567477260;0.472941312551876,0.493371639202522,0.781868040225753;0.470045385697765,0.490587966890448,0.780669512974247;0.467149458843655,0.487804294578374,0.779470985722740;0.464253531989544,0.485020622266300,0.778272458471233;0.461357605135433,0.482236949954226,0.777073931219726;0.458461678281323,0.479453277642152,0.775875403968219;0.455565751427212,0.476669605330078,0.774676876716712;0.452669824573102,0.473885933018004,0.773478349465205;0.449773897718991,0.471102260705930,0.772279822213699;0.446877970864881,0.468318588393856,0.771081294962192;0.443982044010770,0.465534916081782,0.769882767710685;0.441086117156659,0.462751243769707,0.768684240459178;0.438190190302549,0.459967571457633,0.767485713207671;0.435294263448438,0.457183899145559,0.766287185956164;0.432398336594328,0.454400226833485,0.765088658704658;0.429502409740217,0.451616554521411,0.763890131453151;0.426606482886107,0.448832882209337,0.762691604201644;0.423710556031996,0.446049209897263,0.761493076950137;0.420814629177886,0.443265537585189,0.760294549698630;0.418470628988527,0.440732921207386,0.759218425802940;0.416126628799169,0.438200304829584,0.758142301907251;0.413782628609811,0.435667688451781,0.757066178011561;0.411438628420453,0.433135072073978,0.755990054115871;0.409094628231095,0.430602455696175,0.754913930220182;0.406750628041736,0.428069839318373,0.753837806324492;0.404406627852378,0.425537222940570,0.752761682428802;0.402062627663020,0.423004606562767,0.751685558533113;0.399816804347890,0.420423883516593,0.750596017131279;0.397570981032759,0.417843160470419,0.749506475729445;0.395325157717628,0.415262437424245,0.748416934327610;0.393079334402498,0.412681714378070,0.747327392925776;0.390833511087367,0.410100991331896,0.746237851523942;0.388587687772236,0.407520268285722,0.745148310122108;0.386341864457106,0.404939545239548,0.744058768720274;0.384096041141975,0.402358822193373,0.742969227318440;0.381850217826844,0.399778099147199,0.741879685916606;0.379604394511714,0.397197376101025,0.740790144514772;0.377358571196583,0.394616653054851,0.739700603112938;0.375112747881452,0.392035930008677,0.738611061711104;0.372866924566322,0.389455206962502,0.737521520309269;0.370621101251191,0.386874483916328,0.736431978907435;0.368375277936060,0.384293760870154,0.735342437505601;0.366129454620930,0.381713037823980,0.734252896103767;0.363883631305799,0.379132314777805,0.733163354701933;0.361637807990669,0.376551591731631,0.732073813300099;0.359391984675538,0.373970868685457,0.730984271898265;0.357146161360407,0.371390145639283,0.729894730496431;0.354900338045277,0.368809422593109,0.728805189094597;0.352654514730146,0.366228699546934,0.727715647692762;0.350408691415015,0.363647976500760,0.726626106290928;0.348162868099885,0.361067253454586,0.725536564889094;0.345917044784754,0.358486530408412,0.724447023487260;0.343671221469623,0.355905807362238,0.723357482085426;0.341425398154493,0.353325084316063,0.722267940683592;0.339179574839362,0.350744361269889,0.721178399281758;0.336933751524231,0.348163638223715,0.720088857879924;0.334687928209101,0.345582915177541,0.718999316478090;0.332442104893970,0.343002192131366,0.717909775076255;0.330196281578840,0.340421469085192,0.716820233674421;0.327950458263709,0.337840746039018,0.715730692272587;0.325704634948578,0.335260022992844,0.714641150870753;0.323458811633448,0.332679299946670,0.713551609468919;0.321212988318317,0.330098576900495,0.712462068067085;0.318967165003186,0.327517853854321,0.711372526665251;0.316721341688056,0.324937130808147,0.710282985263417;0.314475518372925,0.322356407761973,0.709193443861583;0.312229695057794,0.319775684715798,0.708103902459748;0.309983871742664,0.317194961669624,0.707014361057914;0.307738048427533,0.314614238623450,0.705924819656080;0.305492225112402,0.312033515577276,0.704835278254246;0.303246401797272,0.309452792531102,0.703745736852412;0.301000578482141,0.306872069484927,0.702656195450578;0.300500289241071,0.303436034742464,0.701328097725289;0.300000000000000,0.300000000000000,0.700000000000000]; 
figure
  subplot(2,2,1)
  
  h = heatmap(timeVector,categorical(CR(intensity==50)),meanmCherry(intensity==50,:));
  %sorty(h,'480','ascend')
  ylabel('CR')
  set(gca,'colormap',customColormap)
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
  set(gca,'ColorLimits',[1,1.5])
    subplot(2,2,3)
  
  h = heatmap(timeVector,categorical(CR(intensity==100)),meanmCherry(intensity==100,:));
  ylabel('CR')
  %sorty(h,'480','ascend')
  set(gca,'colormap',customColormap)
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
  set(gca,'ColorLimits',[1,1.5])
  
  subplot(2,2,2)
  h = heatmap(timeVector,categorical(CR(intensity==50)),meanPR(intensity==50,:));
  %sorty(h,'480','ascend')
  set(gca,'colormap',customColormap1)
  xlabel('Time (min)')
  ylabel('CR')
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
  set(gca,'ColorLimits',[-.0003,.0007])
  
  subplot(2,2,4)
  
  h = heatmap(timeVector,categorical(CR(intensity==100)),meanPR(intensity==100,:));
  %sorty(h,'480','ascend')
  set(gca,'colormap',customColormap1)
  xlabel('Time (min)')
  ylabel('CR')
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
  set(gca,'ColorLimits',[-.0003,.0007])
  
  figure
  subplot(2,2,1)
  h = heatmap(timeVector,categorical(CR(intensity==50)),mCherry_CV(intensity==50,:));
  %sorty(h,'480','ascend')
  set(gca,'colormap',customColormap2)
  xlabel('Time (min)')
  ylabel('CR')
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
   set(gca,'ColorLimits',[0,.12])
  
    subplot(2,2,3)
  h = heatmap(timeVector,categorical(CR(intensity==100)),mCherry_CV(intensity==100,:));
  %sorty(h,'480','ascend')
  set(gca,'colormap',customColormap2)
  xlabel('Time (min)')
  ylabel('CR')
  set(gca,'XDisplayLabels',xlabels)
  set(gca,'FontSize',8)
  sgtitle('CV')
  set(gca,'ColorLimits',[0,.12])
   
       subplot(2,2,2)
   %heatmap(timeVector,categorical(CR(intensity==50)),numberOfCells(intensity==50,:))
   set(gca,'colormap',customColormap1)
   xlabel('Time (min)')
   ylabel('CR')
   %set(gca,'XDisplayLabels',xlabel)
   set(gca,'FontSize',8)
   sgtitle('CV')

   
        subplot(2,2,4)
   heatmap(timeVector,categorical(CR(intensity==0)),mCherry_CV(intensity==0,:))
   set(gca,'colormap',customColormap1)
   xlabel('Time (min)')
   ylabel('CR')
   %set(gca,'XDisplayLabels',xlabel)
   set(gca,'FontSize',8)
   sgtitle('CV')
   set(gca,'ColorLimits',[0,.12])

 
 %% Make graphs for average timeOn delays
  figure 
  bar(1:length(CR(intensity==50)),[meanTimeDelayOn(intensity==50),meanTimeDelayOn(intensity==100)]')
 hold on
 errorbar((1:length(CR(intensity==50)))-.15,meanTimeDelayOn(intensity==50),stdTimeDelayOn(intensity==50),'linestyle','none','color','k','linewidth',2,'capsize',0)
 errorbar((1:length(CR(intensity==100)))+.15,meanTimeDelayOn(intensity==100),stdTimeDelayOn(intensity==100),'linestyle','none','color','k','linewidth',2,'capsize',0)
 set(gca,'XTick',1:length(CR(intensity==100)))
 set(gca,'xTickLabel',CR(intensity==100))
 set(gca,'XTickLabelRotation',90)
 ylabel('max fold change')
 legend('int=50%', 'int=100%')
  %ylim([0,500])
  ylabel('Average time delay (min)')


  figure 
  bar(1:length(CR(intensity==50)),[meanTimeDelayOn(intensity==50)]')
 hold on
 errorbar((1:length(CR(intensity==50))),meanTimeDelayOn(intensity==50),stdTimeDelayOn(intensity==50),'linestyle','none','color','k','linewidth',1,'capsize',4)
 set(gca,'XTick',1:length(CR(intensity==50)))
 set(gca,'xTickLabel',CR(intensity==50))
 set(gca,'XTickLabelRotation',90)
 %legend('int=50%')
  %ylim([0,500])
  ylabel('Average time delay (min)')

  figure 
  bar(1:length(CR(intensity==50)),[maxFC(intensity==50)]')
 hold on
 errorbar((1:length(CR(intensity==50))),maxFC(intensity==50),stdFC(intensity==50),'linestyle','none','color','k','linewidth',1,'capsize',4)
 set(gca,'XTick',1:length(CR(intensity==50)))
 set(gca,'xTickLabel',CR(intensity==50))
 set(gca,'XTickLabelRotation',90)
 %legend('int=50%')
  ylim([1,1.5])
  ylabel('Maxmimum FC (a.u.)')

  figure 
  bar(1:length(CR(intensity==50)),[percent_activated(intensity==50)]')
 hold on
 set(gca,'XTick',1:length(CR(intensity==50)))
 set(gca,'xTickLabel',CR(intensity==50))
 set(gca,'XTickLabelRotation',90)
 %legend('int=50%')
  %ylim([0,500])
  ylabel('Fraction Activated')
  %% 
  

 %Make graphs for % activated
  figure
  bar(1:length(CR(intensity==50)),[percent_activated(intensity==50),percent_activated(intensity==100)]')
  ylim([0,1.1])
   set(gca,'XTick',1:length(CR(intensity==100)))
 set(gca,'xTickLabel',CR(intensity==100))
 set(gca,'XTickLabelRotation',90)
  ylabel('Fraction activated')
legend('int=50','int=100')

figure
bar(categorical(CR(intensity==50)),[meanTimeDelayOff(intensity==50),meanTimeDelayOff(intensity==100)]')
ylabel('Average time delay off (min)')

figure
bar(categorical(CR(intensity==50)),[meanTimeDelayOff(intensity==50)]')
ylabel('Average time delay off (min)')
%% 
figure
bar(categorical(CR(intensity==50)),[meanTimeDelayOff(intensity==50),meanTimeDelayOff(intensity==100)]')
ylabel('Average time delay off (min)')



 %% Compare results to endpoint data

  %Step 1: Find maximum FC and PR for each CR and for all CRs
  %Divide this up by intensity
  %Intensity=50%
  
  [maxFC_sorted,I]=sort(maxFC);
  intensity_sorted=intensity(I);
  err=stdFC(I);

 
 figure
     subplot(1,2,1)
         bar(1:length(CR(intensity_sorted==50)),maxFC_sorted(intensity_sorted==50))
         ylim([.9,1.8])
         set(gca,'XTick',1:length(CR(intensity_sorted==50)))
         set(gca,'xTickLabel',CR(I(intensity_sorted==50)))
         set(gca,'XTickLabelRotation',90)
         set(gca,'Yscale','log')
         title('Intensity=50%')
         ylabel('Maximum fold change')
         hold on
         errorbar(1:length(CR(intensity_sorted==50)),maxFC_sorted(intensity_sorted==50),err(intensity_sorted==50),'linestyle','none','color','k','linewidth',2,'capsize',0);
    %Intensity=100%
    subplot(1,2,2)
         bar(1:length(CR(intensity_sorted==100)),maxFC_sorted(intensity_sorted==100))
         ylim([.9,1.8])
         set(gca,'XTick',1:length(CR(intensity_sorted==100)))
         set(gca,'xTickLabel',CR(I(intensity_sorted==100)))
         set(gca,'XTickLabelRotation',90)
         set(gca,'Yscale','log')
         title('Intensity=100%')
         hold on
         errorbar(1:length(CR(intensity_sorted==100)),maxFC_sorted(intensity_sorted==100),err(intensity_sorted==100),'linestyle','none','color','k','linewidth',2,'capsize',0);

 [maxPR_sorted,I]=sort(maxPR);
 intensity_sorted=intensity(I);
 err=stdPR(I);
  figure
    subplot(1,2,1)
  bar(1:length(CR(intensity_sorted==50)),maxPR_sorted(intensity_sorted==50))
  hold on
  errorbar(1:length(CR(intensity_sorted==50)),maxPR_sorted(intensity_sorted==50),err(intensity_sorted==50),'linestyle','none','color','k','linewidth',2,'capsize',0);
ylim([10^-5,3*10^-3])
 set(gca,'XTick',1:length(CR(intensity_sorted==50)))
 set(gca,'xTickLabel',CR(I(intensity_sorted==50)))
 set(gca,'XTickLabelRotation',90)
 set(gca,'Yscale','log')
 title('Intensity=50%')
 ylabel('Maximum PR')
 %Intensity=100%
   subplot(1,2,2)
 bar(1:length(CR(intensity_sorted==100)),maxPR_sorted(intensity_sorted==100))
 hold on
 errorbar(1:length(CR(intensity_sorted==100)),maxPR_sorted(intensity_sorted==100),err(intensity_sorted==100),'linestyle','none','color','k','linewidth',2,'capsize',0);
ylim([10^-5,3*10^-3])
 set(gca,'XTick',1:length(CR(intensity_sorted==100)))
 set(gca,'xTickLabel',CR(I(intensity_sorted==100)))
 set(gca,'XTickLabelRotation',90)
 set(gca,'Yscale','log')
 title('Intensity=100%')
 figure
 bar(1:length(CR(intensity==50)),[maxFC(intensity==50),maxFC(intensity==100)]')
 hold on
 errorbar((1:length(CR(intensity==50)))-.15,maxFC(intensity==50),stdFC(intensity==50),'linestyle','none','color','k','linewidth',2,'capsize',0)
 errorbar((1:length(CR(intensity==100)))+.15,maxFC(intensity==100),stdFC(intensity==100),'linestyle','none','color','k','linewidth',2,'capsize',0)
 set(gca,'XTick',1:length(CR(intensity==100)))
 set(gca,'xTickLabel',CR(intensity==100))
 set(gca,'XTickLabelRotation',90)
 ylabel('max fold change')
 legend('int=50%', 'int=100%')
 set(gca,'Yscale','log')

 figure
 bar(1:length(CR(intensity==50)),[maxPR(intensity==50),maxPR(intensity==100)]')
 hold on
 errorbar((1:length(CR(intensity==50)))-.15,maxPR(intensity==50),stdPR(intensity==50),'linestyle','none','color','k','linewidth',2,'capsize',0)
 errorbar((1:length(CR(intensity==100)))+.15,maxPR(intensity==100),stdPR(intensity==100),'linestyle','none','color','k','linewidth',2,'capsize',0)
 set(gca,'XTick',1:length(CR(intensity==100)))
 set(gca,'xTickLabel',CR(intensity==100))
 set(gca,'XTickLabelRotation',90)
 ylabel('max PR')
 legend('int=50%', 'int=100%')
 %set(gca,'Yscale','log')
 %make max CV figures
 figure
   subplot(2,1,1)
 bar(categorical(CR(intensity==50)),[maxmCherry_CV(intensity==50),maxmCherry_CV(intensity==100)]');
 legend('Int=50','Int=100')
 ylabel('maximum mCherry CV')
   subplot(2,1,2)
 bar(categorical(CR(intensity==50)),[maxPR_CV(intensity==50),maxPR_CV(intensity==100)]');
 legend('Int=50','Int=100')
 ylabel('maximum PR CV')
 figure
 % percent activated vs maxFC
   subplot(1,2,1)
 scatter(maxFC(intensity==50),percent_activated(intensity==50),'k.')
 text(maxFC(intensity==50)+.01,percent_activated(intensity==50),CR(intensity==50),'k')
 ylabel('Fraction activated')
 xlabel('Maximum fold change')
 xlim([0,1.6])
   %splotubplot(1,2,2)
 scatter(maxFC(intensity==100),percent_activated(intensity==100),'k.')
 text(maxFC(intensity==100)+.01,percent_activated(intensity==100),CR(intensity==100),'b')
 xlabel('Maximum fold change')
 xlim([0,1.6])
 figure
 %max PR for activated cells

 bar(1:length(CR(intensity==50)),[maxPR_activated(intensity==50),maxPR_activated(intensity==100)]')
 hold on
  errorbar((1:length(CR(intensity==50)))-.15,maxPR(intensity==50),stdPR(intensity==50),'linestyle','none','color','k','linewidth',2,'capsize',0)
  errorbar((1:length(CR(intensity==100)))+.15,maxPR(intensity==100),stdPR(intensity==100),'linestyle','none','color','k','linewidth',2,'capsize',0)
 set(gca,'XTick',1:length(CR(intensity==100)))
 set(gca,'xTickLabel',CR(intensity==100))
 set(gca,'XTickLabelRotation',90)
 ylabel('Max PR for activated cells')
 legend('int=50%', 'int=100%')
 
 %CVvs FC at 480 min
 figure
   %subplot(1,1)
 fontsize(14,"points")
 scatter(FC_480(intensity==50),CV_480(intensity==50));
 text(FC_480(intensity==50),CV_480(intensity==50),CR(intensity==50))
 ylabel('CV at 480 min')
 xlabel('FC at 480 min')
 ylim([0,0.14])
 title('int=50%')
   %subplot(1,2)
   figure
 scatter(FC_480(intensity==100),CV_480(intensity==100));
 text(FC_480(intensity==100),CV_480(intensity==100),CR(intensity==100))
 title('int=100%')

%% Save all figures
 savestatus=input('Do you want to save all the figures? 1-Yes, 2-No');
  if savestatus==1
      if activation_basis==2 && max_basis==1
  FolderName = 'TimePlots/PR_activation_mean_max';   % Your destination folder
      elseif activation_basis==2 && max_basis==2
          FolderName='TimePlots/PR_activation_SC_max';
      elseif activation_basis==1 && max_basis==2
          FolderName='TimePlots/FC_activation_SC_max';
      elseif activation_basis==1 && max_basis==1
          FolderName='TimePlots/FC_activation_mean_max';
      end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = strcat('Fig',string(iFig),'.fig');
  saveas(FigHandle, fullfile(FolderName, FigName));
end
  else
  end