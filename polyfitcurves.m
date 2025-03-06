[FoldChangeData,names]=xlsread('Doubling time calculations.xlsx');
I50=FoldChangeData(:,19);
I100=FoldChangeData(:,37);
t=FoldChangeData(:,1);
p1 = polyfit(t,I50,7);
p2 = polyfit(t,I100,7);

t1=(0:5:61200);
y1=polyval(p1,t1,7);
y2=polyval(p2,t1,7);

figure(1)
plot(t,I50,'o')
hold on
plot(t1,y1)
hold off

figure(2)
plot (t,I100,'o')
hold on
plot(t1,y2)
hold off
