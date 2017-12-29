figure
hold on
load 'ERKMEKmodel2_4000' NN dosN  %stored data is ERKMEKmodel2 and ERKMEKmodel1
plot(dosN(2:end-1),NN(2:end-1),'-k')
load 'ERKMEKmodel1_4000' NN dosN Ehigh Elow Mhigh Mlow Ehigh Elow Mhigh Mlow x3
plot(dosN(2:end-1),NN(2:end-1),'-r')
[Ehigh Elow Mhigh Mlow]
load 'ERKMEKmodel3_4000' ERKl ERKh MEKl MEKh
figure,plot(ERKl,MEKl,'g.')
hold on, plot(ERKh,MEKh,'m.')
figure,histogram(x3(3,:,1),15:.25:24);hold on;histogram(x3(3,:,2),15:.25:24)
title('N=3')
figure,histogram(x3(4,:,1),15:.25:24);hold on;histogram(x3(4,:,2),15:.25:24)
title('N=5');
figure,histogram(x3(7,:,1),15:.25:24);hold on;histogram(x3(7,:,2),15:.25:24)
title('N=7')
