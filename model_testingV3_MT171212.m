% get dose response behavior of ERK model with and without correlation of
% MEK and ERK
clear all
close all

NumFluc=50;   %number of runs
start_val = 22100000; % Default values for ERK/MEK molecules
%dosN=[10 10.25 10.5 10.75 10.875 11 11.125 11.25 11.375 11.5 11.625 11.75 12 12.25 12.5];
dosN=[11.1];
showpanel=[1];
doses = 2.^(dosN); % Model perturbation: elevated RasGTP (should be in the neighborhood of 2000-20000 molecules
idx = 0;

erk1=0;erk2=0;
ERKh(1:NumFluc)=0;ERKl(1:NumFluc)=0;
mek1=0;mek2=0;
MEKh(1:NumFluc)=0;MEKl(1:NumFluc)=0;

for CoVary=[2]   %CoVary = 1 for covariant MEK/ERK, =2 for random MEK/ERK
    figure
    for i = 1:length(doses)
        mod1 = 0.15*randn(NumFluc,1); % introduce noise for MEK (lognormal)
        mod2 = 0.15*randn(NumFluc,1); % otherwise noise in ERK is independent from MEK (lognormal)
        mod3= 0.1 *randn(NumFluc,1); % Multiplier for upstream noise from ERK/MEK (lognormal)
        idx = idx+1;
        p_mod = [ ];
        
        names = {'ERKpp','MEKpp'};
        options = struct;
        options.DEBUG = 0;
        options.SIM_TIME = 60*80; % Total simulation time (in seconds)
        
        % Simulate all doses (only need to equilibrate on first iteration)
        output = [];
        countN=0;
        
        for k=1:length(mod1)
            init_mod = {'MEK',start_val*exp(mod1(k)); 'ERK', start_val*exp(mod2(k))};   %add noise in real space
            [t,x,simdata] = erkSimulate({'RasGTP',exp(mod3(k))*doses(i)},names, p_mod, init_mod,options); %add noise in real space
            hold on
            xx(:,k)=log2(50000+x(:,1));
%             if k<100   % plot first 100 traces
%                 plot(t,xx,'r')
%             end
           % axis([0 4800 14 25])
        end
    end
end

figure
for i=1:50
    hold on
    plot(t,xx(:,i),'r')
end
%save 'out4' t, xx;

load 'out4', t, xx;
figure
for i=1:30
    hold on
    plot(t,xx(:,i),'r')
end

figure
hold on
load 'ERKMEKmodel2' NN dosN  %stored data is ERKMEKmodel2 and ERKMEKmodel1
plot(dosN(2:end-1),NN(2:end-1),'-r')
load 'ERKMEKmodel1' NN dosN
plot(dosN(2:end-1),NN(2:end-1),'-b')
[Ehigh Elow Mhigh Mlow]
figure,plot(ERKl,MEKl,'g.')
hold on, plot(ERKh,MEKh,'m.')
figure,histogram(x3(5,:,1),15:.5:24);hold on;histogram(x3(6,:,2),15:.5:24)
figure,histogram(x3(9,:,1),15:.5:24);hold on;histogram(x3(9,:,2),15:.5:24)
figure,histogram(x3(10,:,1),15:.5:24);hold on;histogram(x3(10,:,2),15:.5:24)