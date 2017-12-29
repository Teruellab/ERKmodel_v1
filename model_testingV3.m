% get dose response behavior of ERK model with and without correlation of
% MEK and ERK
clear all
close all

NumFluc=25;   %number of runs
start_val = 22100000; % Default values for ERK/MEK molecules
%dosN=[10 10.25 10.5 10.75 10.875 11 11.125 11.25 11.375 11.5 11.625 11.75 12 12.25 12.5];
dosN=[11.125];
showpanel=[1];
doses = 2.^(dosN); % Model perturbation: elevated RasGTP (should be in the neighborhood of 2000-20000 molecules
idx = 0;

erk1=0;erk2=0;
ERKh(1:NumFluc)=0;ERKl(1:NumFluc)=0;
mek1=0;mek2=0;
MEKh(1:NumFluc)=0;MEKl(1:NumFluc)=0;

for CoVary=[1 2]   %CoVary = 1 for covariant MEK/ERK, =2 for random MEK/ERK
    figure
    for i = 1:length(doses)
        
        mod1 = 0.15*randn(NumFluc,1); % introduce noise for MEK (lognormal)
        if CoVary==1
            mod2 = mod1; % introduce same noise for ERK if want covariance (lognormal)
        else
            mod2 = 0.15*randn(NumFluc,1); % otherwise noise in ERK is independent from MEK (lognormal)
        end
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
            if sum(i==showpanel)>0
                subplot(length(showpanel),1,find(i==showpanel))
            end
            xx=log2(50000+x(:,1));
            if k<100   % plot first 100 traces
                plot(t,xx)
            end
            if xx(end)>18   %if ppERK is high
                countN=countN+1;
            end
            if CoVary==2 & i>4 & i <8 & xx(end)>18   % collect high ppERK for random MEK/ERK condition
                erk2=erk2+1;
                ERKh(erk2)=start_val*exp(mod2(k));
                mek2=mek2+1;
                MEKh(mek2)=start_val*exp(mod1(k));
            end
            if CoVary==2 & i>4 & i <8 & xx(end)<=16.5  %collect low ppERK for random MEK/ERK condition  
                erk1=erk1+1;
                ERKl(erk1)=start_val*exp(mod2(k));
                mek1=mek1+1;
                MEKl(mek1)=start_val*exp(mod1(k));
            end
            axis([0 4800 14 25])
            x3(i,k,CoVary)=xx(end);
        end
        NN(i)=countN/length(mod1);   %total number of activated cells divided by number of cells
    end
    
    if CoVary==1
        save 'ERKMEKmodel2' NN dosN
    else
        Ehigh=mean(ERKh(1:erk2));
        Elow=mean(ERKl(1:erk1));
        Mhigh=mean(MEKh(1:mek2));
        Mlow=mean(MEKl(1:mek1));
        save 'ERKMEKmodel1' NN dosN Ehigh Elow Mhigh Mlow Ehigh Elow Mhigh Mlow x3
    end
end
% figure
% hold on
% load 'ERKMEKmodel2' NN dosN  %stored data is ERKMEKmodel2 and ERKMEKmodel1
% plot(dosN(2:end-1),NN(2:end-1),'-r')
% load 'ERKMEKmodel1' NN dosN
% plot(dosN(2:end-1),NN(2:end-1),'-b')
% [Ehigh Elow Mhigh Mlow]
% figure,plot(ERKl,MEKl,'g.')
% hold on, plot(ERKh,MEKh,'m.')
% figure,histogram(x3(6,:,1),15:.5:24);hold on;histogram(x3(6,:,2),15:.5:24)
% figure,histogram(x3(8,:,1),15:.5:24);hold on;histogram(x3(8,:,2),15:.5:24)
% figure,histogram(x3(10,:,1),15:.5:24);hold on;histogram(x3(10,:,2),15:.5:24)