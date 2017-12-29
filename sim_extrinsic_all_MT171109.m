% General modeling parameters
num = 1024; % Number of "single-cell" simulations
sim_time = 120;
all_doses = exp(7:.25:10); % Model perturbation: elevated RasGTP (should be in the neighborhood of 2000-20000 molecules)

% Get all pairwise combinations of steady-state species
all_species = {'Raf','Phase1','MEK','Phase2','ERK','Phase3','Phase4'};
all_init = [40000 40000 21000000 400000 22100000 10000000 40000];
combos = nchoosek(1:length(all_init),2);

all_erk_unpaired = cell(size(combos,1),1);
all_erk_paired = cell(size(combos,1),1);
all_distr1 = cell(size(combos,1),1);
all_distr2 = cell(size(combos,1),1);
all_distr1a = cell(size(combos,1),1);
all_distr1b = cell(size(combos,1),1);
%%
% for z = 1:size(combos,1)
for z=1:1
    p_idx1 = combos(z,1);
    p_idx2 = combos(z,2);
    
    % Generate according distributions for the two species being compared.
    %distr1 = normrnd(log(all_init(p_idx1))/log(2), log(all_init(p_idx1))/log(2)*0.04, [1, num]);
    distr1 = normrnd(log(all_init(p_idx1)), log(all_init(p_idx1))*0.05, [1, num]);
    distr1a = exp(distr1); %Transform back into "real values"
    distr1b = exp.^(distr1);
    figure, hist(distr1,100);
    figure, hist(distr1a,100);
    figure, hist(distr1b,100);
%     distr2 = normrnd(log(all_init(p_idx2))/log(2), log(all_init(p_idx2))/log(2)*0.05, [1, num]);
%     distr2a = exp(distr2); %Transform back into "real values"
end
% 
%     % 1) Simulate UNPAIRED parameters - keep randomized order from distributions.
%     erk{z} = zeros(num,sim_time*60+1,length(all_doses));
%     parfor i = 1:num
%         p_mod = [];
%         names = {'ERKpp'};
%         options = struct;
%         options.DEBUG = 0;
%         options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
%         init_mod = {all_species{p_idx1},distr1(i); all_species{p_idx2}, distr2(i)};
%         % Simulate all doses (only need to equilibrate on first iteration)
%         output = [];
%         for j = 1:length(all_doses)
%             if isempty(output)
%                 [t,x,simdata] = erkSimulate({'RasGTP',all_doses(j)},names, p_mod, init_mod,options);
%             else
%                 options.STEADY_STATE = simdata.STEADY_STATE;
%                 [~,x] = erkSimulate({'RasGTP',all_doses(j)}, names, p_mod, init_mod,options);
%             end
%             output = cat(3,output,x(:)');
%         end
%         erk(i,:,:) = output;
%     end
%     all_erk_unpaired{z} = erk;
%     all_distr1{z} = distr1;
%     all_distr2{z} = distr2;
%     
%     % 2) Simulate PAIRED parameters - sort distr1 and distr2 in ascending order, forcing perfect Spearman correlation.
%     distr1 = sort(distr1,'ascend');
%     distr2 = sort(distr2,'ascend');
%     erk{z} = zeros(num,sim_time*60+1,length(all_doses));
%     parfor i = 1:num
%         p_mod = [];
%         names = {'ERKpp'};
%         options = struct;
%         options.DEBUG = 0;
%         options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
%         init_mod = {all_species{p_idx1},distr1(i); all_species{p_idx2}, distr2(i)};
%         % Simulate all doses (only need to equilibrate on first iteration)
%         output = [];
%         for j = 1:length(all_doses)
%             if isempty(output)
%                 [t,x,simdata2] = erkSimulate({'RasGTP',all_doses(j)},names, p_mod, init_mod,options);
%             else
%                 options.STEADY_STATE = simdata2.STEADY_STATE;
%                 [~,x] = erkSimulate({'RasGTP',all_doses(j)}, names, p_mod, init_mod,options);
%             end
%             output = cat(3,output,x(:)');
%         end
%         erk(i,:,:) = output;
%     end
%     all_erk_paired{z} = erk;    
% end
