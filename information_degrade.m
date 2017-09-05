% Examining information degradation in the MEK/ERK model
%
% Quick definititions: use single-dimension (AUC) readout for activation. Inputs.... we need a lot, 20+ (32?), spanning
% coverage of a pretty wide dose range.
%
% Simulate with gamma distibutions w/ specified noise - this is more valid for a distribution of abundances, rather
% than concentrations, per se, but still is the best-accepted  model we have.

% Generating Gamma distributions from a "known" CV/ mu: Gamma distribution definitions are:
% mu = a/b
% sigma = a/b^2
% 
% NOTE: MATLAB's a/b notation for Gamma actually matches k/theta notation, not alpha/beta notation
% (a = k; b = 1/theta)
%
% MATLAB's constants can be translated:
% a = k = 1/(CV^2)
% b = theta = mu/k = mu * (CV^2)




%% Looking at extrinsic variation: multiple single-cell draws with Gamma-distributed protein levels

cv = 0.05:.05:.4;
CC_cv = zeros(length(cv),5);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
num = 1024; % Number of simulations
sim_time = 120; % Total length of simulation (minutes)
sim_interval = 15; % Time resoluition of simulation (seconds)
ras_doses = exp(linspace(7,13,32)); % Model perturbation: elevated RasGTP (should be in the neighborhood of 2000-20000 molecules)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Generate distributions for each fixed species in the model
%        Species   Mean       CV
inits = {'Raf'     40000      .3
         'Phase1'  40000      .3
         'MEK'     21000000   .3
         'Phase2'  400000     .3
         'ERK'     21000000*2 .3
         'Phase3'  10000000   .3
         'Phase4'  40000      .3
};



for expt = 1:length(cv)
    gam_dist = @(mu,cv,n) random(makedist('Gamma','a',1/(cv^2),'b',mu*(cv^2)),[1 n]);
    init_dist = [];
    for i = 1:size(inits,1)
        init_dist = cat(1,init_dist, gam_dist(inits{i,2},cv(expt),num));
    end



    erk_unpaired = zeros(num,sim_time*60/sim_interval+1,length(ras_doses)); % Variable to hold simulation output
    % Run dose response for each individual
    parfor i = 1:num
        p_mod = [];
        names = {'ERKpp'};
        options = struct;
        options.DEBUG = 0;
        options.DT = sim_interval; % Time resolution 
        options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
        init_mod = inits(:,1:2);
        init_mod(:,2) = num2cell(init_dist(:,i));
            % Simulate all doses (only need to equilibrate on first iteration)
        output = [];
        for j = 1:length(ras_doses)
            if isempty(output)
                [t,x,simdata] = erkSimulate({'RasGTP',ras_doses(j)},names, p_mod, init_mod, options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = erkSimulate({'RasGTP',ras_doses(j)}, names, p_mod, init_mod,options);
            end
            output = cat(3,output,x(:)');
        end
        erk_unpaired(i,:,:) = output;
    end

    % Calculate channel capacity across all measured responses
    all_auc = cell(size(erk_unpaired,3),1);

    for i = 1:length(all_auc)
        all_auc{i} = nansum(erk_unpaired(:,:,i),2)';
    end
    
    CC_cv(expt,:) = [getCC(all_auc,3,0) getCC(all_auc,5,0) getCC(all_auc,7,0) getCC(all_auc,9,0) getCC(all_auc,11,0)];
end


