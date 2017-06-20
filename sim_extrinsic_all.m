%% Looking at extrinsic variation: multiple single-cell draws with log-normally distributed behavior

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
num = 1024; % Number of simulations
sim_time = 120;
ras_doses = exp(7:.25:11); % Model perturbation: elevated RasGTP (should be in the neighborhood of 2000-20000 molecules)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Generate distributions for each fixed species in the model
%        Species   Mean       CV
inits = {'Raf'     40000      .04
         'Phase1'  40000      .04
         'MEK'     21000000   .04
         'Phase2'  400000     .04
         'ERK'     21000000   .04
         'Phase3'  10000000   .04
         'Phase4'  40000      .04
};
init_dist = [];
for i = 1:size(inits,1)
    init_dist = cat(1,init_dist, normrnd(log(inits{i,2})/log(2), log(inits{i,2})/log(2)*inits{i,3}, [1, num]));
end
init_dist = 2.^init_dist; % Undo log transform
% Reorder MEK/ERK levels to correlate them
init_paired = init_dist;
init_paired(5,:) = sort(init_paired(5,:),'ascend');
init_paired(3,:) = sort(init_paired(3,:),'ascend');


erk_unpaired = zeros(num,sim_time*60+1,length(ras_doses)); % Variable to hold simulation output
% Run dose response for each individual
parfor i = 1:num
    p_mod = [];
    names = {'ERKpp'};
    options = struct;
    options.DEBUG = 0;
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



erk_paired = zeros(num,sim_time*60+1,length(ras_doses)); % Variable to hold simulation output
% Run dose response for each individual
parfor i = 1:num
    p_mod = [];
    names = {'ERKpp'};
    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
    init_mod = inits(:,1:2);
    init_mod(:,2) = num2cell(init_paired(:,i));
        % Simulate all doses (only need to equilibrate on first iteration)
    output = [];
    for j = 1:length(ras_doses)
        if isempty(output)
            [t,x,simdata2] = erkSimulate({'RasGTP',ras_doses(j)},names, p_mod, init_mod, options);
        else
            options.STEADY_STATE = simdata2.STEADY_STATE;
            [~,x] = erkSimulate({'RasGTP',ras_doses(j)}, names, p_mod, init_mod,options);
        end
        output = cat(3,output,x(:)');
    end
    erk_paired(i,:,:) = output;
end