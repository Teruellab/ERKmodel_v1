%% Looking at extrinsic variation: 500 single-cell draws
init_val = 22100000*(10^-0.25); % Default values for ERK/MEK molecules - adjusted to match dynamics of MF10A cells
num = 5; % Number of simulations
sim_time = 120;

%        Species   Mean       CV
inits = {'Raf'     40000      .1
         'Phase1'  40000      .1
         'MEK'     21000000   .1
         'Phase2'  400000     .1
         'ERK'     21000000   .1
         'Phase3'  10000000   .1
         'Phase4'  40000      .1
};

init_dist = [];
for i = 1:size(inits,1)
    init_dist = cat(1,init_dist, normrnd(log(inits{i,2})/log(2), log(inits{i,2})/log(2)*inits{i,3}, [1, num]));
end
init_dist = 2.^init_dist; % Undo log transform

[~, pair_order] = sort(init_dist(3,:),'ascend');

init_paired = init_dist;
init_paired(5,:) = sort(init_paired(5,:),'ascend');
init_paired(3,:) = sort(init_paired(3,:),'ascend');



ras_doses = exp(7:.25:9); % Model perturbation: elevated RasGTP (should be in the neihborhood of 2000-20000 molecules)
erk_unpaired = zeros(num,sim_time*60+1,length(ras_doses)); % Variable to hold simulation output

% Run dose response for each individual
for i = 1:num
    p_mod = [];
    names = {'ERKpp'};
    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
    init_mod = inits(:,1:2);
    init_mod(:,2) = mat2cell(init_dist(:,i),ones(size(init_dist(:,1),1),1),1);
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
%%

% Sort all_dist rows
erk_paired = zeros(num,sim_time*60+1,length(ras_doses));

paired_dist = [sort(unpaired_dist(1,:),'ascend');sort(unpaired_dist(2,:),'ascend')];


% Run dose response for each individual
parfor i = 1:num
    p_mod = [];
    names = {'ERKpp'};
    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
    init_mod = {'MEK',paired_dist(1,i); 'ERK', paired_dist(2,i)};
    % Simulate all doses (only need to equilibrate on first iteration)
    output = [];
    for j = 1:length(ras_doses)
        if isempty(output)
            [t,x,simdata2] = erkSimulate({'RasGTP',ras_doses(j)},names, p_mod, init_mod,options);
        else
            options.STEADY_STATE = simdata2.STEADY_STATE;
            [~,x] = erkSimulate({'RasGTP',ras_doses(j)}, names, p_mod, init_mod,options);
        end
        output = cat(3,output,x(:)');
    end
    erk_paired(i,:,:) = output;

end