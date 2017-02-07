%% Looking at extrinsic variation: 500 single-cell draws
init_val = 22100000; % Default values for ERK/MEK molecules
num = 500; % Number of simulations
sim_time = 120;

unpaired_dist = init_val*10.^normrnd(0, 0.2, [2, num]);
all_doses = exp(7:.25:10); % Model perturbation: elevated RasGTP (should be in the neihborhood of 2000-20000 molecules)

erk_unpaired = zeros(num,sim_time*60+1,length(all_doses));

% Run dose response for each individual
parfor i = 1:num
    p_mod = [];
    names = {'ERKpp'};
    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 60*sim_time; % Total simulation time (in seconds)
    init_mod = {'MEK',unpaired_dist(1,i); 'ERK', unpaired_dist(2,i)};
    % Simulate all doses (only need to equilibrate on first iteration)
    output = [];
    for j = 1:length(all_doses)
        if isempty(output)
            [t,x,simdata] = erkSimulate({'RasGTP',all_doses(j)},names, p_mod, init_mod,options);
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = erkSimulate({'RasGTP',all_doses(j)}, names, p_mod, init_mod,options);
        end
        output = cat(3,output,x(:)');
    end
    erk_unpaired(i,:,:) = output;

end


% Sort all_dist rows
erk_paired = zeros(num,sim_time*60+1,length(all_doses));

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
    for j = 1:length(all_doses)
        if isempty(output)
            [t,x,simdata2] = erkSimulate({'RasGTP',all_doses(j)},names, p_mod, init_mod,options);
        else
            options.STEADY_STATE = simdata2.STEADY_STATE;
            [~,x] = erkSimulate({'RasGTP',all_doses(j)}, names, p_mod, init_mod,options);
        end
        output = cat(3,output,x(:)');
    end
    erk_paired(i,:,:) = output;

end