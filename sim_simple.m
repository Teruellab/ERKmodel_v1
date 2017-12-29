% Basic modeling: simulate dose-response behavior of a single set of parameters


ERK_init = 22100000; % Default values for ERK molecules
MEK_init =  22100000; % Default values for MEK molecules


% Specify all model         
p_mod = [ ];
init_mod = {'MEK',MEK_init; 'ERK', ERK_init};
doses = exp(7:.25:10); % Model perturbation: elevated RasGTP (should be in the neihborhood of 2000-20000 molecules)
names = {'ERKpp','MEKpp'}; % Outputs of model: levels of phosphorylated ERK and MEK

% Specify other model options
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 4800; % Total simulation time (in seconds)

% Simulate all doses (only need to get steady-state equilibration on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = erkSimulate({'RasGTP',doses(i)},names, p_mod, init_mod,options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = erkSimulate({'RasGTP',doses(i)}, names, p_mod, init_mod,options);
    end
    output = cat(3,output,x); 
end

% Main model output: doubly-phosphorylated ERK
figure
plot(t,squeeze(output(:,1,:)))
set(gca,'YTickLabel',{}, 'XTickLabel',{})





