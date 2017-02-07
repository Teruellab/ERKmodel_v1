function delta = erkOde(t,x,ode_options,v)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% delta = erkOde(t,x,ode_options,v)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Full ERK ODE model - stiff system, designed to be solved using ode15 (see this function's help file to set
% options).  
% In phase 1 (v.PHASE == 1),simulation is run w/o stimulus, until convergence to initialize steady-state,
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set the persistent variable that discrete-delay reactions refer to (commented out if discrete delays aren't used)
%persistent DELAY;
if isempty(t)
%     % The sim function calls the ode function once without arguments to reset persistent variables
%     sz = 5000; % Size of total delay memory - if DELAYS.idx exceeds this, function will throw error
%     DELAY.t = zeros(sz,1);
%     DELAY.nfkb = zeros(sz,1); % increase if getting errors
%     DELAY.ikbat = zeros(sz,1); % increase if getting errors
%     DELAY.idx     = 1; % index (starts at 1)
    return;
end


% Slice parameters, get previous concentrations
p = v.PARAMS;
delta = zeros(size(x));
