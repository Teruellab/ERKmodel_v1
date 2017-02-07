% Quick test: get dose response behavior of ERK model when start values are perturbed in tandem vs individually

init_val = 22100000; % Default values for ERK/MEK molevules
mod1 = -1:.25:1; % Multiplier for ERK/MEK (log10)

ha = tight_subplot(length(mod1),length(mod1));
idx = 0;
for k = mod1
    for j = mod1
        
        idx = idx+1;
        p_mod = [ ];
        init_mod = {'MEK',init_val*10^j; 'ERK', init_val*10^k};
        
        doses = exp(7:.25:10); % Model perturbation: elevated RasGTP (should be in the neihborhood of 2000-20000 molecules)
        
        names = {'ERKpp','MEKpp'};
        options = struct;
        options.DEBUG = 0;
        options.SIM_TIME = 60*80; % Total simulation time (in seconds)

        % Simulate all doses (only need to equilibrate on first iteration)
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
        plot(ha(idx),t,squeeze(output(:,1,:)))
        set(ha(idx),'YTickLabel',{}, 'XTickLabel',{})
        if j==k
            set(ha(idx),'YTickLabel',{}, 'XTickLabel',{},'Color',[0.9098    0.9137    0.9373])
        end
    end
end




