% computes local sensitivity analysis of all parameters for given sexORrep
% variable
clear all

incr = 1; % increase (1) or decrease (0) parameters by 5%

SSfile = './steady_state_baseline/magnesium_mod_SS_male.mat';

save_res = 1; %input('save? (0/1) ');
%notes = input('notes: ');

dat = load(SSfile);
SS_IG = dat.SS;
pars = dat.param_vals;
param_names = dat.param_names;

SSbase = compute_SS(pars, SS_IG, -1, incr);

sens = zeros(length(pars), 4); % col 1: PTHp, col 2: calcitriol, col3: Mg, col4: Ca
frac_change = zeros(size(sens));
basevals = SSbase(2:5)'; % used to compute fractional changes

for ii = 1:66
    SSsens = compute_SS(pars, SS_IG, ii, incr);
    sens(ii, :) = SSsens(2:5);
    frac_change(ii,:) = 100.0*(sens(ii,:) - basevals)./basevals;
end


%% save results
if save_res
    %notes = input('notes: ');
    if incr
        type_change = 'I';
    else
        type_change = 'D';
    end
    fname = strcat('./results_local_sensitivity/localsensitivity_male_',type_change,'.mat');
    save(fname, 'pars', 'frac_change', 'SSbase', 'sens', 'param_names');
    fprintf('sensitivity analysis results saved to %s \n', fname)
end


%%%
function SS = compute_SS(pars, IC, parchange_ID, incr)
    % pars - par values
    % IG - initial guess
    % parchange_ID - par ID to change, set to -1 for normal SS
    params = pars;
    if parchange_ID > -1
        if incr == 1
            f_change = 1.05;
        else
            f_change = 0.95;
        end
        params(parchange_ID) = f_change * pars(parchange_ID);
    end
    
    tspan = [0 4000];
    options = odeset('RelTol',1.0e-6,'AbsTol',1e-6);
    [~, y] = ode15s(@(t,y) magnesium_mod(t,y,params), tspan, IC, options);
    IG = y(end,:)';
    
    options = optimoptions('fsolve', 'Display', 'off', 'MaxFunEvals', 100000, 'MaxIter', 100000);
    [SS, residual, exitflag, ~] = fsolve(@(y) magnesium_mod(0,y,params), IG, options);
    % check between ODE and fsolve result
    Vp = params(4);
    SS(2:5) = SS(2:5)/Vp; % change to concentration
    consODE = IG(2:5)/Vp;
    consfsolve = SS(2:5);
    fracchange = (consODE - consfsolve)./consODE;
    if max(abs(fracchange)) > 0.1
        fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
        fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
    end

    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
end