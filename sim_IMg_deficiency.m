clear;

format longg
run('read_params.m')

sex = 'male';

SSfile = './steady_state_baseline/magnesium_mod_SS_male.mat';
IC = load(SSfile).SS;
Vp = param_vals(4);

IMg_change = [1, 0.5, 0.25, 0.1];

base_IMg = param_vals(28);

for i = 1:length(IMg_change)
    IMg_change_val = IMg_change(i);
    IMg_new_val = IMg_change_val * base_IMg;
    param_vals(28) = IMg_new_val;

    % solving ODE with new IMg
    tspan = [0 262800]; % 6 months

    options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
    
    [t,y] = ode15s(@(t,y) magnesium_mod(t,y,param_vals), tspan, IC, options);

    IC_new = y(end,:)';

    fprintf('solving fsolve \n')
    opts_fsolve = optimoptions('fsolve', 'Display', 'off', 'MaxFunEvals', 100000, 'MaxIter', 100000);
    [SS, fval, exitflag, output] = fsolve(@(y) magnesium_mod(0,y,param_vals), IC_new, opts_fsolve);
    
    consODE = IC_new(2:5)/Vp;
    consfsolve = SS(2:5)/Vp;
    
    fracchange = (consODE - consfsolve)./consODE;
    if max(abs(fracchange)) > 0.1
        fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
        fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
    end
    
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end
    
    valsSS = compute_fluxes(SS', param_vals);

    % export to .mat files for info
    fname_save = strcat('./results_IMg_deficiency/', 'IMg_', sex, '_IMgchange-', num2str(IMg_change_val), '.mat');
    save(fname_save, 'SS', 'param_vals', 'sex', 'valsSS', 'IMg_new_val', 'IMg_change_val')
end