clear all;

format longg
run('read_params.m')

sex = 'male';

SSfile = './steady_state_baseline/magnesium_mod_SS_male.mat';
IC = load(SSfile).SS;
Vp = param_vals(4);

% Select which simulation to perform
only_low_IMg = 0;
ICa_incr = 1; % 1: increase ICa by 50%; 0: decrease ICa by 50%

base_IMg = param_vals(28);
IMg_new_val = 0.5 * base_IMg;
param_vals(28) = IMg_new_val;

if only_low_IMg == 0
    if ICa_incr == 1
        ICa_change = 1.5;
    else
        ICa_change = 0.5;
    end
    base_ICa = param_vals(24);
    ICa_new_val = ICa_change * base_ICa;
    param_vals(24) = ICa_new_val;
end


% solving ODE with new IMg
tspan = [0 10000];

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
    
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

if only_low_IMg == 1
    sim_type = '0';
elseif ICa_incr == 1
    sim_type = '1.5';
else
    sim_type = '0.5';
end

% export to .mat files for info
fname_save = strcat('./results_IMg_ICa/', 'ICa_', sex, '_ICachange-', sim_type, '.mat');
save(fname_save, 'SS', 'param_vals', 'sex', 'valsSS')
