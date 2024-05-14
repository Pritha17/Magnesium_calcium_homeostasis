clear all;

format longg
run('read_params.m')

sex = 'male';

SSfile = './steady_state_baseline/magnesium_mod_SS_male.mat';
IC = load(SSfile).SS;
Vp = param_vals(4);

ID3_change = (100:-1:0)/100.0;

base_ID3 = param_vals(11);

for i = 1:length(ID3_change)
    ID3_change_val = ID3_change(i);
    ID3_new_val = ID3_change_val * base_ID3;
    param_vals(11) = ID3_new_val;

    % solving ODE with new IMg
    tspan = [0 10000];

    options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
    
    [t,y] = ode15s(@(t,y) magnesium_mod(t,y,param_vals), tspan, IC, options);

    IC_new = y(end,:)';
    %plot(t, y(:,10))
    
    fprintf('solving fsolve \n')
    opts_fsolve = optimoptions('fsolve', 'Display', 'off', 'MaxFunEvals', 100000, 'MaxIter', 100000);
    [SS, fval, exitflag, output] = fsolve(@(y) magnesium_mod(0,y,param_vals), IC_new, opts_fsolve);
    
    consODE = IC_new(2:5)/Vp;
    consfsolve = SS(2:5)/Vp;
    fracchange = (consODE - consfsolve)./consODE;
    if max(abs(fracchange)) > 0.1
        fprintf("%f\n", ID3_change_val)
        fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))
        fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
    end
    
    if exitflag<1
        fprintf('***exitflag indicates error!!*** \n')
    end

    valsSS = compute_fluxes(SS', param_vals);

    % export to .mat files for info
    fname_save = strcat('./results_25D_deficiency/', 'ID3_', sex, '_ID3change-', num2str(ID3_change_val), '.mat');
    save(fname_save, 'SS', 'param_vals', 'sex', 'valsSS', 'ID3_new_val', 'ID3_change_val')
end