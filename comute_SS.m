clear;

format longG
run('read_params.m')

sex = 'male';

do_ode = 1;
dofig = 1;

% set initial conditions
SSfile = './steady_state_baseline/magnesium_mod_SS_male.mat';
SS_IG = load(SSfile).SS;
Vp = param_vals(4);

if do_ode
    tspan = [0 4000];

    IC = SS_IG;
    
    options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
    
    [t,y] = ode15s(@(t,y) magnesium_mod(t,y,param_vals), tspan, IC, options);
    
    if dofig
        figure(1)
        clf
        num_rows = 2; num_cols = 5;
        lw = 2; c1 = [0.3010 0.7450 0.9330];
        c_gray = [0.8500 0.3250 0.0980];
    
        subplot(num_rows,num_cols,1)
        plot(t, y(:,1),'linewidth', lw, 'color', c1)
        %ylim([0, y(end,1)+4])
        ylabel('PTH_g (pmol)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Parathyroid gland PTH', 'FontSize', 14)
        yline(0,'color',c_gray,'linewidth',lw)
        grid on
        
        subplot(num_rows,num_cols,2)
        plot(t, y(:,2)/Vp, 'LineWidth', lw, 'Color', c1)
        ylabel('[PTH]_p (pmol/L)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Plasma PTH', 'FontSize', 14)
        yline(1.5, 'color', c_gray, 'linewidth', lw) % min of normal range
        yline(13, 'color', c_gray,'linewidth', lw) % max of normal range
        grid on
        
        subplot(num_rows,num_cols,3)
        plot(t, y(:,3)/Vp, 'linewidth', lw, 'color', c1)
        ylabel('[1,25(OH)_2D_3]_p (pmol/L)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        yline(70, 'color', c_gray,'linewidth', lw) % min of normal range
        yline(192, 'color', c_gray,'linewidth', lw) % max of normal range
        title('Plasma 1,25(OH)_2D_3', 'FontSize', 14)
        grid on
    
        subplot(num_rows,num_cols,4)
        plot(t, y(:,4)/Vp,'linewidth', lw, 'color', c1)
        %ylim([0.1, 0.9])
        ylabel('[Mg^{2+}]_p (mmol/L)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Plasma Mg^{2+}', 'FontSize', 14)
        yline(0.4, 'color', c_gray,'linewidth', lw) % min of normal range
        yline(0.8, 'color', c_gray,'linewidth', lw) % max of normal range
        grid on
    
        subplot(num_rows,num_cols,5)
        plot(t, y(:,5)/Vp,'linewidth', lw, 'color', c1)
        %ylim([1,1.4])
        ylabel('[Ca^{2+}]_p (mmol/L)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Plasma Ca^{2+}', 'FontSize', 14)
        yline(1.1, 'color', c_gray,'linewidth', lw) % min of normal range
        yline(1.3, 'color', c_gray,'linewidth', lw) % max of normal range
        grid on
    
        subplot(num_rows,num_cols,6)
        plot(t,y(:,6),'linewidth',lw,'color',c1)
        ylabel('NCa_f (mmol)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Fast bone pool (Ca^{2+})', 'FontSize', 14)
        yline(0,'color',c_gray,'linewidth',lw)
        grid on
        
        subplot(num_rows,num_cols,7)
        plot(t,y(:,7),'linewidth',lw, 'color', c1)
        yline(0, 'color',c_gray,'linewidth',lw)
        ylabel('NCa_s (mmol)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Slow bone pool (Ca^{2+})', 'FontSize', 14)
        grid on
        
        subplot(num_rows,num_cols,8)
        plot(t,y(:,8),'linewidth',lw,'color',c1)
        %ylim([0, y(end,9)+0.1])
        ylabel('NMg_f (mmol)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Fast bone pool (Mg^{2+})', 'FontSize', 14)
        yline(0,'color',c_gray,'linewidth',lw)
        grid on
        
        subplot(num_rows,num_cols,9)
        plot(t,y(:,9),'linewidth',lw, 'color', c1)
        yline(0, 'color',c_gray,'linewidth',lw)
        ylabel('NMg_s (mmol)', 'FontSize', 18)
        xlabel('t', 'FontSize', 18)
        title('Slow bone pool (Mg^{2+})', 'FontSize', 14)
        grid on
        
        %sgtitle('Steady states', 'FontSize', 14)
    
    end
    
    SS_IG = y(end,:)';
end

fprintf("fsolve for SS")
options = optimoptions('fsolve', 'Display', 'final', 'MaxFunEvals', 100000, 'MaxIter', 100000);
[SS, residual, exitflag, output] = fsolve(@(y) magnesium_mod(0,y,param_vals), SS_IG, options);

% check between ODE and fsolve result
consODE = SS_IG(:)/Vp;
consfsolve = SS(:)/Vp;
fracchange = (consODE - consfsolve)./consODE;
if max(abs(fracchange)) > 0.1
    fprintf('maximum ODE to fsolve change: %0.3f \n', max(abs(fracchange)))

    fprintf('**WARNING: ODE to fsolve change by more than 10 percent *** \n')
end

if exitflag<1
    fprintf('***exitflag indicates error!!*** \n')
end

fprintf('final steady states \n')
fprintf('           SS\n')
fprintf('PTHg:      %0.3f\n', SS(1))
fprintf('PTHp_con:  %0.3f\n', SS(2)/Vp)
fprintf('D3p_con:   %0.3f\n', SS(3)/Vp)
fprintf('Mgp_con:      %0.3f\n', SS(4)/Vp)
fprintf('Cap_con:      %0.3f\n', SS(5)/Vp)
fprintf('NCaf:      %0.3f\n', SS(6))
fprintf('NCas:      %0.3f\n', SS(7))
fprintf('NMgf:      %0.3f\n', SS(8))
fprintf('NMgs:      %0.3f\n', SS(9))

valsSS = compute_fluxes(SS', param_vals);

save_SS = input('save SS? (0/1) ');
if save_SS
    fname_save = strcat('./steady_state_baseline/','magnesium_mod_SS_', sex, '.mat');
    save(fname_save, 'SS', 'valsSS','param_vals','sex',...
                        'exitflag', 'residual', 'param_names')

    fprintf('results saved to: \n %s \n', fname_save);
end


