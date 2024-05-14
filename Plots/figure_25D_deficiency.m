clear all;

sex = 'male';
ID3_change = (100:-1:0)/100.0;

fig_ID3change = (1 - ID3_change) * 100;

[male_SS, male_vars, male_SS_norm, male_vars_norm] = getvals(ID3_change, 'male');

% make figures
mcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880;... 
    0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0 0 1];

graymap = gray(6);
darkgray = graymap(2,:);

lw = 3.0;
%ls1 = '-';
f_gca = 18;
fleg = 12;
xlab = 'Inhibition of [25(OH)D]_p (%)';

t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
inds = [2 3 4 5];
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
%xline(0.0, 'color', darkgray, 'linewidth', 2.5)

for ii = 1:length(inds)
    ind = inds(ii);
    plot(fig_ID3change, male_SS_norm(:,ind), 'color', mcolors(ii,:),'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized concentration')
ylim([0.0,inf])
legend('', '[PTH]_p', '[1,25(OH)_2D_3]_p', '[Mg^{2+}]_p', '[Ca^{2+}]_p', ...
    'fontsize', fleg, 'location','northwest')
title('(A) Plasma content')
set(gca, 'fontsize', f_gca)
grid on

nexttile
hold on
inds = [5 6];
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
%xline(0.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    ind = inds(ii);
    plot(fig_ID3change, male_vars_norm(:,ind), 'color', mcolors(ii,:),'linewidth',lw)
end
Ca_bp = (male_vars(:,15) + male_vars(:,18)) / (male_vars(1,15) + male_vars(1,18));
Mg_bp = (male_vars(:,16) + male_vars(:,20)) / (male_vars(1,16) + male_vars(1,20));

plot(fig_ID3change, Ca_bp, 'color', mcolors(3,:),'linewidth',lw)
plot(fig_ID3change, Mg_bp, 'color', mcolors(4,:),'linewidth',lw)
xlabel(xlab)
ylabel('Normalized fluxes')

legend('', 'Gut absorption (Ca^{2+})', 'Gut absorption (Mg^{2+})',...
    'Bone to plasma (Ca^{2+})', 'Bone to plasma (Mg^{2+})',...
    'fontsize', fleg, 'location', 'southwest')

title('(B) Fluxes into the plasma')
set(gca, 'fontsize', f_gca)
grid on

nexttile
hold on
inds = [11 12 17 19];
yline(1.0, 'color', darkgray, 'linewidth', 2.5)
%xline(0.0, 'color', darkgray, 'linewidth', 2.5)
for ii = 1:length(inds)
    ind = inds(ii);
    plot(fig_ID3change, male_vars_norm(:,ind), 'color', mcolors(ii,:),'linewidth',lw)
end
xlabel(xlab)
ylabel('Normalized fluxes')
legend('', 'Urine excretion (Ca^{2+})', 'Urine excretion (Mg^{2+})', ...
    'Plasma to bone ((Ca^{2+})', 'Plasma to bone ((Mg^{2+})',...
    'fontsize', fleg, 'location', 'southwest')
title('(C) Fluxes from the plasma')
set(gca, 'fontsize', f_gca)
grid on


function [SSvals, SSvars_vals, SSvals_norm, SSvars_norm] = getvals(ID3_change, sex)
% put data into matrices
SSvals = zeros(length(ID3_change), 9);
SSvars_vals = zeros(length(ID3_change), 29);

for ii = 1:length(ID3_change)
    ID3_change_val = ID3_change(ii);
    fname = strcat('./results_25D_deficiency/', 'ID3_', sex, '_ID3change-', num2str(ID3_change_val), '.mat');
    dat = load(fname);
    SSvals(ii,:) = dat.SS;
    vals = dat.valsSS;
    SSvars_vals(ii,:) =  [vals.PTHg_synthesis;              % 1
                                vals.PTHg_exocytosis;       % 2
                                vals.D3_influx              % 3
                                vals.deg_D3_act             % 4
                                vals.Gut_absorption_Ca      % 5
                                vals.Gut_absorption_Mg;     % 6
                                vals.Renal_filtration_Ca;   % 7
                                vals.Renal_filtration_Mg    % 8
                                vals.Renal_frac_reab_Ca;    % 9
                                vals.Renal_frac_reab_Mg     % 10
                                vals.Urine_excretion_Ca;    % 11
                                vals.Urine_excretion_Mg     % 12
                                vals.Bone_accretion_Ca;     % 13
                                vals.Bone_accretion_Mg      % 14
                                vals.Bone_resorption_Ca;    % 15
                                vals.Bone_resorption_Mg;    % 16
                                vals.Plasma_to_FastPool_Ca; % 17
                                vals.FastPool_to_Plasma_Ca; % 18
                                vals.Plasma_to_FastPool_Mg; % 19
                                vals.FastPool_to_Plasma_Mg; % 20
                                vals.F_Ca_Mg;               % 21
                                vals.Renal_filtration_Mg;   % 22
                                vals.Lambda_PT_Mg;          % 23
                                vals.Lambda_TAL_Mg;         % 24
                                vals.Lambda_DCT_Mg;         % 25
                                vals.Renal_filtration_Ca;   % 26
                                vals.Lambda_PT_Ca;          % 27
                                vals.Lambda_TAL_Ca;         % 28
                                vals.Lambda_DCT_Ca]';       % 29       
    fclose('all');
end

% make normalized data
SSbase = SSvals(1,:);
SSvarsbase = SSvars_vals(1,:);

SSvals_norm = SSvals./SSbase;
SSvars_norm = SSvars_vals./SSvarsbase;
end % getvarsvals