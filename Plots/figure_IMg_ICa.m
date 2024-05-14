clear all;

sex = 'male';
ICa_change = [0.0 0.5 1.5];

[male_SS, male_vars, male_SS_norm, male_vars_norm] = getvals(ICa_change, 'male');

% make figures
mcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880;... 
    0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0 0 1];

graymap = gray(6);
darkgray = graymap(2,:);

lw = 3.0;
%ls1 = '-';
f_gca = 18;
fleg = 12;

UP  = char(8593);
DOWN  = char(8595);
w = 0.9;

t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
yline(0.0, 'color', darkgray, 'linewidth', 2.5)
x_SS_1 = categorical(["[PTH]_p"; "[1,25(OH)_2D_3]_p"; "[Mg^{2+}]_p"; "[Ca^{2+}]_p"]);

y_SS_1 = [male_SS_norm(2:end,2).'-1; male_SS_norm(2:end,3).'-1; male_SS_norm(2:end,4).'-1; male_SS_norm(2:end,5).'-1];

bar(y_SS_1, w)
set(gca, 'XTickLabel',x_SS_1, 'XTick',1:numel(x_SS_1))

ylabel('Fractional change in concentration')

legend({'', ['50%' DOWN 'IMg + 50%' DOWN 'ICa'],...
        ['50%' DOWN 'IMg + 50%' UP 'ICa']},...
        'fontsize', fleg, 'location','northeast')
title('(A) Plasma content')
set(gca, 'fontsize', f_gca)
grid on

nexttile
hold on
yline(0.0, 'color', darkgray, 'linewidth', 2.5)
x_SS_2 = categorical(["Gut absorption (Ca^{2+})"; "Gut absorption (Mg^{2+})";...
                    "Bone to plasma (Ca^{2+})"; "Bone to plasma (Mg^{2+})"]);

Ca_bp = (male_vars(2:end,15) + male_vars(2:end,18)) / (male_vars(1,15) + male_vars(1,18));
Mg_bp = (male_vars(2:end,16) + male_vars(2:end,20)) / (male_vars(1,16) + male_vars(1,20));

y_SS_2 = [male_vars_norm(2:end,5).'-1; male_vars_norm(2:end,6).'-1; Ca_bp.'-1; Mg_bp.'-1];
bar(y_SS_2, w)
set(gca, 'XTickLabel',x_SS_2, 'XTick',1:numel(x_SS_2))

ylabel('Fractional change in fluxes')

title('(B) Fluxes into the plasma')
set(gca, 'fontsize', f_gca)
grid on

nexttile
hold on
yline(0.0, 'color', darkgray, 'linewidth', 2.5)
x_SS_3 = categorical(["Urine excretion (Ca^{2+})"; "Urine excretion (Mg^{2+})";...
                       "Plasma to bone ((Ca^{2+})"; "Plasma to bone ((Mg^{2+})"]);

y_SS_3 = [male_vars_norm(2:end,11).'-1; male_vars_norm(2:end,12).'-1; male_vars_norm(2:end,17).'-1; male_vars_norm(2:end,19).'-1];
bar(y_SS_3, w)
set(gca, 'XTickLabel',x_SS_3, 'XTick',1:numel(x_SS_3))

ylabel('Fractional change in fluxes')

title('(C) Fluxes from the plasma')
set(gca, 'fontsize', f_gca)
grid on


function [SSvals, SSvars_vals, SSvals_norm, SSvars_norm] = getvals(ICa_change, sex)
% put data into matrices
SSvals = zeros(length(ICa_change), 9);
SSvars_vals = zeros(length(ICa_change), 29);

for ii = 1:length(ICa_change)
    ICa_change_val = ICa_change(ii);
    fname = strcat('./results_IMg_ICa/', 'ICa_', sex, '_ICachange-', num2str(ICa_change_val), '.mat');
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
end %ii

% make normalized data
SSbase = SSvals(1,:); %101
SSvarsbase = SSvars_vals(1,:);

SSvals_norm = SSvals./SSbase;
SSvars_norm = SSvars_vals./SSvarsbase;
end % getvarsvals