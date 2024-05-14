% postprocess output from all the local sensitivity analysis
% output from compute_localsensitivity.m
clear all;

incr = 1; % increase (1) or decrease (0) parameters by 5%

% male
if incr
    type_change = 'I';
else
    type_change = 'D';
end
fname = strcat('./results_local_sensitivity/localsensitivity_male_',type_change,'.mat');
male_dat = load(fname);

[male_xlabs, male_vals] = getdata(male_dat);

xlabels = male_xlabs;

% Make figures
temp = male_vals;

clims = [min(temp, [], 'all'), max(temp, [], 'all')];
ylabels = {'[PTH]_p', '[1,25(OH)_2D_3]_p', '[Mg^{2+}]_p', '[Ca^{2+}]_p'};

fsize = 13;
colmap = parula;
cmissdat = 'w';
labmissdat = '<1.0%';

% figure without removing 
fig_noremove = 0;
if fig_noremove
figure(1)
clf
h1 = heatmap(xlabels, ylabels, male_vals,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'off');
h1.FontSize = fsize;
end

% remove Nan values
Nan_male = find_all_Nan(male_vals);

[xlabs_rm, male_vals_rm] = removeAllNan(Nan_male, xlabels, male_vals);

% figure with removing 
fig_remove = 1;
if fig_remove
f2 = figure(2);
f2.Position = [1 1 3000 600];
clf
h1 = heatmap(xlabs_rm, ylabels, male_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'ColorLimits', clims, 'colorbarvisible', 'on');
h1.FontSize = fsize;

sgtitle('(A)', 'fontsize', 24)
end

%------------------
% functions used
%------------------
function [xlabels, round_data] = getdata(dat)
    frac_sens = dat.frac_change;
    labels = cell(size(dat.param_names));
    for ii = 1:size(dat.param_names, 2)
        labels{ii} = param_name_conversion(string(dat.param_names{ii}));
    end
    xlabels = labels;
    round_data = round(frac_sens', 2, 'significant');

    [r,c] = find(abs(round_data) <= 1.0);% r - row values, c - column value
    for ii = 1:length(r)
        round_data(r(ii),c(ii)) = NaN;
    end
end

function AllNan_vals = find_all_Nan(round_data)
    % finds which indices are all Nan values
    PTH_Nan = find(isnan(round_data(1,:)));
    D3_Nan  = find(isnan(round_data(2,:)));
    temp1 = intersect(PTH_Nan, D3_Nan);
    Mg_Nan  = find(isnan(round_data(3,:)));
    temp2 = intersect(temp1, Mg_Nan);
    Ca_Nan  = find(isnan(round_data(4,:)));
    AllNan_vals = intersect(temp2, Ca_Nan);
end

function [rmNan_xlabels, rmNan_round_data] = removeAllNan(AllNan_vals, xlabels, round_data)
    % removes columns listed in indices from AllNan_vals
    rmNan_round_data = round_data;
    rmNan_xlabels = xlabels;
    
    rmNan_round_data(:,AllNan_vals) = [];
    for ii = length(AllNan_vals):-1:1
        rmNan_xlabels(AllNan_vals(ii)) = [];
    end
end