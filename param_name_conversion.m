% make a file to convert param names to figure appropriate
function n = param_name_conversion(param_name)
    %display(isstring(param_name))
    if strcmp(param_name, 'beta_exo_PTHg ')
        n = '\beta_{exo}^{PTHg}';
    elseif strcmp(param_name, 'k_prod_PTHg ')
        n = 'k_{prod}^{PTHg}';
    elseif strcmp(param_name, 'gamma_deg_PTHp ')
        n = '\gamma_{deg}^{PTHp}';
    elseif strcmp(param_name, 'gamma_prod_D3 ')
        n = '\gamma_{prod}^{D3}';
    elseif strcmp(param_name, 'D3_inact ')
        n = 'D_{3}^{inact}';
    elseif strcmp(param_name, 'k_PTHg_deg ')
        n = 'k_{deg}^{PTHg}';   
    elseif strcmp(param_name, 'k_PTHp_deg ')
        n = 'k_{deg}^{PTHp}';    
    elseif strcmp(param_name, 'D3_act_peak ')
        n = 'D_{3}^{Mg-peak}';
    elseif strcmp(param_name, 'D3_act_s ')
        n = 'D_{3}^{Mg-std}';
    elseif strcmp(param_name, 'k_min_25D_prod ' )
        n = 'k_{25D-prod}^{min}';  
    elseif strcmp(param_name, 'max_Mg_effect_25D ')
        n = 'V_{Mg-25D}^{max}';
    elseif strcmp(param_name, 'K_Mg_25D ')
        n = 'K_{Mg-25D}';
    elseif strcmp(param_name, 'K_D3 ')
        n = 'K_{D3}';
    elseif strcmp(param_name, 'nD3 ')
        n = 'n_{D3}';
    elseif strcmp(param_name, 'k_deg_D3 ')
        n = 'k_{deg}^{D3}';
    elseif strcmp(param_name, 'K_conv_PTH ')
        n = 'K_{conv}^{PTH}';
    elseif strcmp(param_name, 'k_min_conv ')
        n = 'k_{conv}^{min}';
    elseif strcmp(param_name, 'Lambda_PT0_Mg ')
        n = '\lambda_{PT0}^{Mg}';
    elseif strcmp(param_name, 'Lambda_PT0_Ca ')
        n = '\lambda_{PT0}^{Ca}';
    elseif strcmp(param_name, 'Lambda_TAL0_Mg ')
        n = '\lambda_{TAL0}^{Mg}';
    elseif strcmp(param_name, 'Lambda_TAL0_Ca ')
        n = '\lambda_{TAL0}^{Ca}';
    elseif strcmp(param_name, 'Lambda_DCT0_Mg ')
        n = '\lambda_{DCT0}^{Mg}';
    elseif strcmp(param_name, 'Lambda_DCT0_Ca ')
        n = '\lambda_{DCT0}^{Ca}';
    elseif strcmp(param_name, 'k_pf_Mg ')
        n = 'k_{pf}^{Mg}';
    elseif strcmp(param_name, 'k_pf_Ca ')
        n = 'k_{pf}^{Ca}';
    elseif strcmp(param_name, 'k_fp_Mg ')
        n = 'k_{fp}^{Mg}';
    elseif strcmp(param_name, 'k_fp_Ca ')
        n = 'k_{fp}^{Ca}';
    elseif strcmp(param_name, 'Gamma_res_min ')
        n = '\Gamma_{res}^{min}';
    elseif strcmp(param_name, 'delta_res_max ')
        n = '\delta_{res}^{max}';
    elseif strcmp(param_name, 'K_D3p_res ')
        n = 'K_{D3}^{res}';
    elseif strcmp(param_name, 'gamma_p ')
        n = '\gamma^p';
    elseif strcmp(param_name, 'Cm_half ')
        n = 'C_{m1/2}';
    elseif strcmp(param_name, 'Gamma_ac_Ca ')
        n = '\gamma_{ac}^{Ca}';
    elseif strcmp(param_name, 'gamma_conv_D3 ')
        n = '\gamma_{conv}^{D3}';
    elseif strcmp(param_name, 'delta_max_D3act ')
        n = '\delta_{max}^{D3}';
    elseif strcmp(param_name, 'ICa ')
        n = 'I_{Ca}';
    elseif strcmp(param_name, 'IMg ')
        n = 'I_{Mg}';
    elseif strcmp(param_name, 'K1_Mg_act ')
        n = 'K_{Mg-act}^1';
    elseif strcmp(param_name, 'gamma_Ca ')
        n = '\gamma_{Ca}';
    elseif strcmp(param_name, 'delta_Mg_act ')
        n = '\delta_{Mg-act}';
    elseif strcmp(param_name, 'Gamma_abs0 ')
        n = '\Gamma_{abs0}';
    elseif strcmp(param_name, 'delta_abs_D3 ')
        n = '\delta_{abs}^{D3}';
    else
        n = param_name;
    end
end