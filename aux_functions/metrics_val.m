function [type1_size, type1_size_mocba, type2_size, type2_size_mocba, obs_hv, obs_hv_mocba, ...
    obs_igd, obs_igd_mocba, prec, rec, f1, prec_m, rec_m, f1_m] = metrics_val( Global, Population, Mat_Obj, Mat_Var)


%% Retrieve data
%[size_set, ~, ~, PF_obs, ~, ~, ~, ~, ~, ~, type1_size, ~, type2_size, design_nondom_size, ~,~,~, prec, rec, f1 ] = ...
%    sets(Global, Population, Mat_Obj);

[~, ~, ~, PF_obs, ~, ~, ~, ~, ~, ~, PF_mocba, ~, ~, ~, type1_size, type1_size_mocba, ~, ~, type2_size, type2_size_mocba, ...
    ~, ~, ~, ~, ~, prec, rec, f1, prec_m, rec_m, f1_m ] = sets(Global, Population, Mat_Obj, Mat_Var);

obs_hv = HV(PF_obs,Global.PF);
obs_hv_mocba = HV(PF_mocba,Global.PF);

obs_igd = IGD(PF_obs,Global.PF);
obs_igd_mocba = IGD(PF_mocba,Global.PF);

end

