function metrics( Global, Population, Mat_Obj, Mat_Var)

%% Retrieve data
%[~, ~, ~, PF_obs, ~, ~, ~, ~, ~, ~, type1_size, ~, type2_size, ~, NonDominated_obs_size, PS_sampled_size, PS_s, ...
%   prec, rec, f1 ] = sets(Global, Population, Mat_Obj);

[~, ~, ~, PF_obs, ~, ~, ~, ~, ~, ~, PF_mocba, ~, ~, ~, type1_size, type1_size_mocba, ~, ~, type2_size, type2_size_mocba, ...
    ~, NonDominated_obs_size, NonDominated_obs_size_mocba, PS_sampled_size, PS_s, prec, rec, f1, prec_m, rec_m, f1_m ] = ...
    sets(Global, Population, Mat_Obj, Mat_Var);


    %% Quality indicators
     %Hypervolume of observed front
     obs_hv = HV(PF_obs,Global.PF);
     %obs_hv2 = stk_dominatedhv(PF_obs,ref); 
     obs_hv_mocba = HV(PF_mocba,Global.PF);

     %IGD of observed front
     obs_igd = IGD(PF_obs,Global.PF);
     obs_igd_mocba = IGD(PF_mocba,Global.PF);

     %Print indicators
     hv_print = ['HV = ', num2str(obs_hv)];
     disp(hv_print);
     hv_print2 = ['HV_mocba = ', num2str(obs_hv_mocba)];
     disp(hv_print2);
     igd_print = ['IGD = ', num2str(obs_igd)];
     disp(igd_print);
     igd_print_mocba = ['IGD_mocba = ', num2str(obs_igd_mocba)];
     disp(igd_print_mocba);
     %pause;
    %}
    
    test = ['Number of points in observed PF = ',num2str(NonDominated_obs_size)];
    disp(test); 
    test1 = ['Number of points in observed PF_mocba = ',num2str(NonDominated_obs_size_mocba)];
    disp(test1); 
    
    print1 = ['Number of true Pareto points sampled = ',num2str(PS_sampled_size)];
    disp(print1);
    
    result1 = ['PS_s = ',num2str(PS_s) ,'%'];
    disp(result1);
    
    prt_t1 = ['Number of MCE error points = ',num2str(type1_size)];
    disp(prt_t1);
    prt_t11 = ['Number of MCE error points (mocba) = ',num2str(type1_size_mocba)];
    disp(prt_t11);
    
    prt_t2 = ['Number of MCI error points = ',num2str(type2_size)];
    disp(prt_t2);
    prt_t22 = ['Number of MCI error points (mocba) = ',num2str(type2_size_mocba)];
    disp(prt_t22);
    
    metric2 = ['precision = ', num2str(prec)];
    disp(metric2);
    
    metric3 = ['recall = ', num2str(rec)];
    disp(metric3);
    
    metric4 = ['f1 = ', num2str(f1)];
    disp(metric4);

    metric22 = ['precision (mocba) = ', num2str(prec_m)];
    disp(metric22);
    
    metric33 = ['recall (mocba) = ', num2str(rec_m)];
    disp(metric33);
    
    metric44 = ['f1 (mocba) = ', num2str(f1_m)];
    disp(metric44);
end

