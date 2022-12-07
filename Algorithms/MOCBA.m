function MOCBA(Global)

    %% Parameters, function and scenario
    test_function = func2str(Global.problem);
    if Global.M == 2
        if strcmp(test_function,'DTLZ7')
            filename = 'DTLZ7_M2_D2_50_100.mat';  
        elseif strcmp(test_function,'ZDT1')
            filename = 'ZDT1_M2_D10_80_100.mat';  
        else
            filename = 'WFG4_M2_D5_20_100.mat';
            %filename = 'WFG4_M2_D5_2.mat';
        end
    else
        if strcmp(test_function,'DTLZ7')
            filename = 'DTLZ7_M3_D5_150_200.mat';
        else
            filename = 'WFG4_M3_D5_100_200.mat';
        end
    end

    load(filename, 'Population'); 
    ref = Global.PF;
    tstart = tic; 
    %Noise level: 1: low, 2: high, 3:large range
    level = 2;
    %Noise case: 1: best case, 2: worst case
    caso = 1;
    %Replication budget
    B = 5;
    
    %% Initialization 
    size_set = size(Population.decs,1);
    %Pop_init = Population;
   
    %Heterogeneous noise
    [constants] = heter_noise_SK(Global, Population, level, caso,B);
    
    %True PS
    design_nondom = NDSort(Population.objs,1);
    design_nondom_size = sum(design_nondom(:) == 1);
    print = ['Size of design space 1st non-dominated set = ',num2str(design_nondom_size)];
    disp(print);
        
    %Initialize objective and variance matrices
    Cell_Obj_rep = cell(size_set,1); %Save each replication during search
    Mat_Obj = zeros(size_set,Global.M);
    Mat_Var = zeros(size_set,Global.M);
    
    %True HV
    true_pareto_set = Population(design_nondom == 1);          
    true_PF = true_pareto_set.objs;
    true_hv = HV(true_PF,ref);
     
    %Simulate (expensive) objectives on all design points
    for i = 1 : size_set
        det_obj = Population(1, i).obj;
        Mat_Obj_rep_i = zeros(B,Global.M);
        for k = 1 : Global.M
            a = constants(k,1);
            t = constants(k,2);
            for j = 1 : B
                Mat_Obj_rep_i(j,k) = det_obj(k)+normrnd(0,(a*det_obj(k)+a*t));
            end
            Mat_Obj(i,k) = mean(Mat_Obj_rep_i(:,k)); %Save mean of the response
            Mat_Var(i,k) = var(Mat_Obj_rep_i(:,k))/B; %Save variance
        end
        %disp('mat i'); disp(Mat_Obj_rep_i);pause;
        Cell_Obj_rep{i} = Mat_Obj_rep_i;
    end
    %
    %disp('cell'); disp(Cell_Obj_rep);pause;
    %disp('obj'); disp(Mat_Obj);pause;
    %disp('var'); disp(Mat_Var);pause;
    
    %% Main loop
    Global.evaluated = 0;
    iter = Global.evaluation+1;
    iteration = 1;
    pop_size_new = size_set;
    while iter > 0
        
              
        %% Save data 
        iter = iter-1;
        %Retreive data with sample means and variances
        [type1_size, type1_size_mocba, type2_size, type2_size_mocba, obs_hv, obs_hv_mocba, ...
        obs_igd, obs_igd_mocba, prec, rec, f1, prec_m, rec_m, f1_m] = metrics_val( Global, Population, Mat_Obj, Mat_Var);

        siz(iteration) = pop_size_new; %#ok<AGROW>
        true_HV(iteration) = true_hv; %#ok<AGROW>

        %Results using sample means
        ET1(iteration) = type1_size; %#ok<AGROW>
        ET2(iteration) = type2_size; %#ok<AGROW>
        ET1_m(iteration) = type1_size_mocba; %#ok<AGROW>
        ET2_m(iteration) = type2_size_mocba; %#ok<AGROW>
        
        obs_IGD(iteration) = obs_igd; %#ok<AGROW>
        obs_IGD_m(iteration) = obs_igd_mocba; %#ok<AGROW>
        obs_HV(iteration) = obs_hv; %#ok<AGROW>
        obs_HV_m(iteration) = obs_hv_mocba; %#ok<AGROW>
        mean_HV(iteration) = mean(obs_HV); %#ok<AGROW>
        mean_HV_m(iteration) = mean(obs_HV_m); %#ok<AGROW>
        gap_hv(iteration) = abs(true_hv-obs_hv); %#ok<AGROW> 

        PREC(iteration) = prec; %#ok<AGROW>
        REC(iteration) = rec; %#ok<AGROW>
        F1(iteration) = f1; %#ok<AGROW>

        PREC_m(iteration) = prec_m; %#ok<AGROW>
        REC_m(iteration) = rec_m; %#ok<AGROW>
        F1_m(iteration) = f1_m; %#ok<AGROW>

        %Performance metrics
        prt = ['Metrics and plots for obs means at iteration ', num2str(Global.evaluation-iter)];
        disp(prt);
        metrics(Global, Population, Mat_Obj, Mat_Var);
        %plots(Global,Population, test_function, Mat_Obj, Mat_Var, 1); 
        fprintf('\n');
        iterations = ['Number of iterations performed = ',num2str(Global.evaluation-iter)];
        disp(iterations);
        %pause;

        %Check termination condition
        if ((type1_size == 0) && (type2_size == 0)) || (iteration == 30)
            elapsed_time = toc(tstart);
            save('saved_file','ET1','ET1_m','ET2','ET2_m','PREC','PREC_m',...
                'REC','REC_m','F1','F1_m','obs_HV','obs_HV_m','true_HV','gap_hv',...
                'mean_HV','mean_HV_m','obs_IGD','obs_IGD_m','elapsed_time','Cell_Obj_rep',...
                'Mat_Obj','Mat_Var','siz','iteration');
            break
        end
        iteration = iteration + 1;

        %Simulate extra reps according to mocba
        pop_new = Population;
        size_set = size(Population,2);
        total_budget = B*size_set;
        alphas = mocba(Mat_Obj',Mat_Var');
        %[dom, nondom] = mocba_separate_designs(Mat_Obj',Mat_Var');
        budget = round(alphas*total_budget);
        if all(budget == 0)
            [~,max_ind] = max(alphas);
            budget(max_ind) = 1;
        end
        for i = 1:size_set
            rep = budget(i);
            pos = i;
            det_obj = pop_new(1, pos).obj;
            obj_temp = zeros(rep,Global.M);
            b_i = size(Cell_Obj_rep{pos},1);
            Mat_Obj_i = [Cell_Obj_rep{pos};obj_temp];    
            %disp('first mat'); disp(Cell_Obj_rep{pos}); pause;
            for k = 1 : Global.M
                a = constants(k,1);
                t = constants(k,2); 
                r1 = b_i+1;
                r2 = b_i+rep;
                for j = r1 : r2
                    Mat_Obj_i(j,k) = det_obj(k)+normrnd(0,(a*det_obj(k)+a*t));
                end
                Mat_Obj(pos,k) = mean(Mat_Obj_i(:,k)); %Save new mean response
                Mat_Var(pos,k) = var(Mat_Obj_i(:,k))/r2; %Save new variance
            end
            %disp('half mat'); disp(Mat_Obj_i); pause;

            Cell_Obj_rep{pos} = Mat_Obj_i;
            %disp('sec mat'); disp(Cell_Obj_rep); pause;
        end
    end
    close all
end
