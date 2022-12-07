function MOCBASK(Global)

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
    level = 3;
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
        %% Fit a SK metamodel to each objective
        PDec = Population.decs;
        %Cell array of SK model objects
        SK = cell(1,Global.M);
        for s = 1:Global.M
            SK{s} = SKfit_new(PDec,Mat_Obj(:,s),ones(size_set,1),Mat_Var(:,s),2,3);
        end

        y = zeros(size_set, Global.M);
        mse = zeros(size_set, Global.M);
        %Predict values of design points
        for i = 1:size_set    
            for m = 1:Global.M
                [y(i,m), mse(i,m)] = SKpredict_new(SK{m},Population(1,i).dec,1);
            end
        end
              
        %% Save data 
        iter = iter-1;
        %Retreive data with sample means and variances
        [type1_size, type1_size_mocba, type2_size, type2_size_mocba, obs_hv, obs_hv_mocba, ...
        obs_igd, obs_igd_mocba, prec, rec, f1, prec_m, rec_m, f1_m] = metrics_val( Global, Population, Mat_Obj, Mat_Var);

        %Retreive data with sk predictions and mse
        [type1_size_sk, type1_size_mocba_sk, type2_size_sk, type2_size_mocba_sk, obs_hv_sk, obs_hv_mocba_sk, ...
        obs_igd_sk, obs_igd_mocba_sk, prec_sk, rec_sk, f1_sk, prec_m_sk, rec_m_sk, f1_m_sk] = metrics_val( Global, Population, y, mse);

        %[type1_size, type2_size,~,~,obs_hv,obs_igd, prec, rec, f1] = metrics_val(Global, Population, y, mse);

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

        %Results using predicted means
        ET1_sk(iteration) = type1_size_sk; %#ok<AGROW>
        ET2_sk(iteration) = type2_size_sk; %#ok<AGROW>
        ET1_m_sk(iteration) = type1_size_mocba_sk; %#ok<AGROW>
        ET2_m_sk(iteration) = type2_size_mocba_sk; %#ok<AGROW>
        
        obs_IGD_sk(iteration) = obs_igd_sk; %#ok<AGROW>
        obs_IGD_m_sk(iteration) = obs_igd_mocba_sk; %#ok<AGROW>
        obs_HV_sk(iteration) = obs_hv_sk; %#ok<AGROW>
        obs_HV_m_sk(iteration) = obs_hv_mocba_sk; %#ok<AGROW>
        mean_HV_sk(iteration) = mean(obs_HV_sk); %#ok<AGROW>
        mean_HV_m_sk(iteration) = mean(obs_HV_m_sk); %#ok<AGROW>
        gap_hv_sk(iteration) = abs(true_hv-obs_hv_sk); %#ok<AGROW> 

        PREC_sk(iteration) = prec_sk; %#ok<AGROW>
        REC_sk(iteration) = rec_sk; %#ok<AGROW>
        F1_sk(iteration) = f1_sk; %#ok<AGROW>

        PREC_m_sk(iteration) = prec_m_sk; %#ok<AGROW>
        REC_m_sk(iteration) = rec_m_sk; %#ok<AGROW>
        F1_m_sk(iteration) = f1_m_sk; %#ok<AGROW>

        %Performance metrics
        prt = ['Metrics and plots for obs means at iteration ', num2str(Global.evaluation-iter)];
        disp(prt);
        metrics(Global, Population, Mat_Obj, Mat_Var);
        %plots(Global,Population, test_function, Mat_Obj, Mat_Var, 1); 
        fprintf('\n');
        prt = ['Metrics and plots for pred means at iteration ', num2str(Global.evaluation-iter)];
        disp(prt);
        metrics(Global, Population, y, mse);
        %plots(Global,Population, test_function, y, mse, 1); 
        iterations = ['Number of iterations performed = ',num2str(Global.evaluation-iter)];
        disp(iterations);
        %pause;

        %Check termination condition
        %if ((type1_size == 0) && (type2_size == 0)) || ((type1_size_sk == 0) && (type2_size_sk == 0)) || (iteration == 30)
        if  ((type1_size_sk == 0) && (type2_size_sk == 0)) || (iteration == 30)
            elapsed_time = toc(tstart);
            save('saved_file','ET1','ET1_sk','ET1_m','ET1_m_sk','ET2','ET2_sk','ET2_m','ET2_m_sk','PREC','PREC_sk','PREC_m','PREC_m_sk',...
                'REC','REC_sk','REC_m','REC_m_sk','F1','F1_sk','F1_m','F1_m_sk','obs_HV','obs_HV_sk','obs_HV_m','obs_HV_m_sk','true_HV','gap_hv','gap_hv_sk',...
                'mean_HV','mean_HV_sk','mean_HV_m','mean_HV_m_sk','obs_IGD','obs_IGD_sk','obs_IGD_m','obs_IGD_m_sk','elapsed_time','Cell_Obj_rep',...
                'Mat_Obj','Mat_Var','y','mse','siz','iteration');
            break
        end
        iteration = iteration + 1;

        %Simulate extra reps according to mocba
        pop_new = Population;
        size_set = size(Population,2);
        total_budget = B*size_set;
        alphas = mocba(Mat_Obj',Mat_Var');
        %alphas = mocba(y',mse');
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
