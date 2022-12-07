function SKMORS_sm(Global)

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
    Pop_init = Population;
   
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

        %Print tables
        %{
        true_obj = Population.objs;
        for j = 1:Global.M
            prt = ['Objective ',num2str(j)];
            disp(prt);
            prt1 = ['  True Obj.  ','  Mean  ','      Var    ',' Prediction','    MSE   '];
            disp(prt1);
            pred = [true_obj(:,j),Mat_Obj(:,j),Mat_Var(:,j),y(:,j),mse(:,j)];
            disp(pred); 
        end 
        pause;
        %}


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
        if ((type1_size == 0) && (type2_size == 0)) || ((type1_size_sk == 0) && (type2_size_sk == 0)) || (iteration == 30)
            elapsed_time = toc(tstart);
            save('saved_file','ET1','ET1_sk','ET1_m','ET1_m_sk','ET2','ET2_sk','ET2_m','ET2_m_sk','PREC','PREC_sk','PREC_m','PREC_m_sk',...
                'REC','REC_sk','REC_m','REC_m_sk','F1','F1_sk','F1_m','F1_m_sk','obs_HV','obs_HV_sk','obs_HV_m','obs_HV_m_sk','true_HV','gap_hv','gap_hv_sk',...
                'mean_HV','mean_HV_sk','mean_HV_m','mean_HV_m_sk','obs_IGD','obs_IGD_sk','obs_IGD_m','obs_IGD_m_sk','elapsed_time','Cell_Obj_rep',...
                'Mat_Obj','Mat_Var','y','mse','siz','iteration');
            break
        end
        iteration = iteration + 1;


        %% Screening heuristic
        %Temporary sets
        pop_temp = Population;
        size_temp = size(Population,2);
        predictions = y;
        mses = mse;
        means = Mat_Obj;
        vars = Mat_Var;

         %Observed means 
         nondom_means = NDSort(means,1);
         nondom_means_pop = pop_temp(nondom_means == 1);
         dom_means_pop = pop_temp(nondom_means ~= 1);

          %Observed predictions 
         nondom_pred = NDSort(predictions,1);
         nondom_pred_pop = pop_temp(nondom_pred == 1);
         dom_pred_pop = pop_temp(nondom_pred ~= 1);

         %PF observed 
         PF_obs = zeros(length(nondom_means_pop),Global.M+1);
         nonPF_obs = zeros(length(dom_means_pop),Global.M+1);
         j = 1;
         k=1;
         for i = 1:size(means,1)
             if nondom_means(1,i) == 1
                 PF_obs(j,1) = i; %index
                 PF_obs(j,2:end) = means(i,:);
                 j = j+1;
             else
                 nonPF_obs(k,1) = i; %index
                 nonPF_obs(k,2:end) = means(i,:);
                 k = k+1;
             end
         end

         %disp('PF_obs'); disp(PF_obs);
         %disp('nonPF_obs'); disp(nonPF_obs);

         %PF predicted 
         PF_pred = zeros(length(nondom_pred_pop),Global.M+1);
         nonPF_pred = zeros(length(dom_pred_pop),Global.M+1);
         j = 1;
         k=1;
         for i = 1:size(predictions,1)
             if nondom_pred(1,i) == 1
                 PF_pred(j,1) = i; %index
                 PF_pred(j,2:end) = predictions(i,:);
                 j = j+1;
             else
                 nonPF_pred(k,1) = i; %index
                 nonPF_pred(k,2:end) = predictions(i,:);
                 k = k+1;
             end
         end

         %disp('PF_pred'); disp(PF_pred); 
         %disp('nonPF_pred'); disp(nonPF_pred); pause;

         %UCB for PF_obs points
         ucb_PF_obs = PF_obs;
         for i = 1:size(ucb_PF_obs,1)
             for j = 1:Global.M
                 k = fix(PF_obs(i,1));
                 %ucb_PF_obs(i,j+1) = PF_obs(i,j+1)+omega*sqrt(vars(k,j));
                 %
                 if iteration <= 5
                     ucb_PF_obs(i,j+1) = PF_obs(i,j+1)+3*sqrt(vars(k,j));
                 elseif (iteration > 5) && (iteration <= 20)
                     ucb_PF_obs(i,j+1) = PF_obs(i,j+1)+2*sqrt(vars(k,j));
                 else
                     ucb_PF_obs(i,j+1) = PF_obs(i,j+1)+1*sqrt(vars(k,j));
                 end
                 %}
             end
         end


         %UCB for PF_pred points
         ucb_PF_pred = PF_pred;
         for i = 1:size(ucb_PF_pred,1)
             for j = 1:Global.M
                 k = fix(PF_pred(i,1));
                 %ucb_PF_pred(i,j+1) = PF_pred(i,j+1)+omega*sqrt(mses(k,j));
                 %
                 if iteration <= 5
                     ucb_PF_pred(i,j+1) = PF_pred(i,j+1)+3*sqrt(mses(k,j));
                 elseif (iteration > 5) && (iteration <= 20)
                     ucb_PF_pred(i,j+1) = PF_pred(i,j+1)+2*sqrt(mses(k,j));    
                 else
                     ucb_PF_pred(i,j+1) = PF_pred(i,j+1)+1*sqrt(mses(k,j));
                 end
                 %}
             end
         end

         %LCB for nonPF_obs points
         lcb_nonPF_obs = nonPF_obs;
         for i = 1:size(lcb_nonPF_obs,1)
             for j = 1:Global.M
                 k = fix(nonPF_obs(i,1));
                 %lcb_nonPF_obs(i,j+1) = nonPF_obs(i,j+1)-omega*sqrt(vars(k,j));
                 %
                 if iteration <= 5
                     lcb_nonPF_obs(i,j+1) = nonPF_obs(i,j+1)-3*sqrt(vars(k,j));
                 elseif (iteration > 5) && (iteration <= 20)
                     lcb_nonPF_obs(i,j+1) = nonPF_obs(i,j+1)-2*sqrt(vars(k,j));
                 else
                     lcb_nonPF_obs(i,j+1) = nonPF_obs(i,j+1)-1*sqrt(vars(k,j));
                 end
                 %}
             end
         end

         %LCB for nonPF_pred points
         lcb_nonPF_pred = nonPF_pred;
         for i = 1:size(lcb_nonPF_pred,1)
             for j = 1:Global.M
                 k = fix(nonPF_pred(i,1));
                 %lcb_nonPF_pred(i,j+1) = nonPF_pred(i,j+1)-omega*sqrt(mses(k,j));
                 %
                 if iteration <= 5
                     lcb_nonPF_pred(i,j+1) = nonPF_pred(i,j+1)-3*sqrt(mses(k,j));
                 elseif (iteration > 5) && (iteration <= 20)
                     lcb_nonPF_pred(i,j+1) = nonPF_pred(i,j+1)-2*sqrt(mses(k,j));
                 else
                    lcb_nonPF_pred(i,j+1) = nonPF_pred(i,j+1)-1*sqrt(mses(k,j));
                 end
                 %}
             end
         end

        %{
        new_means = [ucb_PF_obs;lcb_nonPF_obs];
        new_means = sortrows(new_means);

        new_preds = [ucb_PF_pred;lcb_nonPF_pred];
        new_preds = sortrows(new_preds);

        disp('originals'); 
        prt = [Mat_Obj,Mat_Var];
        disp(prt);
        disp('new'); disp(new_means); pause;

        disp('originals'); 
        prt = [y,mse];
        disp(prt);
        disp('new'); disp(new_preds); pause;
        %}

        
        %% Method 1: Distance to worst UB of current PF
        %
        %Screen out bad points
        bad_points_mean = zeros(1,1);
        bad_points_pred = zeros(1,1);
        k=1;
        for i = 1:size(lcb_nonPF_obs,1)
            for j = 1:Global.M
                if max(ucb_PF_obs(:,j+1)) < lcb_nonPF_obs(i,j+1) 
                    bad_points_mean(k,1) = lcb_nonPF_obs(i,1); 
                    k = k+1;
                    break;
                end
            end
        end

        k = 1;
        for i = 1:size(lcb_nonPF_pred,1)
            for j = 1:Global.M
                if max(ucb_PF_pred(:,j+1)) < lcb_nonPF_pred(i,j+1)
                    bad_points_pred(k,1) = lcb_nonPF_pred(i,1); 
                    k = k+1;
                    break;
                end
            end
        end
        
        %% Method 2: Closest UB of any member in PF
        %{
        %Screen out bad points
        bad_points_mean = zeros(1,1);
        bad_points_pred = zeros(1,1);
        %Distances to the closest point in current front
        %index = (1:size(nonPF_obs,1))';
        dist_obs = pdist2(PF_obs,nonPF_obs);
        [~, ind_obs] = min(dist_obs,[],1);
        %disp('size PF obs'); disp(size(PF_obs,1));
        %disp('obs'); disp(ind_obs); pause;

        dist_pred = pdist2(PF_pred,nonPF_pred);
        [~, ind_pred] = min(dist_pred,[],1);
        %disp('size PF pred'); disp(size(PF_pred,1));
        %disp('pred'); disp(ind_pred);pause;

        k=1;
        for i = 1:size(lcb_nonPF_obs,1)
            for j = 1:Global.M
                if ucb_PF_obs(ind_obs(1,i),j+1) < lcb_nonPF_obs(i,j+1) 
                    bad_points_mean(k,1) = lcb_nonPF_obs(i,1); 
                    k = k+1;
                    break;
                end
            end
        end

        k=1;
        for i = 1:size(lcb_nonPF_pred,1)
            for j = 1:Global.M
                if ucb_PF_pred(ind_pred(1,i),j+1) < lcb_nonPF_pred(i,j+1)
                    bad_points_pred(k,1) = lcb_nonPF_pred(i,1); 
                    k = k+1;
                    break;
                end
            end
        end
        %}
        
        
        %% Screen out bad points
        bad_points = intersect(bad_points_mean,bad_points_pred);
        good_points = setdiff(1:size_temp,(bad_points)'); 
        %disp(good_points); pause;
        %disp(bad_points); pause;

        %Keep only points with good bounds
        pop_new = pop_temp(good_points);
        means_new = means((good_points)',:);
        vars_new = vars((good_points)',:);
        preds_new = predictions((good_points)',:);
        mses_new = mses((good_points)',:); %#ok<NASGU>

        pop_size_new = size(pop_new,2);
        disp(pop_size_new);

        %{
        plots(Global,pop_new, test_function, means_new);
        metrics(Global, pop_new, means_new);
        %pause;
        plots(Global,pop_new, test_function, preds_new);
        metrics(Global, pop_new, preds_new);
        %pause;
        %}
        %}
              
        
        %All points
        if mod(iteration,5) == 0
            pop_new = Population;
            means_new = means;
            vars_new = vars;
            preds_new = predictions;
            mses_new = mses; %#ok<NASGU>
            pop_size_new = size(pop_new,2);
            disp(pop_size_new);
            good_points = (1:size_temp);
        end
        
        
        %% MORS procedure

        %Compute confidence bounds of predictions (hat(y))
         cb_pred = preds_new;
         for i = 1:pop_size_new
             for j = 1:Global.M
                 if preds_new(i,j) <= means_new(i,j)
                     cb_pred(i,j) = preds_new(i,j)-sqrt(vars_new(i,j));
                 else
                     cb_pred(i,j) = preds_new(i,j)+sqrt(vars_new(i,j));
                 end
             end
         end

         %Compute distances and normalize
         d = diag(pdist2(cb_pred,means_new));
         d_norm = d;
         d_min = min(d);
         d_max = max(d);
         for i =1:pop_size_new
             d_norm(i) = (d(i)-d_min)/(d_max-d_min);
         end

         %disp('dist norm'); disp(d_norm); pause;

         %% EHV Change

         %Compute HV of current front 
         obs_HV_tot = HV(means_new,ref);  
         %pred_HV_tot = HV(preds_new,ref); 
         hv_diff = zeros(length(cb_pred),1);

         %
         parfor i = 1:pop_size_new
             temp_means = means_new;
             temp_means(i,:) = preds_new(i,:);
             %Compute new HV with prediction value
             obs_HV = HV(temp_means,ref);
             hv_diff(i) = abs(obs_HV-obs_HV_tot);
         end
         %

         %OR

         %{
         parfor i = 1:pop_size_new
             temp_preds = preds_new;
             temp_preds(i,:) = means_new(i,:);
             %Compute new HV with prediction value
             pred_HV = HV(temp_preds,ref);
             hv_diff(i) = abs(pred_HV-pred_HV_tot);
         end
         %}

         %disp('hv diff'); disp(hv_diff); pause;

         %% Normalize HV differences
         hv_diff_norm = hv_diff;
         hv_min = min(hv_diff);
         hv_max = max(hv_diff);
         for i =1:size(hv_diff,1)
             hv_diff_norm(i) = (hv_diff(i)-hv_min)/(hv_max-hv_min);
         end

         %disp('hv diff norm'); disp(hv_diff_norm); pause;  
      
         
         %% Rank and select according to a non-dominated sort
         pf = [hv_diff_norm,d_norm];
         %disp(pf); pause;
         max_pf = NDSort(-pf,1);
         %disp(max_pf); pause;
         top_points = find(max_pf == 1);
         %disp(top_points);pause;
         %plot_mors(2,-pf); pause;
         
         
         %Retrieve points
         best_points = good_points(top_points); %#ok<FNDSB> 
         best_points_size = length(best_points);  
            

        %% Simulate extra reps on best points
        Budget = B*size(Pop_init,2);
        for i = 1:best_points_size
            %rep = round(Budget*best_points_rep(i));
            rep = round(Budget/best_points_size);
            pos = best_points(i);
            det_obj = Population(1, pos).obj;
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
        
        
        %Plots visualization
        %{
        figure(10000);
        plots_vis(Global,Population, test_function, y);
        pause(0.01);
        %}
        %{
        figure(10000);
        plots_vis(Global,pop_new, test_function, preds_new);
        pause(0.01);
        %}
        %
        
        %{
        figure(10000);
        plot(obs_HV);
        hold on
        plot(mean_HV);
        hold on
        plot(true_HV);
        hold off
        pause(0.01);
        %}
    end

    %close all; 

end
