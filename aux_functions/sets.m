%% Copyright 2020 Sebastian Rojas Gonzalez <srojas33@gmail.com><sebastian.rojasgonzalez@ugent.be>
%{
Source code for the paper: 
%--------------------------------------------------------------------------
Gonzalez, S. R., Jalali, H., & Van Nieuwenhuyse, I. (2020). A multiobjective 
stochastic simulation optimization algorithm. European Journal of Operational 
Research, 284(1), 212-226.
%--------------------------------------------------------------------------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the above copyright notice is
retained.
%}


function [pop_set, true_pareto_set, true_PF, PF_obs, nonPF_obs, PF_obs_true, PS_obs, PS_mocba, PF_pop, nonPF_pop, PF_mocba, nonPF_mocba, ...
    type1_obs_obj, type1_obs_obj_mocba, type1_size, type1_size_mocba, type2_obs_obj, type2_obs_obj_mocba, type2_size, type2_size_mocba, ...
    design_nondom_size, NonDominated_obs_size, NonDominated_obs_size_mocba, PS_sampled_size, PS_s, prec, rec, f1, prec_m, rec_m, f1_m ] = ...
    sets(Global, Population, Mat_Obj, Mat_Var)

    %% Design space      
    %design_space = Population.decs;
    %candid_set = size(design_space,1);
    %pause;

    %True PS
    design_nondom = NDSort(Population.objs,1);
    design_nondom_size = sum(design_nondom(:) == 1);
    true_pareto_set = Population(design_nondom == 1);

    %True PF             
    true_PF = true_pareto_set.objs;
     
    %% Observed sets
    pop_set = size(Population.decs,1);
    
    %Observed means and points
     NonDominated_obs = NDSort(Mat_Obj,1);
     NonDominated_obs_size = sum(NonDominated_obs(:) == 1);
     
     %PF observed 
     PF_obs = zeros(NonDominated_obs_size,Global.M);
     nonPF_obs = zeros(pop_set-NonDominated_obs_size,Global.M);
     j = 1;
     k=1;
     for i = 1:size(Mat_Obj,1)
         if NonDominated_obs(1,i) == 1
             PF_obs(j,:) = Mat_Obj(i,:);
             j = j+1;
         else
             nonPF_obs(k,:) = Mat_Obj(i,:);
             k = k+1;
         end
     end
     

    %Non-dominated set of population (not entirely identified)
    PS_sampled = NDSort(Population.objs,1);
    PS_sampled_size = sum(PS_sampled(:) == 1);


    %Intersect PS sampled with PS of design space          
    int_points1 = intersect(Population.decs, true_pareto_set.decs, 'rows');
    int_points1_size = size(int_points1,1);
    PS_s = int_points1_size/design_nondom_size*100;


     % True objectives of PFobs
     PF_obs_true = zeros(NonDominated_obs_size,Global.M);
     j = 1;
     for i = 1:size(Mat_Obj,1)
         if NonDominated_obs(1,i) == 1
             PF_obs_true(j,:) = Population(:,i).obj;
             j = j+1;
         end
     end

     %PS observed
     PS_obs = zeros(NonDominated_obs_size,Global.D);
     j = 1;
     for i = 1:size(Mat_Obj,1)
         if NonDominated_obs(1,i) == 1
             PS_obs(j,:) = Population(1,i).dec;
             j = j+1;
         end
     end

    %True values of PF sampled (not entirely identified) 
    PF_pop = zeros(PS_sampled_size,Global.M);
    nonPF_pop = zeros(pop_set-PS_sampled_size,Global.M);
     j = 1;
     k = 1;
     for i = 1:size(Mat_Obj,1)
         if PS_sampled(1,i) == 1
             PF_pop(j,:) = Population(1,i).obj;
             j = j+1;
         else
             nonPF_pop(k,:) = Population(1,i).obj;
             k = k+1;
         end
     end

     %% Sets according to MOCBA
     [mocba_dom_ind, mocba_nondom_ind] = mocba_separate_designs(Mat_Obj',Mat_Var');
     mocba_dom_ind = mocba_dom_ind + 1;
     mocba_nondom_ind = mocba_nondom_ind + 1;
     %disp(mocba_dom_ind); pause;
     %PFmocba
     PF_mocba = Mat_Obj(mocba_nondom_ind',:);
     %disp(PF_mocba);pause
     NonDominated_obs_size_mocba = size(PF_mocba,2);
     nonPF_mocba = Mat_Obj(mocba_dom_ind',:);
     PS_mocba = Population(mocba_nondom_ind').decs;

     
     
     %% Misclassification errors 
   
     %% Misclassification by exclusion (ET1)
     %First, remove set found from true PS
     int_points_t1 = true_pareto_set.decs;
     int_points_t1(ismember(int_points_t1,PS_obs,'rows'),:)=[];
     %Second, remove set found from all points sampled
     int_points_t1_1 = Population.decs;
     int_points_t1_1(ismember(int_points_t1_1,PS_obs,'rows'),:)=[];
     %Third, intersect both sets for points error Type 1
     type1 = intersect(int_points_t1, int_points_t1_1, 'rows');
     type1_size = size(type1,1);
     %MCE errors
     [~,type1_ind] = intersect(Population.decs,type1,'rows');
     type1_obs_obj = Mat_Obj(type1_ind,:);

     %According to MOCBA
     %First, remove set found from true PS
     int_points_t1_mocba = true_pareto_set.decs;
     int_points_t1_mocba(ismember(int_points_t1_mocba,PS_mocba,'rows'),:)=[];
     %Second, remove set found from all points sampled
     int_points_t1_1_mocba = Population.decs;
     int_points_t1_1_mocba(ismember(int_points_t1_1_mocba,PS_mocba,'rows'),:)=[];
     %Third, intersect both sets for points error Type 1
     type1_mocba = intersect(int_points_t1_mocba, int_points_t1_1_mocba, 'rows');
     type1_size_mocba = size(type1_mocba,1);
     %MCE errors
     [~,type1_ind_mocba] = intersect(Population.decs,type1_mocba,'rows');
     type1_obs_obj_mocba = Mat_Obj(type1_ind_mocba,:);

     %% Misclassification by inclusion (ET2)
     %First remove true points from set found
     type2 = PS_obs;
     type2(ismember(type2,true_pareto_set.decs,'rows'),:)=[];
     %MCI errors
     type2_size = size(type2,1);
     [~,type2_ind] = intersect(Population.decs,type2,'rows');
     type2_obs_obj = Mat_Obj(type2_ind,:);

     %According to MOCBA
     %First remove true points from set found
     type2_mocba = PS_mocba;
     type2_mocba(ismember(type2_mocba,true_pareto_set.decs,'rows'),:)=[];
     %MCI errors
     type2_size_mocba = size(type2_mocba,1);
     [~,type2_ind_mocba] = intersect(Population.decs,type2_mocba,'rows');
     type2_obs_obj_mocba = Mat_Obj(type2_ind_mocba,:);



     %% Error metrics

     %Intersect PS sampled with PS identified 
     %int_points2 = intersect(int_points1, PS_obs, 'rows');
     %int_points2_size = size(int_points2,1); 
     %APS = (int_points1_size+int_points2_size)/(2*design_nondom_size+type2_size)*100;
    
    %APS = 1-(type1_size+type2_size)/pop_set;
    prec = design_nondom_size/(design_nondom_size+type2_size);
    rec = design_nondom_size/(design_nondom_size+type1_size);
    f1 = 2*prec*rec/(prec+rec);

    prec_m = design_nondom_size/(design_nondom_size+type2_size_mocba);
    rec_m = design_nondom_size/(design_nondom_size+type1_size_mocba);
    f1_m = 2*prec_m*rec_m/(prec_m+rec_m);

end

