%function [size_set, true_pareto_set, true_PF, PF_obs, nonPF_obs, PF_obs_true, PS_obs, PF_pop, nonPF_pop, type1_obs_obj, type1_size, ...
%    type2_obs_obj, type2_size, design_nondom_size, NonDominated_obs_size, PS_sampled_size, PS_s, PS_ID, APS ] = sets_vis(Global, ...
%    Population, Mat_Obj)

function [PF_obs, nonPF_obs, type1_obs_obj, type1_size, type2_obs_obj, type2_size] = sets_vis(Global, Population, Mat_Obj)

    %%Design space      
    design_space = Population.decs;
    size_set = size(design_space,1);
    %pause;

%True PS
    %True PS
    design_nondom = NDSort(Population.objs,1);
    design_nondom_size = sum(design_nondom(:) == 1);
    true_pareto_set = Population(design_nondom == 1);

    %True PF             
     %true_PF = true_pareto_set.objs;
     
      %Observed means and points
     NonDominated_obs = NDSort(Mat_Obj,1);
     NonDominated_obs_size = sum(NonDominated_obs(:) == 1);
     
     

     %PF observed 
     PF_obs = zeros(NonDominated_obs_size,Global.M);
     nonPF_obs = zeros(size_set-NonDominated_obs_size,Global.M);
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
    %PS_sampled = NDSort(Population.objs,1);
    %PS_sampled_size = sum(PS_sampled(:) == 1);


    %Intersect PS sampled with PS of design space          
    %int_points1 = intersect(Population.decs, true_pareto_set.decs, 'rows');
    %int_points1_size = size(int_points1,1);
    %PS_s = int_points1_size/design_nondom_size*100;


     % True objectives of PFobs
     %{
     PF_obs_true = zeros(NonDominated_obs_size,Global.M);
     j = 1;
     for i = 1:size(Mat_Obj,1)
         if NonDominated_obs(1,i) == 1
             PF_obs_true(j,:) = Population(:,i).obj;
             j = j+1;
         end
     end
     %}

     %PS observed
     %
     PS_obs = zeros(NonDominated_obs_size,Global.D);
     jj = 1;
     for ii = 1:size(Population.decs,1)
         if NonDominated_obs(1,ii) == 1
             PS_obs(jj,:) = Population(1,ii).dec;
             jj = jj+1;
         end
     end
     %

    %True values of PF sampled (not entirely identified) 
    %{
    PF_pop = zeros(PS_sampled_size,Global.M);
    nonPF_pop = zeros(size_set-PS_sampled_size,Global.M);
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
    %}

    %Error Type I points
     %First remove set found from true PS
     int_points_t1 = true_pareto_set.decs;
     int_points_t1(ismember(int_points_t1,PS_obs,'rows'),:)=[];
     %Second remove set found from all points sampled
     int_points_t1_1 = Population.decs;
     int_points_t1_1(ismember(int_points_t1_1,PS_obs,'rows'),:)=[];
     %Third intersect both sets for points error Type 1
     type1 = intersect(int_points_t1, int_points_t1_1, 'rows');
     type1_size = size(type1,1);
     %Retrieve observed means of type1 points
     %if type1_size ~= 0
        [~,type1_ind] = intersect(Population.decs,type1,'rows');
        type1_obs_obj = Mat_Obj(type1_ind,:);
     %end

     %Error Type II points
     %First remove true points from set found
     type2 = PS_obs;
     type2(ismember(type2,true_pareto_set.decs,'rows'),:)=[];
     %Points error Type 2
     type2_size = size(type2,1);
     %if type2_size ~= 0
        [~,type2_ind] = intersect(Population.decs,type2,'rows');
        type2_obs_obj = Mat_Obj(type2_ind,:);
     %end
     
     %{
     %Percentage of PS that was sampled and correctly identified (PS_ID)
     PS_ID = (int_points1_size-type1_size)/design_nondom_size*100;


     %Intersect PS sampled with PS identified 
    int_points2 = intersect(int_points1, PS_obs, 'rows');
    int_points2_size = size(int_points2,1);
    %APS = (int_points1_size+int_points2_size)/(2*design_nondom_size+type2_size)*100;
    APS = 1-(type1_size+type2_size)/size_set;
    %}


end

