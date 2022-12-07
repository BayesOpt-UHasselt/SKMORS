%function [size_set, true_pareto_set, true_PF, PF_obs, nonPF_obs, PF_obs_true, PS_obs, PF_pop, nonPF_pop, type1_obs_obj, type1_size, ...
%    type2_obs_obj, type2_size, design_nondom_size, NonDominated_obs_size, PS_sampled_size, PS_s, PS_ID, APS ] = sets_vis(Global, ...
%    Population, Mat_Obj)

function [PF_obs, nonPF_obs] = sets_simple(num_obj,Mat_Obj)

     size_set = size(Mat_Obj,1);
     
      %Observed means and points
     NonDominated_obs = NDSort(Mat_Obj,1);
     NonDominated_obs_size = sum(NonDominated_obs(:) == 1); 
     

     %PF observed 
     PF_obs = zeros(NonDominated_obs_size,num_obj);
     nonPF_obs = zeros(size_set-NonDominated_obs_size,num_obj);
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
     

    

end

