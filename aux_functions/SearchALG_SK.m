function [PopDec] = SearchALG_SK(PCheby,Dec,model, candid_set)
% Solution update. Returns the point with the best modified expected improvement

    N     = size(candid_set,1);
    %disp(size(candid_set));
    %pause;
    %D     = size(candid_set,2);
    
    %Calculate MEI
   [Gbest,Gbest_ind] = min(PCheby);
   xcurr_best=Dec(Gbest_ind,:);
    
   Z_min=SKpredict(model,xcurr_best,ones(1,1)); % the predicted response at the sampled point with 
                                                                                % lowest sampled mean (xcurr_best)                                                                                                                       
   MEI  = zeros(N,1); 
   EI  = zeros(N,1); 
    
   Mat = zeros(N,3); %Result matrix (MEI, response and EI)
                 
   %Iterate over all search space
   for i = 1 : N
       [y,~,mse] = SKpredict_new(model,candid_set(i,:),1); %mse is the deterministic kriging variance
       mse = sqrt(abs(mse));
       %s = sqrt(abs(s));
            
      MEI(i) = (Z_min-y) * normcdf((Z_min-y)/mse,0,1) + mse * normpdf((Z_min-y)/mse,0,1);
      %EI(i)     = (Gbest-y)*normcdf((Gbest-y)/s)+s*normpdf((Gbest-y)/s);
        
      Mat(i,1) = MEI(i); %Save MEI
      Mat(i,2) = y; %Save response
      %Mat(i,3) = EI(i); %Save EI
   end
        
   [~,index] = sortrows(Mat,'descend'); %Sort based on MEI (first column of Mat)
   %[~,index] = sortrows(Mat,3,'descend'); %Sort based on EI (third column of Mat)
   Best_point = candid_set(index(1),:); 
    
   PopDec = Best_point;
end
