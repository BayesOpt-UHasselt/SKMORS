function [linear_const] = heter_noise_SK(Global, Candidates, level, caso, reps)

    %level of noise: 1: low, 2: high
    %noise case: 1: best case, 2: worst case

    candid_objs = Candidates.objs;
    linear_const = zeros(Global.M,2);
    noise_sd_range = zeros(Global.M,2);
    
        for i = 1:Global.M
            %Objective i
            obj_i = candid_objs(:,i);
            
            [rank_obj_i] = sort(obj_i,'descend');
            range_obj_i = abs(rank_obj_i(1)-rank_obj_i(end));
            range = ['Objective range', num2str(i), ' = ',num2str(range_obj_i)];
            disp(range);
            
            %Bounds for noise s.d.
            if level == 1
                lower_obj_i = 0.001/sqrt(reps)*range_obj_i;
                upper_obj_i = 0.5/sqrt(reps)*range_obj_i;
            elseif level == 2
                lower_obj_i = 1/sqrt(reps)*range_obj_i;
                upper_obj_i = 2/sqrt(reps)*range_obj_i;
            else
                lower_obj_i = 0.01/sqrt(reps)*range_obj_i;
                upper_obj_i = 1/sqrt(reps)*range_obj_i;
            end
                
            %Noise cases
            if caso == 1 %Best case
                b_obj_i = (rank_obj_i(1,1)*lower_obj_i-rank_obj_i(end,1)*upper_obj_i)/(upper_obj_i-lower_obj_i);
                a_obj_i = lower_obj_i/(rank_obj_i(end,1)+b_obj_i);
            else %Worst case
                b_obj_i = (rank_obj_i(end,1)*lower_obj_i-rank_obj_i(1,1)*upper_obj_i)/(upper_obj_i-lower_obj_i);
                a_obj_i = lower_obj_i/(rank_obj_i(1,1)+b_obj_i);
            end
            
            min_noise = a_obj_i*rank_obj_i(end,1)+a_obj_i*b_obj_i;      
            max_noise = a_obj_i*rank_obj_i(1,1)+a_obj_i*b_obj_i;
            
            linear_const(i,1) = a_obj_i;
            linear_const(i,2) = b_obj_i;
            noise_sd_range(i,1) = min_noise;
            noise_sd_range(i,2) = max_noise;
        end
        disp('a and b constants for the noise s.d.');
        disp(linear_const);
        disp('min and max s.d. for the noise s.d.');
        disp(noise_sd_range);
        %pause;
end

