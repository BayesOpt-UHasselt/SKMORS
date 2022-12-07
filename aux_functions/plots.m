function plots(Global,Population,test_function, Mat_Obj, Mat_Var, num)

%Retrieve data
%[~,~, ~, PF_obs, nonPF_obs, ~, ~, PF_pop, nonPF_pop, type1_obs_obj, type1_size, type2_obs_obj, type2_size,~,~,~,~,~,~ ] = sets(Global, ...
%    Population, Mat_Obj,Mat_Var);

[~, ~, ~, PF_obs, nonPF_obs, ~, ~, ~, PF_pop, nonPF_pop, PF_mocba, nonPF_mocba,type1_obs_obj, type1_obs_obj_mocba, type1_size, ...
    type1_size_mocba, type2_obs_obj, type2_obs_obj_mocba, type2_size, type2_size_mocba, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~ ] = ...
    sets(Global, Population, Mat_Obj, Mat_Var);

%Axis
if strcmp(test_function,'DTLZ7')
    size_axis = [0 0.91 2.3 4.1];
elseif (strcmp(test_function,'WFG1')) || (strcmp(test_function,'WFG3')) || (strcmp(test_function,'WFG4'))
    size_axis = [0 3 0 5];
else
    size_axis = [-0.1 1.1 -0.1 1.1];
end

if num == 1 %Normal errors
         %Observed performance 
         figure;
         %plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0,0,0]);
         hold on
         plot(nonPF_obs(:,1),nonPF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         if type1_size ~= 0
            plot(type1_obs_obj(:,1),type1_obs_obj(:,2),'^','MarkerSize',24,'MarkerEdgeColor','b','LineWidth',1,'MarkerFaceColor','b');
         end
         if type2_size ~= 0
            plot(type2_obs_obj(:,1),type2_obs_obj(:,2),'s','MarkerSize',24,'MarkerEdgeColor','r','LineWidth',1,'MarkerFaceColor','r');
         end
         %{
         titl = ['PF vs. MCE and MCI on ',test_function];
         hTitle = title(titl);
         set(hTitle,'FontSize',40)
         %}
         set(gca,'FontSize',30);
         axis(size_axis);
         legend('PF','Dominated points', 'MCE', 'MCI');
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;

         %True performance
         %
         figure;
         plot(PF_pop(:,1),PF_pop(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0,0,0]);
         hold on
         plot(nonPF_pop(:,1),nonPF_pop(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         axis(size_axis);
         %{
         titl = ['True performance on ',test_function];
         hTitle = title(titl);
         set(hTitle,'FontSize',40)
         %}
         set(gca,'FontSize',30);
         legend('Pareto front','Dominated points','Location','northeast');
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;
         %
elseif num == 2 %MOCBA errors 
         %Observed performance 
         figure;
         %plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         plot(PF_mocba(:,1),PF_mocba(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0,0,0]);
         hold on
         plot(nonPF_mocba(:,1),nonPF_mocba(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         if type1_size_mocba ~= 0
            plot(type1_obs_obj_mocba(:,1),type1_obs_obj_mocba(:,2),'^','MarkerSize',24,'MarkerEdgeColor','b','LineWidth',1,'MarkerFaceColor','b');
         end
         if type2_size_mocba ~= 0
            plot(type2_obs_obj_mocba(:,1),type2_obs_obj_mocba(:,2),'s','MarkerSize',24,'MarkerEdgeColor','r','LineWidth',1,'MarkerFaceColor','r');
         end
         titl = ['MOCBA: PF vs. MCE and MCI on ',test_function];
         hTitle = title(titl);
         set(hTitle,'FontSize',40)
         set(gca,'FontSize',30);
         axis(size_axis);
         legend('PF','Dominated points', 'MCE', 'MCI');
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;

         
end


end

