function plot_mors(num_obj,Mat_Obj)

%Retrieve data
%[~,~, ~, PF_obs, nonPF_obs, ~, ~, ~, ~, type1_obs_obj, type1_size, type2_obs_obj, type2_size,~,~,~,~,~,~ ] = sets_vis(Global, ...
%    Population, Mat_Obj);

[PF_obs, nonPF_obs] = sets_simple(num_obj,Mat_Obj);

if num_obj == 2
        %Axis
        %{
        if strcmp(test_function,'DTLZ7')
            size_axis = [0 0.91 2.3 4.1];
        elseif (strcmp(test_function,'WFG1')) || (strcmp(test_function,'WFG3')) || (strcmp(test_function,'WFG4'))
            size_axis = [0 3 0 5];
        else
            size_axis = [-0.1 1.1 -0.1 1.1];
        end
        %}
        
        

         % Plot observed means 
         %{
         figure;
         Draw(PF_obs);
         axis(size_axis);
         title(sprintf('%s on %s(observed means)',func2str(Global.algorithm),func2str(Global.problem)),'Interpreter','none');
         pause;
         %}

         %PFobs + type1 + type2 
         %figure;
         plot(-PF_obs(:,1),-PF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','k');
         hold on
         plot(-nonPF_obs(:,1),-nonPF_obs(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         %hTitle = title('PF');
         %set(hTitle,'FontSize',20)
         set(gca,'FontSize',30);
         %axis(size_axis);
         legend('Pareto front','Dominated points');
         xlabel('EHVD'); ylabel('PD');
         hold off
         %pause;

         %PF sampled vrs PF obs 
         %{
         figure;
         plot(PF_pop(:,1),PF_pop(:,2),'o','MarkerSize',22,'MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0,0,0]);
         hold on
         plot(nonPF_pop(:,1),nonPF_pop(:,2),'o','MarkerSize',22,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         axis(size_axis);
         hTitle = title('True values of all points ');
         set(hTitle,'FontSize',40)
         set(gca,'FontSize',30);
         legend('Pareto front','Dominated points','Location','southeast');
         xlabel('f_1'); ylabel('f_2');
         hold off
         %pause;
         %}
         
         

         %True PF vrs PF sampled vrs PFobs 
         %{
         figure;
         plot(PF_pop(:,1),PF_pop(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','k');
         hold on
         %plot(true_PF(:,1),true_PF(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',8,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         title('True PF vrs PF sampled vrs PFobs');
         axis(size_axis);
         legend('True PF', 'Sampled PF, PF_{obs}');
         hold off
         pause;
         %}

         %PF_obs vrs true obj values 
         %{
         figure;
         plot(PF_obs_true(:,1),PF_obs_true(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','k');
         hold on
         %plot(true_PF(:,1),true_PF(:,2),'o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
         plot(PF_obs(:,1),PF_obs(:,2),'o','MarkerSize',8,'MarkerEdgeColor',[.4 .4 .4],'LineWidth',1,'MarkerFaceColor',[0.7,0.7,0.7]);
         title('True values of PF_{obs} vrs. PF_{obs}');
         axis(size_axis);
         legend('True values of PF_{obs}','PF_{obs}');
         xlabel('f^1'); ylabel('f^2');
         hold off
         pause;
         %}
end


end

