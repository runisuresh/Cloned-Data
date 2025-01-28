function fig = plot_hemo(hemo, n_orders, fig_num)
    
 fig = figure(fig_num);
 subplot(2,2,1)
 sp1 = bar(1:n_orders, hemo.P_order(:,1) / (10 * 133.33),'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
 xlabel('Vessel Order (-)'); ylabel('Pressure (mmHg)')
 hold on
 er = errorbar(1:n_orders, hemo.P_order(:,1)/ (10 * 133.33), 0 * hemo.P_order(:,2)/ (10 * 133.33), hemo.P_order(:,2)/ (10 * 133.33));
 er.Color = [0 0 0];
 er.LineStyle = 'none';
%  ylim([0 20])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,2)
 sp2 = bar(1:n_orders, hemo.Q_order(:,1) * 60,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
 xlabel('Vessel Order (-)'); ylabel('Flow (mL/min)')
 hold on
 er = errorbar(1:n_orders, hemo.Q_order(:,1) * 60, 0 * hemo.Q_order(:,2), hemo.Q_order(:,2) * 60);
 er.Color = [0 0 0];
 er.LineStyle = 'none';
 set(gca, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,3)
 sp3 = bar(1:n_orders, hemo.Sigma_order(:,1) / 10000,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
 xlabel('Vessel Order (-)'); ylabel('Circ. Stress (kPa)')
 hold on
 er = errorbar(1:n_orders, hemo.Sigma_order(:,1) / 10000, 0 * hemo.Sigma_order(:,2) / 10000, hemo.Sigma_order(:,2) / 10000);
 er.Color = [0 0 0];
 er.LineStyle = 'none';
%  ylim([0 20])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,4)
 sp4 = bar(1:n_orders, hemo.WSS_order(:,1),'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
 xlabel('Vessel Order (-)'); ylabel('Wall Shear Stress (dyne/cm^2)')
 hold on
 er = errorbar(1:n_orders, hemo.WSS_order(:,1), 0 * hemo.WSS_order(:,2), hemo.WSS_order(:,2));
 er.Color = [0 0 0];
 er.LineStyle = 'none';
%  ylim([0 20])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 set(gcf,'Units','inches','Position',[.25 .25 12 8],'PaperUnits','inches','PaperPosition',[.25 .25 12 8])

end