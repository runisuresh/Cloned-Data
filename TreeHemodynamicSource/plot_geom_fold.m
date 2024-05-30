function fig = plot_geom_fold(order_geom_orig, order_geom_new, n_orders, fig_num)

 fold_r = order_geom_new(1,:) ./ order_geom_orig(1,:);
 fold_h = order_geom_new(2,:) ./ order_geom_orig(2,:);
 
 fig = figure(fig_num);
 subplot(1,2,1)
 sp1 = bar(1:n_orders, fold_r,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('Fold Radius (-)')
 hold on
%  er = errorbar(1:n_orders, hemo.P_order(:,1)/ (10 * 133.33), 0 * hemo.P_order(:,2)/ (10 * 133.33), hemo.P_order(:,2)/ (10 * 133.33));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(1,2,2)
 sp2 = bar(1:n_orders, fold_h,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('Fold Thickness (-)')
 hold on
%  er = errorbar(1:n_orders, hemo.Q_order(:,1), 0 * hemo.Q_order(:,2), hemo.Q_order(:,2));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 set(gcf,'Units','inches','Position',[.25 .25 12 8],'PaperUnits','inches','PaperPosition',[.25 .25 12 8])
 

end