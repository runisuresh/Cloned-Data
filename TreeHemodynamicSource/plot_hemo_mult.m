function fig = plot_hemo_mult(hemo_all, n_orders, fig_num)

 [temp, num_hemo] = size(hemo_all);
 n_orders_all = ones(n_orders, num_hemo);
 n_orders_all = n_orders_all .* [1:n_orders]';
 
 P_all = zeros(n_orders, num_hemo);
 Q_all = zeros(n_orders, num_hemo);
 Sigma_all = zeros(n_orders, num_hemo);
 WSS_all = zeros(n_orders, num_hemo);
 for i = 1:num_hemo
     P_all(:,i) = hemo_all(i).P_order(:,1) / (10 * 133.33);
     Q_all(:,i) = hemo_all(i).Q_order(:,1) * 60;
     Sigma_all(:,i) = hemo_all(i).Sigma_order(:,1) / 10000;
     WSS_all(:,i) = hemo_all(i).WSS_order(:,1);
 end
 
 fig = figure(fig_num);
 subplot(2,2,1)
 sp1 = bar(n_orders_all, P_all,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('Pressure (mmHg)')
 hold on
%  er = errorbar(1:n_orders, hemo.P_order(:,1)/ (10 * 133.33), 0 * hemo.P_order(:,2)/ (10 * 133.33), hemo.P_order(:,2)/ (10 * 133.33));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
% sp1(1).FaceColor = 'k';
% sp1(2).FaceColor = 'r';
% sp1(3).FaceColor = 'b';
 ylim([0 20])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,2)
 sp2 = bar(n_orders_all, Q_all,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('Flow (mL/min)')
 hold on
%  er = errorbar(1:n_orders, hemo.Q_order(:,1), 0 * hemo.Q_order(:,2), hemo.Q_order(:,2));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
% sp2(1).FaceColor = 'k';
% sp2(2).FaceColor = 'r';
% sp2(3).FaceColor = 'b';
 set(gca, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,3)
 sp3 = bar(n_orders_all, Sigma_all,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('IMS (kPa)')
 hold on
%  er = errorbar(1:n_orders, hemo.WSS_order(:,1), 0 * hemo.WSS_order(:,2), hemo.WSS_order(:,2));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
% sp3(1).FaceColor = 'k';
% sp3(2).FaceColor = 'r';
% sp3(3).FaceColor = 'b';
% legend('Homeostatic', 'Incr Flow', 'After G&R')
%  ylim([5 20])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,4)
 sp4 = bar(n_orders_all, WSS_all,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 xlabel('Vessel Order (-)'); ylabel('WSS (dyne/cm^2)')
 hold on
%  er = errorbar(1:n_orders, hemo.WSS_order(:,1), 0 * hemo.WSS_order(:,2), hemo.WSS_order(:,2));
%  er.Color = [0 0 0];
%  er.LineStyle = 'none';
% sp4(1).FaceColor = 'k';
% sp4(2).FaceColor = 'r';
% sp4(3).FaceColor = 'b';
%  ylim([40 125])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 set(gcf,'Units','inches','Position',[.25 .25 12 8],'PaperUnits','inches','PaperPosition',[.25 .25 12 8])

end