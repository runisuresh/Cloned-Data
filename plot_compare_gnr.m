function fig = plot_compare_gnr(max_days, n_orders, step_size, fig_num)

colors = pink(n_orders + 1);

% colors = [0.0 0.0 0.0
%           0.4 0.4 0.4
%           0.8 0.8 0.8];
      
%ones(5,3).*(0.8:-0.2:0)';
%          [0.75 0.75 0.75
%           0.6 0.6 0.6
%           0.5 0.5 0.5
%           0.25 0.25 0.25
%           0.0 0.0 0.0];
data_files = cell(1,n_orders);
for ord = 1:n_orders
    data_files{ord} = "GnR_out_ord" + ord;
end

%Parameterized inputs
n_par = length(data_files);
time = 0:step_size:max_days;

for i = 1:n_orders
    
    data = load(data_files{i});
    
    a = data(1:length(time),1);
    h = data(1:length(time),2);
    rhoR = data(1:length(time),3);
    
    P = data(1:length(time),9);
    P_h = data(1:length(time),10);
    
    sigma_t = P.*a./h;
    sigma_t_h = P_h(1)*a(1)/h(1);
    delta_sigma_t = (sigma_t - sigma_t_h) / sigma_t_h;
    
    bar_tauw = data(1:length(time),7);
    bar_tauw_h = data(1,8);
    delta_tauw = (bar_tauw - bar_tauw_h) / bar_tauw_h;
    
    rhoR_p1 = data(1:length(time),4);
    rhoR_p2 = data(1:length(time),5);
    rhoR_i = data(1:length(time),6);

    fig = figure(fig_num);
    subplot(2,2,1)
    hold on
    plot(time / 7, a/a(1), 'LineWidth', 2.0, 'Color', colors(end - i,:))
    % plot(time(end), a_e/a(1), 's')
    xlabel('Time (weeks)'); ylabel('radius (-)')
%     axis([ 0 10 -inf inf])
    set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold', 'Xtick', 0:2:10)

    subplot(2,2,2)
    hold on
    plot(time / 7, h/h(1), 'LineWidth', 2.0, 'Color', colors(end - i,:))
    % plot(time(end), h_e/h(1), 's')
    xlabel('Time (weeks)'); ylabel('thickness (-)')
%     axis([ 0 10 -inf inf])
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold', 'Xtick', 0:2:10)

%     subplot(2,2,3)
%     hold on
%     plot(time / 7, rhoR, 'LineWidth', 1.5, 'Color', colors(end - i,:))
%     % plot(time(end), h_e/h(1), 's')
%     xlabel('Time (weeks)'); ylabel('rhoR (kg/m^3)')
%     % axis([ -50 560 0.8 2.0])
%     set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

    subplot(2,2,3)
    hold on
    plot(time / 7, delta_sigma_t, 'LineWidth', 2.0, 'Color', colors(end - i,:))
    % plot(time(end), h_e/h(1), 's')
    xlabel('Time (weeks)'); 
    ylabel('Change in IMS (-)')
%     ylabel('\Delta\sigma_{\theta} (-)')
%     axis([ 0 10 0 0.4])
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold', 'Xtick', 0:2:10)

    subplot(2,2,4)
    hold on
    plot(time / 7, delta_tauw, 'LineWidth', 2.0, 'Color', colors(end - i,:))
    % plot(time(end), h_e/h(1), 's')
    xlabel('Time (weeks)');
    ylabel('Change in WSS (-)')
%     ylabel('\Delta\tau_{w} (-)')
%     axis([ 0 10 0 0.4])
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold', 'Xtick', 0:2:10)

%     subplot(2,3,6)
%     hold on
%     plot(time / 7, rhoR_i, 'LineWidth', 1.5, 'Color', colors(i,:))
%     % plot(time(end), h_e/h(1), 's')
%     xlabel('Time (weeks)'); ylabel('rho_{i,smc} (kg/m^3)')
%     % axis([ -50 560 0.8 2.0])
%     set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

end
set(gcf,'Units','inches','Position',[.25 .25 12 8],'PaperUnits','inches','PaperPosition',[.25 .25 12 8])

end

