function fig = plot_hemo_mult_init(hemo_all, n_orders, fig_num)

    [temp, num_hemo] = size(hemo_all);
    n_orders_all = ones(n_orders, num_hemo);
    n_orders_all = n_orders_all .* [1:n_orders]';

    P_all = zeros(n_orders, num_hemo);
    P_all_err = zeros(n_orders, num_hemo);
    Q_all = zeros(n_orders, num_hemo);
    Q_all_err = zeros(n_orders, num_hemo);
    Sigma_all = zeros(n_orders, num_hemo);
    WSS_all = zeros(n_orders, num_hemo);
    for i = 1:num_hemo
        P_all(:,i) = hemo_all(i).P_order(:,1) / (10 * 133.33);
        P_all_err(:,i) = hemo_all(i).P_order(:,2) / (10 * 133.33);
        Q_all(:,i) = hemo_all(i).Q_order(:,1) * 60;
        Q_all_err(:,i) = hemo_all(i).Q_order(:,2) * 60;
        Sigma_all(:,i) = hemo_all(i).Sigma_order(:,1) / 10000;
        WSS_all(:,i) = hemo_all(i).WSS_order(:,1);
    end

    fig = figure(fig_num);
    for hemo_n = 1:num_hemo

    subplot(1,2,1)
    sp1 = bar(n_orders_all, P_all, 1.0,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);

    % Get the x coordinate of the bars
    x = nan(num_hemo, n_orders);
    for i = 1:num_hemo
        x(i,:) = sp1(i).XEndPoints;
    end

    xlabel('Vessel Order (-)'); ylabel('Pressure (mmHg)')
    hold on
    er = errorbar(x', P_all, 0 * P_all_err, P_all_err, 0 * P_all_err, 0 * P_all_err, 'k','linestyle','none','LineWidth',1.5);

    sp1(1).FaceColor = 'k';
    sp1(2).FaceColor = [0.5 0.5 0.5];
    % sp1(3).FaceColor = 'b';
    ylim([0 15])
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold', 'Ytick', 0:3:15)

    subplot(1,2,2)
    sp2 = bar(n_orders_all, Q_all, 1.0,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);

    % Get the x coordinate of the bars
    x = nan(num_hemo, n_orders);
    for i = 1:num_hemo
        x(i,:) = sp2(i).XEndPoints;
    end

    xlabel('Vessel Order (-)'); ylabel('Flow (mL/min)')
    hold on
    er2 = errorbar(x', Q_all, 0 * Q_all_err, Q_all_err, 0 * Q_all_err, 0 * Q_all_err,'k','linestyle','none','LineWidth',1.5);

    sp2(1).FaceColor = 'k';
    sp2(2).FaceColor = [0.5 0.5 0.5];
    ylim([10^-5 10])
    set(gca, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

    end
    set(gcf,'Units','inches','Position',[.25 .25 12 4],'PaperUnits','inches','PaperPosition',[.25 .25 12 4])
end