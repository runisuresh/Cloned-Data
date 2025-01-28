function mat_tree = prescr_mat_tree2(par, mat, n_orders)
%Currently modifies collagen mass fraction, smooth muscle passive
%behavior, smooth muscle active behavior, elastin material parameter, and 
%elastin fraction according to vector scaling factor par.
    
    %Material property changes down the tree
    mat_order_top = mat;
    %Decrease collagen content as we go down the tree
    % mat_order_phic = mat_order_top.phic * linspace(1, 0.05, n_orders); %Decreasing collagen content to 5% of original amount
    mat_order_phic = zeros(1, n_orders);
    mat_order_phic(n_orders) = mat_order_top.phic;
    %Change SMC passive behavior
    mat_order_cm1 = zeros(1, n_orders);
    mat_order_cm1(n_orders) = mat_order_top.cm(1);
%     mat_order_cm2 = zeros(1, n_orders);
%     mat_order_cm2(n_orders) = mat_order_top.cm(2);
    %Change active tone
    % mat_order_T_act = flip(mat_order_top.T_act * [1.0 0.8 0.75 0.7 0.65 0.4 0.3 0.2 ]); % mat_order_top.T_act * ones(1, n_orders); %Decreasing active stress to 10% of original amount
    mat_order_T_act = zeros(1, n_orders);
    mat_order_T_act(n_orders) = mat_order_top.T_act;
    %Change elastin stiffness
    % mat_order_ce = flip(mat_order_top.ce * linspace(1, 0.5, n_orders)); % mat_order_top.ce * ones(1, n_orders); Decreasing elastin material parameter to 10% of original amount
    mat_order_ce = zeros(1, n_orders);
    mat_order_ce(n_orders) = mat_order_top.ce;
    %Change elastin mass fraction
    mat_order_phie = zeros(1, n_orders);
    mat_order_phie(n_orders) = mat_order_top.phie;
    
    %Loop through lower orders
    for ord = n_orders - 1:-1:1
        mat_order_phic(ord) = mat_order_phic(n_orders) * exp(par(1) *  (ord - n_orders));
        mat_order_cm1(ord) = mat_order_cm1(ord + 1) * par(2);
%         mat_order_cm2 = mat_order_cm2(ord + 1) * par(2);
        %mat_order_T_act(ord) = mat_order_T_act(ord + 1) * par(3);
        mat_order_T_act(ord) = 0.9 * mat_order_T_act(n_orders) * exp(par(3) *  (ord - n_orders)) + 0.1 * mat_order_T_act(n_orders);% * (1 - par(3) * (n_orders - ord));
        mat_order_ce(ord) = mat_order_ce(n_orders) * exp(par(4) *  (ord - n_orders));
        mat_order_phie(ord) = mat_order_phie(ord + 1) * par(5); % exp( *  (ord - n_orders));
    end

    mat_order_phic_diff = mat_order_top.phic - mat_order_phic; %Find col frac to be added to smc fractions
%     mat_order_phie_diff = mat_order_top.phie - mat_order_phie; %Find eln frac to be added to smc fractions
%     mat_order_phim = mat_order_top.phim + mat_order_phic_diff + mat_order_phie_diff; % + mat_order_phie_diff;
%     eln_vs_smc_frac = 0.0; %1.0 for all eln, 0.0 for all smc
    mat_order_phie = mat_order_top.phie + mat_order_phic_diff / 2; % eln_vs_smc_frac * mat_order_phic_diff; %+ mat_order_phic_diff + mat_order_phie_diff;
%     mat_order_phim = mat_order_top.phim + (1 - eln_vs_smc_frac) * mat_order_phic_diff; % + mat_order_phie_diff;
    mat_order_phim = mat_order_top.phim + mat_order_phic_diff / 2; % + mat_order_phie_diff; % + mat_order_phie_diff;
    
    phi_alpha_order = [mat_order_phie; mat_order_phim; mat_order_phic]';
    
    mat_tree.phi_alpha_order = phi_alpha_order;
    mat_tree.cm1_order = mat_order_cm1;
    mat_tree.mat_order_T_act = mat_order_T_act;
    mat_tree.mat_order_ce = mat_order_ce;

end