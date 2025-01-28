function J = ord_diam_obj(par, input_pars)

    global diameter thickness length new_conn
    
    n_orders = input_pars.n_orders;
    
    %Diam ratio is diam order n / n - 1 -> d_n-1 = d_n / d_ratio
    %Thick ratio is thick order n / diam order n -> h_n-1 = d_n * h_ratio
    %Len ratio is len order n / n - 1 -> l_n-1 = l_n / l_ratio
    diam_ratio_const = par;
    order_diam_ratio = diam_ratio_const * ones(n_orders - 1,1); %can vary along the orders, but constant for now, reasonably linear in Jiang
    thick_ratio_const = input_pars.h_iv_LPA / (2 * input_pars.r_iv_LPA); %calculated from in vivo gemoetry w/ active stress at 15 mmHg from data in Ramachandra and Humphrey 2019, J Biomech
    order_thick_ratio = thick_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, some evidence of constant in humans after infancy, see Figure 1 Haworth et al., 1983 Am J Cardiol
    len_ratio_const = input_pars.len_ratio; %from regression in Figure 5, Jiang et al., 1994 J App Physio
    order_len_ratio = len_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, reasonably linear in Jiang

    %Initialize highest order values
    order_geom(1, n_orders) = input_pars.r_iv_LPA;
    order_geom(2, n_orders) = input_pars.h_iv_LPA;
    order_geom(3, n_orders) = input_pars.l_iv_LPA;
    order_geom(4, n_orders) = input_pars.V_LPA;
    %Loop through lower orders
    for ord = n_orders - 1: -1: 1
        order_geom(1, ord) = order_geom(1, ord + 1) / order_diam_ratio(ord);
        order_geom(2, ord) = order_thick_ratio(ord) * 2 * order_geom(1, ord);
        order_geom(3, ord) = order_geom(3, ord + 1) / order_len_ratio(ord);
        order_geom(4, ord) = pi * ( (order_geom(1, ord) + order_geom(2, ord))^2 - order_geom(1, ord)^2 ) * order_geom(3, ord);
    end

    diameter = 2 * order_geom(1, :) * 1D2; %convert to cm
    thickness = 2 * order_geom(2, :) * 1D2; %convert to cm
    length = order_geom(3, :) * 1D2; %convert to cm

    %Calculate hemodynamics for given orders
    if (isfile('morphometric_tree.mat') && new_conn == 0)
        gen = 0;
    elseif (new_conn == 1) 
        gen = 1;
        new_conn = 0;
    else
        gen = 1;
    end
    
    order = n_orders;
    order_remodeling = n_orders;
    alpha = 1;
    beta = 1;
    PVR_LPA_tree = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
    
    J = PVR_LPA_tree - input_pars.PVR_LPA;

end

