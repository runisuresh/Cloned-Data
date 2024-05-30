function order_geom_actual = calc_tree_equil2(order_geom, mat_tree, mat_temp, n_orders, hemo)
%Find true loaded geometries with structure of variable material properties
%applied to mat_temp for each vessel order.
    
    order_geom_actual = zeros(4, n_orders);
    %order_stress
    for ord = n_orders: -1: 1

        %Assign new material properties to mat_temp for each order
        %Different matrix fractions
        mat_temp.phie = mat_tree.phi_alpha_order(ord, 1);
        mat_temp.phim = mat_tree.phi_alpha_order(ord, 2);
        mat_temp.phic = mat_tree.phi_alpha_order(ord, 3);
        %Different passive SMC contribution
        mat_temp.cm(1) = mat_tree.cm1_order(ord);
        %Different active stress contribution
        mat_temp.T_act = mat_tree.mat_order_T_act(ord);
        %Different elastin material parameter contribution
        mat_temp.ce = mat_tree.mat_order_ce(ord);

        P_iv = hemo.P_order(ord,1) / 10;
        ri_curr = order_geom(1, ord);
        h_curr = order_geom(2, ord);
        ro_curr = ri_curr + h_curr;
        L_curr = order_geom(3, ord);
        V_curr = order_geom(4, ord);
        lz_curr = 1.0;
        F_iv = ones(3,1);
        objective = @(ro)loaded_obj(ro, mat_temp, P_iv, ri_curr, h_curr, V_curr, L_curr, lz_curr, F_iv);
        ro_g = ro_curr;
        ro_actual = fzero(objective, ro_g);
        if ro_actual ~= ro_actual
            ro_actual = ro_curr;
            display("WARNING: Previous iter geom for order " + num2str(ord))
        end
        ri_actual = sqrt(ro_actual^2 - V_curr/(pi * lz_curr * L_curr));
        h_actual = ro_actual - ri_actual;

        order_geom_actual(1, ord) = ri_actual;
        order_geom_actual(2, ord) = h_actual;
        order_geom_actual(3, ord) = L_curr;
        order_geom_actual(4, ord) = V_curr;
        %t_actual = calc_stress(ro_act, mat, ri_act, h_act, V, L_ivref, lz_ivref, F_ivref);

    end

end