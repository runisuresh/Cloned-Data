function J = tree_mat_obj3(par, order_geom, mat_temp, n_orders, hemo)
%Objective function for finding variable material properties down the PA
%tree. Modifies according to mat_tree.
    
    mat_tree = prescr_mat_tree2(par, mat_temp, n_orders); %exp changes in mat props 
    %mat_tree = prescr_mat_tree3(par, mat_temp, n_orders); %linear changes in mat props
    
%     cmd_flags = "-u 0 --simulate_equil 0";
%     output_flags = 0;
%     order_geom_actual = calc_tree_equil_wgnr(order_geom, mat_tree, mat, n_orders, hemo, cmd_flags, output_flags);

    J = zeros(1, n_orders);
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
        t2_LP = P_iv * ri_curr / h_curr;
        ro_curr = ri_curr + h_curr;
        L_curr = order_geom(3, ord);
        V_curr = order_geom(4, ord);
        lz_curr = 1.0;
        F_iv = ones(3,1);        
        t = calc_stress(ro_curr, mat_temp, ri_curr, h_curr, V_curr, L_curr, lz_curr, F_iv);
        
        J(ord) = (t(2) - t2_LP) / t2_LP;
        
    end
    


end