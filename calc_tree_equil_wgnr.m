function order_geom_actual = calc_tree_equil_wgnr(order_geom, mat_tree, mat_temp, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags, params)
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
        %Different active stress contribution
        mat_temp.T_act = mat_tree.mat_order_T_act(ord);
        %Different elastin material parameter contribution
        mat_temp.ce = mat_tree.mat_order_ce(ord);
        ri_iv = order_geom(1, ord) * 1D3; %convert m to mm inputs
        h_iv = order_geom(2, ord) * 1D3; %convert m to mm inputs
        P_iv = hemo.P_order(ord,1) / 10; %%convert dynes to Pa inputs
        Q_iv = hemo.Q_order(ord,1);
        
        %Write native_in file for current order
        write_Native_in(mat_temp, ord, ri_iv, h_iv, P_iv, Q_iv, output_flags, params);
        
        command = "./gnr " + cmd_flags + " --gamma_p " + gamma_p_ord(ord) + " --gamma_q " + gamma_q_ord(ord) + " --gamma_act " + gamma_act_ord(ord);
        if output_flags == 2 || output_flags == 4
            command = command + " -n ord" + ord;
        end
        
        [temp1, temp2] = system(command);
        if (temp1 > 0)
           disp("Warning: G&R model did not converge") 
        end
        
        switch output_flags
            case 0
                gnr_out = load("GnR_out_");
                
                order_geom_actual(1, ord) = gnr_out(1,1);
                order_geom_actual(2, ord) = gnr_out(1,2);
            case 1
                gnr_out = load("Equil_GnR_out_");
                if abs(gnr_out(6) - 1.0) > 0.01
                   disp('equilibrated gnr solution did not converge') 
                end
                
                order_geom_actual(1, ord) = gnr_out(1,1);
                order_geom_actual(2, ord) = gnr_out(1,2);
            case 2
                gnr_out = load("GnR_out_ord" + ord);
                
                order_geom_actual(1, ord) = gnr_out(end,1);
                order_geom_actual(2, ord) = gnr_out(end,2);
            case 3
                gnr_out = load("GnR_out_");
                
                order_geom_actual(1, ord) = gnr_out(end,1);
                order_geom_actual(2, ord) = gnr_out(end,2);
            case 4
                gnr_out = load("Exp_out_ord" + ord);
                
                order_geom_actual(1, ord) = gnr_out(end,1);
                order_geom_actual(2, ord) = gnr_out(end,2);
        end
        
        
        order_geom_actual(3, ord) = order_geom(3, ord);
        order_geom_actual(4, ord) = order_geom(4, ord);
        
    end

end