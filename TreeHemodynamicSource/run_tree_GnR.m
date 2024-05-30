close all; clear all; clc; clf;

global connectivity diameter thickness length new_conn

%--------------------------------------------------------------------------
%G&R simulation settings
%--------------------------------------------------------------------------
tree_diamrat_prev_flag = 1; %Use diameter ratio from a previous solve 1
tree_solve_prev_flag = 1; %Use down the tree material parameters from a previous solve 1
tree_gen_prev_flag = 1; %Use a morphometric tree from a previous solve 1
single_vessels_flag = 0; %Disable hemodynamic feedback 1
mech_infl_flag = " --mech_infl_flag 0"; %Use inflammation driven by WSS 1
dQ = 0.00; %Pertubation in inflow
dP = 0.00; %Pertubation in outlet pressure
k_ramp = 1/4; %days to ramp hemodynamic changes over
step_size = 1.0; %size of step in days
max_days = 180; %maximum days to run the G&R
save_steps = 1; %number of steps between hemo->G&R feedback, min 1
other = "uniforminfl_test2";

%Derived from simulation settings
save_dir = "max" + max_days + "save" + save_steps + "dQ" + dQ + "dP" + dP + other;
ts = 0: step_size: max_days / step_size - 1; %Time steps for the simulation
dQ_s = dQ * (1 - exp( -k_ramp * ts)); %Ramping function for flow %* ones(size(ts));
dP_s = dP * (1 - exp( -k_ramp * ts)); %Ramping function for pressure %* ones(size(ts));

%Check if current run exists
if(exist(save_dir, 'dir'))
    answer = questdlg('Save directory already exists, overwrite?', 'Overwrite?', 'Yes', 'No', 'Yes');
    switch answer
        case 'Yes'
            
        case 'No'
            disp('Stopping exceution, rename save directory')
            return
    end
end

%Update global flag for generating conn matrix
if (tree_gen_prev_flag == 1)
    new_conn = 0;
else
    new_conn = 1;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Setting up morphometric trees for hemodynamic evaluation
%--------------------------------------------------------------------------
%For rat from Jiang et al 1994
connectivity=[...
    0.19 4.06 2.48 1.36 0 0 0 0 0 0 0;
    0 0.12 1.55 0.93 0 0 0 0 0 0 0;
    0 0 0.17 2.26 0.31 0.05 0.11 0.13 0 0 0;
    0 0 0 0.05 2.00 0.70 0.62 0.43 0.39 0 0;
    0 0 0 0 0.18 1.92 1.56 0.85 0.83 0.50 0;
    0 0 0 0 0 0.18 2.05 1.15 1.00 0.67 0;
    0 0 0 0 0 0 0.05 1.55 1.06 0.67 0;
    0 0 0 0 0 0 0 0.08 1.67 0.83 6;
    0 0 0 0 0 0 0 0 0.06 1.67 4;
    0 0 0 0 0 0 0 0 0 0.17 2;
    0 0 0 0 0 0 0 0 0 0 0;
    ];

%To reduce connectivity matrix for mouse, decreasing total vessel numbers
connectivity = connectivity(1:10, 1:10);

%Initialize vessel order metrics
n_orders = size(connectivity, 1); %Number of vessel orders, vessels decrease in size with decreasing order
orders = 1:n_orders;

%Biomechanical data for mouse LPA, from Ramachandra and Humphrey 2019, J Biomech
load('Ramachandra_LPA_fit_mm4ff.mat');

%Extract RPA geometry from structure;
r_iv_LPA = kinem.ri_iv_act;
h_iv_LPA = kinem.h_iv_act;
l_iv_LPA = kinem.F_iv(3) * kinem.L;
V_LPA = kinem.V;

%Prescribed hemodynamics from the literature 
CO = 10.4; %ml/min, From Champion et al., 2000 Am J Heart Circ Physiol
FlowSplit = 0.30; %From Razavi et al., 2011 Am J Physiol Lung Cell Mol Physiol
Q_LPA = CO / 60 * FlowSplit; %ml/s 
mPAP = kinem.P_iv / 133.33; %mmHg, From Champion et al., 2000 Am J Heart Circ Physiol
PCWP =  4.8; %mmHg, From Champion et al., 2000 Am J Heart Circ Physiol
PVR_LPA = (mPAP - PCWP) / Q_LPA * 1333.33;
P_term = PCWP; %mmHg, terminal pressure of the tree

input_pars.n_orders = n_orders;
input_pars.r_iv_LPA = r_iv_LPA;
input_pars.h_iv_LPA = h_iv_LPA;
input_pars.l_iv_LPA = l_iv_LPA;
input_pars.V_LPA = V_LPA;
input_pars.PVR_LPA = PVR_LPA;
input_pars.len_ratio = 1.60;

if tree_diamrat_prev_flag == 0
    lb = 1.0;
    ub = 2.0;
    par0 = 1.4;
    options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-12,'TolX',1e-12,...
                           'Display','iter');
    current_objective = @(par)ord_diam_obj(par, input_pars);
    diam_ratio_const = lsqnonlin(current_objective, par0, lb, ub, options);
else
    temp_geom = load('hemo_and_geom.mat');
    diam_ratio_const = temp_geom.geom_save.order_geom_prescr(1,2) / ...
                       temp_geom.geom_save.order_geom_prescr(1,1);
end


%Diam ratio is diam order n / n - 1 -> d_n-1 = d_n / d_ratio
%Thick ratio is thick order n / diam order n -> h_n-1 = d_n * h_ratio
%Len ratio is len order n / n - 1 -> l_n-1 = l_n / l_ratio
order_diam_ratio = diam_ratio_const * ones(n_orders - 1,1); %can vary along the orders, but constant for now, reasonably linear in Jiang
thick_ratio_const = h_iv_LPA / (2 * r_iv_LPA); %calculated from in vivo gemoetry from data in Ramachandra and Humphrey 2019, J Biomech
order_thick_ratio = thick_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, some evidence of constant in humans after infancy, see Figure 1 Haworth et al., 1983 Am J Cardiol
thick_rat_incr = 1.00; %increasing thickness ratio for each order
len_ratio_const = 1.60; %from regression in Figure 5, Jiang et al., 1994 J App Physio
order_len_ratio = len_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, reasonably linear in Jiang

%Create radius, thickness, length, volume array for each order 
order_geom = zeros(4, n_orders);

%Initialize highest order values
order_geom(1, n_orders) = r_iv_LPA;
order_geom(2, n_orders) = h_iv_LPA;
order_geom(3, n_orders) = l_iv_LPA;
order_geom(4, n_orders) = V_LPA;
%Loop through lower orders
for ord = n_orders - 1: -1: 1
    order_geom(1, ord) = order_geom(1, ord + 1) / order_diam_ratio(ord);
    order_geom(2, ord) = order_thick_ratio(ord) * 2 * order_geom(1, ord) * thick_rat_incr^(n_orders - ord);
    order_geom(3, ord) = order_geom(3, ord + 1) / order_len_ratio(ord);
    order_geom(4, ord) = pi * ( (order_geom(1, ord) + order_geom(2, ord))^2 - order_geom(1, ord)^2 ) * order_geom(3, ord);
end
order_geom_prescr = order_geom;

diameter = 2 * order_geom(1, :) * 1D2; %convert to cm
thickness = order_geom(2, :) * 1D2; %convert to cm
length = order_geom(3, :) * 1D2; %convert to cm

%Calculate hemodynamics for given orders
order = n_orders;
order_remodeling = n_orders;
alpha = 1;
beta = 1;
gen = 0;
resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
hemo = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);
hemo_prescr = hemo;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Determining the solid material parameters down the tree to maintain
%prescribed geometry
%--------------------------------------------------------------------------
%Check if using previous material parameter solution
if tree_solve_prev_flag == 0
    
    %Iteratively update the material parameters for the calculated hemodynamics
    %for the tree
    diff = 1;
    diff_max = 0.01; %tolerance for pressure diff with updates

    count = 0;
    count_max = 10;  %max number of hemodynamic updates

    %Intialize preliminary initial guess for material parameter variations down
    %the tree
    par0 =  0.8 + [rand rand rand rand rand] * (1 - 0.8) - [rand rand rand rand rand] * (1 - 0.95); %[0.8 1.0 0.8 0.8 1.0]; %
    
    while (diff > diff_max && count < count_max)
        %Recalculate the loaded geometry based on hemodynamics from the tree
        %assessment and find variable wall properties that support the geometry.
        if (count > 0)
            par0 = par_est;
        end

        lb = [0.0 1.0 0.00 0.0 1.0];
        ub = [1.0 1.0 1.0 1.0 1.0];

        options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-12,'TolX',1e-12,...
                               'Display','iter');
        current_objective = @(par)tree_mat_obj3(par, order_geom, mat, n_orders, hemo);
        par_est = lsqnonlin(current_objective, par0, lb, ub, options);
%         par_est =  [0 1 0 0 1]; %  uncomment for same material params down the tree
        
        mat_tree = prescr_mat_tree2(par_est, mat, n_orders);       
        order_geom_actual = calc_tree_equil2(order_geom, mat_tree, mat, n_orders, hemo);

        diameter = 2 * order_geom_actual(1, :) * 1D2; %convert to cm
        thickness = order_geom_actual(2, :) * 1D2; %convert to cm

        %Recalculate hemodynamics for given orders
        order = n_orders;
        order_remodeling = n_orders;
        alpha = 1;
        beta = 1;
        gen = 0;
        resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
        hemo_temp = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);
        
%         fig_num = 1;
%         fig1 = plot_hemo(hemo_temp, n_orders, fig_num);
        
%         %Insert breakpoint on next line to get plot of hemo with const material properties
        diff = max(abs((hemo_temp.P_order(:,1) - hemo.P_order(:,1))./ hemo.P_order(:,1)));
        
        hemo = hemo_temp;
        count = count + 1;
    end
    if count == count_max
       disp("Warning: Parameter estimation to match prescribed geometries failed") 
    end
else
    load("mat_tree_props.mat");
    order_geom_actual = calc_tree_equil2(order_geom, mat_tree, mat, n_orders, hemo);

    diameter = 2 * order_geom_actual(1, :) * 1D2; %convert to cm
    thickness = order_geom_actual(2, :) * 1D2; %convert to cm

    %Recalculate hemodynamics for given orders
    order = n_orders;
    order_remodeling = n_orders;
    alpha = 1;
    beta = 1;
    gen = 0;
    resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
    hemo = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);
    
end

%Use GnR code to prescribe initial tree geometries
cmd_flags = "-u 0 --simulate_equil 0";
output_flags = 0;
order_geom_actual_temp = order_geom;
count = 0;
count_max = 25;
diff = 1;
diff_max = 0.00001; %tolerance for pressure diff with updates
gamma_p_ord = zeros(1,n_orders);
gamma_q_ord = zeros(1,n_orders);
gamma_act_ord = zeros(1,n_orders);
while (diff > diff_max && count < count_max)
    order_geom_actual = calc_tree_equil_wgnr(order_geom_actual_temp, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);
    
    diameter = 2 * order_geom_actual(1, :) * 1D2; %convert to cm
    thickness = order_geom_actual(2, :) * 1D2; %convert to cm
    
    %Recalculate hemodynamics for given orders
    order = n_orders;
    order_remodeling = n_orders;
    alpha = 1;
    beta = 1;
    gen = 0;
    resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
    hemo = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);

    diff = max((order_geom_actual(:,1) - order_geom_actual_temp(:,1)) ./ order_geom_actual_temp(:,1));
    order_geom_actual_temp = order_geom_actual;
    count = count + 1;
end
if count == count_max
   disp("Warning: Initial deformed geometries did not converge") 
end

% %Uncomment to look at difference in prescribed and actual initial tree
% %geometries
% ylabels = {'Radius (mm)', 'Radius Error (-)'};
% scaling = [1D3];
% figure(15)
% subplot(1, 2, 1)
% semilogy(orders, order_geom(1, :) * scaling(1), 'o','LineWidth', 2.0, 'Color', [0.0 0.0 0.0], 'MarkerSize',8)
% hold on
% semilogy(orders, order_geom_actual(1, :) * scaling(1), '*', 'LineWidth', 2.0, 'Color', [0.5 0.5 0.5], 'MarkerSize',8)
% xlabel('Vessel Order (-)'); ylabel(ylabels{1});
% xlim([0 11])
% xticks(1:1:10)
% ylim([10^-2 1])
% set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
% 
% subplot(1, 2, 2)
% plot(orders, (order_geom_actual(1, :) - order_geom(1, :)) ./ order_geom(1, :), 'o', 'LineWidth', 2.0, 'Color', [0.0 0.0 0.0], 'MarkerSize',8)
% xlabel('Vessel Order (-)'); ylabel(ylabels{2});
% axis([0 11 -0.05 0.05])
% xticks(1:1:10)
% yticks(-0.05:0.025:0.05)
% set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
% set(gcf,'Units','inches','Position',[.25 .25 12 4],'PaperUnits','inches','PaperPosition',[.25 .25 12 4])
% 
% fig_num = 16;
% hemo_all = [hemo_prescr hemo];
% fig16 = plot_hemo_mult_init(hemo_all, n_orders, fig_num);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Performing G&R simulations to adapt tree to perturbations
%--------------------------------------------------------------------------
%Run initial G&R to reset outputs and establish history
cmd_flags = "-r 0 -m " + max_days + " -s " + save_steps + " -d " + step_size + " -u 1 --simulate_equil 0";
output_flags = 2;
order_geom = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);

diameter = 2 * order_geom(1, :) * 1D2; %convert to cm
thickness = order_geom(2, :) * 1D2; %convert to cm

%Add in pertubation to start G&R
Q_LPA_new = (1.0 + dQ_s(1)) * Q_LPA; %ml / s 
P_term_new = (1.0 + dP_s(1)) * P_term; %mmHg

order = n_orders;
order_remodeling = n_orders;
alpha = 1;
beta = 1;
gen = 0;
resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
hemo_perturb_rigid = calculate_hemo_morphometric_tree_termBC(Q_LPA_new, P_term_new);
hemo_perturb = hemo_perturb_rigid;

gamma_p_ord = (hemo_perturb_rigid.P_order(:,1) ./ hemo.P_order(:,1) - 1.0);
gamma_q_ord = (hemo_perturb_rigid.Q_order(:,1) ./ hemo.Q_order(:,1) - 1.0);


if single_vessels_flag == 0
    %Use GnR code to update tree geometries with elastic deformation only
    count = 1;
    diff = ones(count_max + 1,1);
    diff_max = 0.005; %tolerance for pressure diff with updates

    while (diff(count) > diff_max && count < count_max)

        cmd_flags = "-r 0 -u 1 -m 2 -s 2 -d 1 --simulate_equil 0";
        output_flags = 3;
        order_geom_perturb = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);

        diameter = 2 * order_geom_perturb(1, :) * 1D2; %convert to cm
        thickness = order_geom_perturb(2, :) * 1D2; %convert to cm

        %Recalculate hemodynamics for given orders
        order = n_orders;
        order_remodeling = n_orders;
        alpha = 1;
        beta = 1;
        gen = 0;
        resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
        hemo_temp = calculate_hemo_morphometric_tree_termBC(Q_LPA_new, P_term_new);

        diff(count + 1) = max(abs((hemo_temp.P_order(:,1) - hemo_perturb.P_order(:,1))./ hemo_perturb.P_order(:,1)));
        hemo_perturb = hemo_temp;
        gamma_p_ord = hemo_perturb.P_order(:,1) ./ hemo.P_order(:,1) - 1.0;
        gamma_q_ord = hemo_perturb.Q_order(:,1) ./ hemo.Q_order(:,1) - 1.0;

        count = count + 1;
    end
    if count == count_max
       disp("Warning: Elastic perturbation adaptation did not converge") 
    end

    %Order commands for debugging
    % -r 0 -m max_days -s save_days -u 1 --simulate_equil 0 --gamma_p gamma_p_ord(ord) --gamma_q gamma_q_ord(ord) -n ord + ord
    tic
    %Step through gnr time course
    for t = 1:(max_days / (save_steps * step_size) - 1)
        disp("G&R time: " + t * save_steps + " days")
        Q_LPA_new = (1.0 + dQ_s(t + 1)) * Q_LPA; %ml / s 
        P_term_new = (1.0 + dP_s(t + 1)) * P_term; %mmHg
        
        %Iteratively update hemo at current time step
        count = 1;
        diff = ones(count_max + 1,1);
        diff_max = 0.001; %tolerance for pressure diff with updates
        hemo_new = hemo_perturb;
        while (diff(count) > diff_max && count < count_max)

            %Run gnr initially to establish history
            if count == 1
               gnr_iter_flag = " --gnr_iter_flag 0";
               gnr_out_flag = " --gnr_out_flag 0";
               output_flags = 4;
            else
               gnr_iter_flag = " --gnr_iter_flag 1";
               gnr_out_flag = " --gnr_out_flag 0";
               output_flags = 4;
            end
            cmd_flags = "-r 1 -m " + max_days + " -s " + save_steps + " -d " + step_size + " -u 1 --simulate_equil 0" + gnr_iter_flag + gnr_out_flag + mech_infl_flag;
            order_geom_new = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);

            diameter = 2 * order_geom_new(1, :) * 1D2; %convert to cm
            thickness = order_geom_new(2, :) * 1D2; %convert to cm

            %Recalculate hemodynamics for all orders
            order = n_orders;
            order_remodeling = n_orders;
            alpha = 1;
            beta = 1;
            gen = 0;
            resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
            hemo_temp = calculate_hemo_morphometric_tree_termBC(Q_LPA_new, P_term_new);

            diff(count + 1) = max(abs((hemo_temp.P_order(:,1) - hemo_new.P_order(:,1))./ hemo_new.P_order(:,1)));
            hemo_new = hemo_temp;
            gamma_p_ord = hemo_new.P_order(:,1) ./ hemo.P_order(:,1) - 1.0;
            gamma_q_ord = hemo_new.Q_order(:,1) ./ hemo.Q_order(:,1) - 1.0;

            count = count + 1;
        end

        %Run gnr initially to establish history
        gnr_iter_flag = " --gnr_iter_flag 1";
        gnr_out_flag = " --gnr_out_flag 1";
        cmd_flags = "-r 1 -m " + max_days + " -s " + save_steps + " -d " + step_size + " -u 1 --simulate_equil 0" + gnr_iter_flag + gnr_out_flag + mech_infl_flag;
        output_flags = 2;
        order_geom_new = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);

        diameter = 2 * order_geom_new(1, :) * 1D2; %convert to cm
        thickness = order_geom_new(2, :) * 1D2; %convert to cm

        %Recalculate hemodynamics for all orders
        order = n_orders;
        order_remodeling = n_orders;
        alpha = 1;
        beta = 1;
        gen = 0;
        resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
        hemo_temp = calculate_hemo_morphometric_tree_termBC(Q_LPA_new, P_term_new);

        hemo_new = hemo_temp;
        gamma_p_ord = hemo_new.P_order(:,1) ./ hemo.P_order(:,1) - 1.0;
        gamma_q_ord = hemo_new.Q_order(:,1) ./ hemo.Q_order(:,1) - 1.0;

        if t == (max_days / (save_steps * step_size) - 1)
            n_iter = 11;
            resistance_decr_iter = zeros(n_iter, 1);
            PVR_clinical_act = zeros(n_iter, 1);
            order_geom_active = order_geom_actual;
            gamma_act_ord = 0 * gamma_p_ord;
            gamma_p_ord_act = zeros(n_orders, n_iter);
            gamma_q_ord_act = zeros(n_orders, n_iter);
            hemo_act = hemo_new;

%             for iter = 1:n_iter
% 
%                 %set active stress for current order to zero
%                 gamma_act_ord(:) = gamma_act_ord(:) - (iter > 1) * 1 / (n_iter - 1);
% 
%                 %Use GnR code to update tree geometries with elastic deformation due to
%                 %active stress change only
%                 count = 1;
%                 diff = ones(count_max + 1,1);
%                 diff_max = 0.001; %tolerance for pressure diff with updates
% 
%                 while (diff(count) > diff_max && count < count_max)
% 
%                     %Run gnr initially to establish history
%                     gnr_iter_flag = " --gnr_iter_flag 1";
%                     gnr_out_flag = " --gnr_out_flag 0";
%                     mech_exp_flag = " --mech_exp_flag 1";
%                     cmd_flags = "-r 1 -m " + max_days + " -s " + save_steps + " -d " + step_size + " -u 1 --simulate_equil 0" + gnr_iter_flag + gnr_out_flag + mech_infl_flag + mech_exp_flag;
%                     output_flags = 4;
%                     order_geom_active = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);
% 
%                     diameter = 2 * order_geom_active(1, :) * 1D2; %convert to cm
%                     thickness = order_geom_active(2, :) * 1D2; %convert to cm
%                     
%                     %Recalculate hemodynamics for all orders
%                     order = n_orders;
%                     order_remodeling = n_orders;
%                     alpha = 1;
%                     beta = 1;
%                     gen = 0;
%                     resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
%                     hemo_temp = calculate_hemo_morphometric_tree_termBC(Q_LPA_new, P_term_new);
% 
%                     diff(count + 1) = max(abs((hemo_temp.P_order(:,1) - hemo_act.P_order(:,1))./ hemo_act.P_order(:,1)));
%                     hemo_act = hemo_temp;
%                     gamma_p_ord = hemo_act.P_order(:,1) ./ hemo.P_order(:,1) - 1.0;
%                     gamma_q_ord = hemo_act.Q_order(:,1) ./ hemo.Q_order(:,1) - 1.0;
%                     count = count + 1;
%                 end
%                 if count == count_max
%                    disp(iter) 
%                    disp("Warning: Elastic perturbation adaptation did not converge") 
%                 end
% 
%                 PVR_clinical_act(iter) = (hemo_act.P_order(n_orders,1) - hemo_act.P_order(1,1)) / Q_LPA_new;
%                 gamma_p_ord_act(:,iter) = gamma_p_ord;
%                 gamma_q_ord_act(:,iter) = gamma_q_ord;
%                 
%             end

        end

    end
    
    gnrtimer = toc;
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %Saving outputs from each set of modeling
    %--------------------------------------------------------------------------
        fig_num = 1;
        fig1 = plot_hemo(hemo, n_orders, fig_num);
        save_fig1 = "baseline_hemo" + ".png";
        saveas(fig1, save_fig1)

        fig_num = 2;
        fig2 = plot_hemo(hemo_perturb, n_orders, fig_num);
        save_fig2 = "perturb_hemo" + ".png";
        saveas(fig2, save_fig2)

        fig_num = 3;
        fig3 = plot_hemo(hemo_new, n_orders, fig_num);
        save_fig3 = "evolved_hemo" + ".png";
        saveas(fig3, save_fig3)

        fig_num = 4;
        hemo_all = [hemo hemo_perturb_rigid hemo_perturb hemo_new];
        fig4 = plot_hemo_mult(hemo_all, n_orders, fig_num);
        save_fig4 = "compare_hemo" + ".png";
        saveas(fig4, save_fig4)

        % fig_num = 4;
        % hemo_all = [hemo hemo_perturb_rigid hemo_perturb];
        % fig4 = plot_hemo_mult(hemo_all, n_orders, fig_num);
        % save_fig4 = "compare_hemo" + ".png";
        % saveas(fig4, save_fig4)

        fig_num = 5;
        fig5 = plot_geom_fold(order_geom_actual, order_geom_new, n_orders, fig_num);
        save_fig5 = "fold_geom_change" + ".png";
        saveas(fig5, save_fig5)

        fig_num = 6;
        fig6 = plot_compare_gnr(max_days, n_orders, step_size, fig_num);
        save_fig6 = "gnr_compare" + ".png";
        saveas(fig6, save_fig6)

        fig_num = 7;
        fig7 = plot_act_PVR_GnR(PVR_clinical_act, fig_num);
        save_fig7 = "pvr_gnr" + ".png";
        saveas(fig7, save_fig7)

        mkdir(save_dir)
        copyfile("morphometric_tree.mat", save_dir);

        save_figs = {save_fig1, save_fig2, save_fig3, save_fig4, save_fig5, save_fig6, save_fig7};
        num_save_figs = size(save_figs);
        for fs = 1:num_save_figs(2)
            copyfile(save_figs{fs}, save_dir);
        end
        for ord = 1:n_orders
            copyfile("GnR_out_ord" + ord, save_dir);
            copyfile("Native_in_ord" + ord, save_dir);
            copyfile("Vs_out_ord" + ord, save_dir);
        end

        if tree_solve_prev_flag == 0
            save("mat_tree_props.mat", "mat_tree");
        end
        copyfile("mat_tree_props.mat", save_dir);

        hemo_save.hemo_prescr = hemo_prescr;
        hemo_save.hemo = hemo;
        hemo_save.hemo_perturb_rigid = hemo_perturb_rigid;
        hemo_save.hemo_perturb = hemo_perturb;
        hemo_save.hemo_new = hemo_new;
        hemo_save.hemo_active = hemo_act;
        hemo_save.PVR_clinical_act = PVR_clinical_act;

        geom_save.order_geom_prescr = order_geom_prescr;
        geom_save.order_geom = order_geom;
        geom_save.order_geom_perturb_rigid = order_geom;
        geom_save.order_geom_perturb = order_geom_perturb;
        geom_save.order_geom_new = order_geom_new;
        geom_save.order_geom_active = order_geom_active;
        geom_save.gnrtime = gnrtimer;

        save("hemo_and_geom.mat","geom_save","hemo_save")
        copyfile("hemo_and_geom.mat", save_dir);
        %--------------------------------------------------------------------------

else
    %Run G&R for single vessels only
    nts_remain = (max_days / step_size) - 1;
    gnr_iter_flag = " --gnr_iter_flag 0";
    gnr_out_flag = " --gnr_out_flag 1";
    cmd_flags = "-r 1 -m " + max_days + " -s " + nts_remain + " -d " + step_size + " -u 1 --simulate_equil 0" + gnr_iter_flag + gnr_out_flag + mech_infl_flag;
    output_flags = 2;
    order_geom_new = calc_tree_equil_wgnr(order_geom_actual, mat_tree, mat, n_orders, hemo, gamma_p_ord, gamma_q_ord, gamma_act_ord, cmd_flags, output_flags);
    
    fig_num = 5;
    fig5 = plot_geom_fold(order_geom_actual, order_geom_new, n_orders, fig_num);
    save_fig5 = "fold_geom_change" + ".png";
    saveas(fig5, save_fig5)
    
    fig_num = 6;
    fig6 = plot_compare_gnr(max_days, n_orders, step_size, fig_num);
    save_fig6 = "gnr_compare" + ".png";
    saveas(fig6, save_fig6)
    
    mkdir(save_dir)
    copyfile("morphometric_tree.mat", save_dir);

    save_figs = {save_fig5, save_fig6};
    num_save_figs = size(save_figs);
    for fs = 1:num_save_figs(2)
        copyfile(save_figs{fs}, save_dir);
    end
    for ord = 1:n_orders
        copyfile("GnR_out_ord" + ord, save_dir);
        copyfile("Native_in_ord" + ord, save_dir);
        copyfile("Vs_out_ord" + ord, save_dir);
    end

    if tree_solve_prev_flag == 0
        save("mat_tree_props.mat", "mat_tree");
        copyfile("mat_tree_props.mat", save_dir);
    end
    
end


