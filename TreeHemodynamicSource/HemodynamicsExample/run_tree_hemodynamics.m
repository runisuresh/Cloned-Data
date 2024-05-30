close all; clear all; clc; clf;

hemo_timer = tic;

global connectivity diameter thickness length new_conn

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

new_conn = 1;

%Initialize vessel order metrics
n_orders = size(connectivity, 1); %Number of vessel orders, vessels decrease in size with decreasing order
orders = 1:n_orders;

% %--------------------------------------------------------------------------
% %Uncomment for Previous Mouse tree generation
% %--------------------------------------------------------------------------
% 
% %Biomechanical data for mouse LPA, from Ramachandra and Humphrey 2019, J Biomech
% load('Ramachandra_LPA_fit_mm4ff.mat');
% 
% %Extract RPA geometry from structure;
% r_iv_LPA = kinem.ri_iv_act;
% h_iv_LPA = kinem.h_iv_act;
% l_iv_LPA = kinem.F_iv(3) * kinem.L;
% V_LPA = kinem.V;
% 
% %Prescribed hemodynamics from the literature 
% CO = 10.4; %ml/min, From Champion et al., 2000 Am J Heart Circ Physiol
% FlowSplit = 0.30; %From Razavi et al., 2011 Am J Physiol Lung Cell Mol Physiol
% Q_LPA = CO / 60 * FlowSplit; %ml/s 
% mPAP = kinem.P_iv / 133.33; %mmHg, From Champion et al., 2000 Am J Heart Circ Physiol
% PCWP =  4.8; %mmHg, From Champion et al., 2000 Am J Heart Circ Physiol
% PVR_LPA = (mPAP - PCWP) / Q_LPA * 1333.33;
% P_term = PCWP; %mmHg, terminal pressure of the tree
% 
% input_pars.n_orders = n_orders;
% input_pars.r_iv_LPA = r_iv_LPA;
% input_pars.h_iv_LPA = h_iv_LPA;
% input_pars.l_iv_LPA = l_iv_LPA;
% input_pars.V_LPA = V_LPA;
% input_pars.PVR_LPA = PVR_LPA;
% input_pars.len_ratio = 1.60;
% 
% lb = 1.0;
% ub = 2.0;
% par0 = 1.4;
% options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-12,'TolX',1e-12,...
%                        'Display','iter');
% current_objective = @(par)ord_diam_obj(par, input_pars);
% diam_ratio_const = lsqnonlin(current_objective, par0, lb, ub, options);
% 
% %Diam ratio is diam order n / n - 1 -> d_n-1 = d_n / d_ratio
% %Thick ratio is thick order n / diam order n -> h_n-1 = d_n * h_ratio
% %Len ratio is len order n / n - 1 -> l_n-1 = l_n / l_ratio
% order_diam_ratio = diam_ratio_const * ones(n_orders - 1,1); %can vary along the orders, but constant for now, reasonably linear in Jiang
% thick_ratio_const = h_iv_LPA / (2 * r_iv_LPA); %calculated from in vivo gemoetry from data in Ramachandra and Humphrey 2019, J Biomech
% order_thick_ratio = thick_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, some evidence of constant in humans after infancy, see Figure 1 Haworth et al., 1983 Am J Cardiol
% thick_rat_incr = 1.00; %increasing thickness ratio for each order
% len_ratio_const = 1.60; %from regression in Figure 5, Jiang et al., 1994 J App Physio
% order_len_ratio = len_ratio_const * ones(n_orders - 1, 1); %can vary along the orders, but constant for now, reasonably linear in Jiang
% 
% %Create radius, thickness, length, volume array for each order 
% order_geom = zeros(4, n_orders);
% 
% %Initialize highest order values
% order_geom(1, n_orders) = r_iv_LPA;
% order_geom(2, n_orders) = h_iv_LPA;
% order_geom(3, n_orders) = l_iv_LPA;
% order_geom(4, n_orders) = V_LPA;
% %Loop through lower orders
% for ord = n_orders - 1: -1: 1
%     order_geom(1, ord) = order_geom(1, ord + 1) / order_diam_ratio(ord);
%     order_geom(2, ord) = order_thick_ratio(ord) * 2 * order_geom(1, ord) * thick_rat_incr^(n_orders - ord);
%     order_geom(3, ord) = order_geom(3, ord + 1) / order_len_ratio(ord);
%     order_geom(4, ord) = pi * ( (order_geom(1, ord) + order_geom(2, ord))^2 - order_geom(1, ord)^2 ) * order_geom(3, ord);
% end
% order_geom_prescr = order_geom;
% 
% diameter = 2 * order_geom(1, :) * 1D2; %convert to cm
% thickness = order_geom(2, :) * 1D2; %convert to cm
% length = order_geom(3, :) * 1D2; %convert to cm
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Uncomment for Rat tree generation
%--------------------------------------------------------------------------
diameter_LPA_rat = sqrt(2.6 / pi) / 10 * 2; % from 2.6mm^2 CSA to diam from Razavi et al., 2011 Am J Physiol Lung Cell Mol Physiol
h_LPA_rat = 1.0;
l_iv_LPA_rat = 4.21 / 10; % from Purcell et al., 2016 Am J Physiol Lung Cell Mol Physiol

diam_ratio_const = 1.58;
order_diam_ratio = diam_ratio_const * ones(n_orders - 1,1);
len_ratio_const = 1.60;
order_len_ratio = len_ratio_const * ones(n_orders - 1, 1);

%Create radius, thickness, length, volume array for each order 
order_geom = zeros(3, n_orders);
order_geom(1, n_orders) = diameter_LPA_rat;
order_geom(2, n_orders) = 1;
order_geom(3, n_orders) = l_iv_LPA_rat;

%Initialize highest order values
for ord = n_orders - 1: -1: 1
    order_geom(1, ord) = order_geom(1, ord + 1) / order_diam_ratio(ord);
    order_geom(2, ord) = 1;
    order_geom(3, ord) = order_geom(3, ord + 1) / order_len_ratio(ord);
end

% diameter = [13.3; 31.7; 44.4; 61.5; 88.1; 152; 266; 417; 602; 929; 1639;] * 1D-4; %convert to cm
%length = [0.05; 0.15; 0.20; 0.27; 0.40; 0.72; 1.24; 1.33; 1.74; 2.64; 18.11] * 1D-1; %convert to cm
% thickness = ones(n_orders,1);

diameter = order_geom(1, :); %convert to cm
thickness = order_geom(2, :); %convert to cm
length = order_geom(3, :); %convert to cm

%Prescribed hemodynamics from the literature 
CO = 125; %ml/min, %From Razavi et al., 2011 Am J Physiol Lung Cell Mol Physiol
FlowSplit = 0.30; %From Razavi et al., 2011 Am J Physiol Lung Cell Mol Physiol
Q_LPA = CO / 60 * FlowSplit; %ml/s
mPAP = 16; %mmHg, From van der Feen et al., 2017 JOVE
PCWP =  1.6; %mmHg, From van der Feen et al., 2017 JOVE
PVR_LPA = (mPAP - PCWP) / Q_LPA * 1333.33;
P_term = PCWP; %mmHg, terminal pressure of the tree

%--------------------------------------------------------------------------

%Calculate hemodynamics for given orders
order = n_orders;
order_remodeling = n_orders;
alpha = 1;
beta = 1;
gen = 1;
resistance = find_opt_morphometric_tree(order, order_remodeling, alpha, beta, gen);
hemo = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Saving outputs from each set of modeling
%--------------------------------------------------------------------------
    
fig_num = 1;
fig1 = plot_hemo(hemo, n_orders, fig_num);
save_fig1 = "baseline_hemo" + ".png";
saveas(fig1, save_fig1)

%--------------------------------------------------------------------------

hemo_timer = toc;

