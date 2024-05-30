clc; clear all; close all;

%Set options for run
test_ind = 1:3;
par_est_flag = 0;
par_wr_flag = 0;
plot_flag = 1;

%Refitting data from Ramachandra and Humphrey, J Biomech, 2019 with a
%microstructurally motivated model
%Uses an in vivo reference configuration with active deformation assumed,
%is pulled back to the traction-free config., and then pushed into the
%loaded passive configurations from the biaxial mechanical testing protocol
%for parameter estimation

%load mechanical testing data
load('Ramachandra_LPA_exp_data.mat');

%IV geometry assumed from 15 mmHg
%Could have this read in from generate_data file
lz_iv_pas = lz_iv;
lt_iv_pas = (ri_iv + h_iv / 2) / (Ri + H / 2);
lr_iv_pas = 1 / (lt_iv_pas * lz_iv_pas);
F_iv_pas = [lr_iv_pas; lt_iv_pas; lz_iv_pas];

%From Figure 3B
lz_iv_act = 1.00;
lt_iv_act = 0.90;
lr_iv_act = 1 / (lz_iv_act * lt_iv_act);
F_iv_act = [lr_iv_act; lt_iv_act; lz_iv_act];

%Combined deformations to traction-free to in-vivo state
F_iv =  F_iv_pas .* F_iv_act;

if par_est_flag == 1
    %Material Parameters
    %Eln, circ, axial, diag1, diag2
    mat.ge = F_iv;
    mat.gm = (F_iv(2) + F_iv(3)) / 2;
    mat.gc = (F_iv(2) + F_iv(3)) / 2;
    mat.phik = [0.25 0.25 0.25 0.25];
    mat.beta = [90 0 35.7797 360 - 35.7797] * pi/180;
    mat.phie = 0.255; %elastin
    mat.phim = 0.180; %smooth muscle
    mat.phic = 1 - mat.phie - mat.phim; % collagen
    mat.T_act = 00000;
    mat.l0 = 0.4;
    mat.lm = 1.1;
    mat.CB = 0.8326;
    mat.CS = mat.CB / 2;
    
    par0 = [rand(1,5) ...
            rand*pi/180 2*rand(1,4) 0.5*rand(1,2)];
        
    lb = [0 0 0 0 0 ...
          0 1.05 1.05 1.05 1.05 0 0];
      
    ub = [1e+7 1e+7 1e+7 1e+7 1e+7 ...
         pi/2 2.00 2.00 2.0 2.0 0.5 0.5];

    % Minimization of the objective function
    options = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-12,'TolX',1e-12,...
                       'Display','iter');
    current_objective = @(par)par_est_obj(par, mat, test_ind, exp_data, Ri, H, V, L, F_iv);
    par_est = lsqnonlin(current_objective, par0, lb, ub, options);

    mat.ce = par_est(1);
    mat.cm = par_est(2:3);
    mat.ck = [par_est(4:5) par_est(4:5) par_est(4:5) par_est(4:5)];
    mat.beta(3) = par_est(6);
    mat.beta(4) = 2 * pi - par_est(6);
    mat.ge = [1 / (par_est(7) * par_est(8)); par_est(7); par_est(8)];
    mat.gm = par_est(9);
    mat.gc = par_est(10);
    mat.phik = [par_est(11) par_est(12) (1 - par_est(11) - par_est(12))/2 (1 - par_est(11) - par_est(12))/2];
    %mat.phik = [0 par_est(11) (1 - par_est(11))/2 (1 - par_est(11))/2];

%     save('Ramachandra_LPA_fit_mm4ff.mat','mat');

else
    load('Ramachandra_LPA_fit_mm4ff.mat');
end

sim_data = 0 * exp_data;
ndata = length(squeeze(exp_data(:,1,1))) * length(test_ind);
all_P_err = zeros(ndata,1);
all_L_err = zeros(ndata,1);
datan = 1;
for testn = test_ind

    %loading inputs
    P_range = squeeze(exp_data(:,1,testn));
    lz_range = squeeze(exp_data(:,5,testn));
    mean_P = mean(P_range);

    %experimental comparisons
    ro_range = squeeze(exp_data(:,2,testn));
    ri_range = 0 * squeeze(exp_data(:,2,testn));
    L_range = squeeze(exp_data(:,3,testn));
    mean_L = mean(L_range);

    %number of points to simulate
    nsim = length(P_range);

    %simulated comparisons
    P_range_sim = zeros(nsim,1);
    L_range_sim = zeros(nsim,1);
    sigmat_range_sim = zeros(nsim,1);
    sigmaz_range_sim = zeros(nsim,1);
    
    for i = 1:nsim

        ro = ro_range(i);
        lz = lz_range(i);

        ri = sqrt(ro^2 - V/(pi * lz * L));
        ri_range(i) = ri;
        t = calc_stress(ro, mat, Ri, H, V, L, lz, F_iv);
        P_range_sim(i) = t(2) * (ro - ri) / ri;
        L_range_sim(i) = calc_axial_force(ro, ri, t(3));
        sigmat_range_sim(i) = t(2);
        sigmaz_range_sim(i) = t(3);
        
        all_P_err(datan) = (P_range_sim(i) - P_range(i)) / mean_P;
        all_L_err(datan) = (L_range_sim(i) - L_range(i)) / mean_L;
        datan = datan + 1;

    end
    
    sim_data(:,1,testn) = P_range_sim;
    sim_data(:,2,testn) = ro_range;
    sim_data(:,3,testn) = L_range_sim;
    sim_data(:,4,testn) = exp_data(:,4,testn);
    sim_data(:,5,testn) = lz_range;
    sim_data(:,6,testn) = sigmat_range_sim;
    sim_data(:,7,testn) = sigmaz_range_sim;
    
end

err = sum(all_P_err.^2);

P_iv_orig = P_iv;
[temp, iv_ind] = min(abs(sim_data(:,1,2) - P_iv_orig));
P_iv = sim_data(iv_ind,1,2);
ro_iv_pas = sim_data(iv_ind,2,2);
ri_iv_pas = sqrt(ro_iv_pas^2 - V/(pi * lz_iv_pas * L));
h_iv_pas = ro_iv_pas - ri_iv_pas;

lz_iv_pas = lz_iv;
lt_iv_pas_new = (ri_iv_pas + h_iv_pas / 2) / (Ri + H / 2);
lr_iv_pas_new = 1 / (lt_iv_pas_new * lz_iv_pas);
F_iv_pas_new = [lr_iv_pas_new; lt_iv_pas_new; lz_iv_pas];

%Check if in vivo ref and traction free ref to in vivo give same stress state
%traction free ref
t_tfref = calc_stress(ro_iv_pas, mat, Ri, H, V, L, F_iv_pas_new(3), F_iv_pas_new);

%in vivo passive ref
lz_ivref = 1.0;
L_ivref = L * F_iv_pas_new(3);
F_ivref = ones(3,1);
t_ivref = calc_stress(ro_iv_pas, mat, ri_iv_pas, h_iv_pas, V, L_ivref, lz_ivref, F_ivref);

if ( abs((t_ivref(2) - t_tfref(2)) / t_tfref(2)) > 0.01)
    disp('Problem with Deformations - different stresses with tf reference and in vivo reference')
end

%Solve for active stress parameters that give true in vivo state
P = P_iv;
r_mid_pas = ri_iv_pas + h_iv_pas / 2;
r_mid_act = r_mid_pas * F_iv_act(2);
h_act = h_iv_pas * F_iv_act(1);
ro_act = r_mid_act + h_act / 2;
ri_act = ro_act - h_act;
T_act_g = 50000;
objective = @(T_act)active_obj(T_act, ro_act, mat, P_iv, ri_iv_pas, h_iv_pas, V, L_ivref, lz_ivref, F_iv_act);
T_act_solve = fzero(objective, T_act_g);
mat.T_act = T_act_solve;

kinem.Ro = Ro;
kinem.H = H;
kinem.Ri = Ri;
kinem.L = L;
kinem.V = V;
kinem.F_iv_pas = F_iv_pas_new;
kinem.F_iv_act = F_iv_act;
kinem.F_iv = F_iv;
kinem.ri_iv_pas = ri_iv_pas;
kinem.h_iv_pas = h_iv_pas;
kinem.ri_iv_act = ri_act;
kinem.h_iv_act = h_act;
kinem.P_iv = P_iv;

if (par_wr_flag == 1)
    save('Ramachandra_LPA_fit_mm4ff.mat','mat','kinem');
end

%Active stress perturbations
% objective = @(ro)loaded_obj(ro, mat, P_iv, ri_iv_pas, h_iv_pas, V, L_ivref, lz_ivref, F_iv_act);
% ro_g = ro_iv_pas*0.9;
% ro_act = fzero(objective, ro_g);
% ri_act = sqrt(ro_act^2 - V/(pi * lz_ivref * L_ivref));
% h_act = ro_act - ri_act;
% t_act = calc_stress(ro_act, mat, ri_act, h_act, V, L_ivref, lz_ivref, F_ivref);
% fprintf('%s %f %s %f \n', 'ri_act: ', ri_act, 'h_act: ', h_act);

%print new native vessel text input file for gnr code
if par_wr_flag == 1
    
    fid = fopen('Native_in.txt', 'w');

    fprintf(fid, '%4.10f %4.10f\n', ri_act, h_act);
    fprintf(fid, '%4.10f\n', 1.00);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ...
                 mat.ce, 0, mat.cm, mat.ck);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', -1, 90, 90, 0, mat.beta(3)*180/pi, mat.beta(4)*180/pi);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', mat.ge);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', 1.0, mat.gm, mat.gc, mat.gc, mat.gc, mat.gc);
    fprintf(fid, '%4.10f\n', 1050);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', mat.phie, mat.phim, mat.phic * mat.phik(1), mat.phic * mat.phik(2), mat.phic * mat.phik(3), mat.phic * mat.phik(4));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', 0, 1.0 / 80.0 * ones(1,5));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ones(1,6));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ones(1,6));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ones(1,6));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', zeros(1,6));
    fprintf(fid, '%4.10f %4.10f\n', P_iv, 20);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', 1.0 / 7.0, 0.4, 1.1);
    fprintf(fid, '%4.10f %4.10f\n', 0.8326, 0.8326 / 2);
    fprintf(fid, '%4.10f %4.10f\n', mat.T_act);

    fclose(fid);
end


if plot_flag == 1

    figure(1)
    subplot(2,1,1)
    hold on
    plot(2 * exp_data(:,2,1)* 1D6, exp_data(:,1,1)/133, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
    plot(2 * exp_data(:,2,2) * 1D6, exp_data(:,1,2)/133.33, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    plot(2 * exp_data(:,2,3) * 1D6, exp_data(:,1,3)/133.33, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
    plot(2 * sim_data(:,2,1)* 1D6, sim_data(:,1,1)/133, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(2 * sim_data(:,2,2) * 1D6, sim_data(:,1,2)/133.33, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(2 * sim_data(:,2,3) * 1D6, sim_data(:,1,3)/133.33, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    axis([600 1800 0 50])
    xlabel('Outer Diameter \mum'); ylabel('Pressure (mmHg)')
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
    subplot(2,1,2)
    hold on
    plot(exp_data(:,4,1), exp_data(:,3,1)*1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
    plot(exp_data(:,4,2), exp_data(:,3,2)*1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    plot(exp_data(:,4,3), exp_data(:,3,3)*1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
    plot(sim_data(:,4,1), sim_data(:,3,1)*1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,2), sim_data(:,3,2)*1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,3), sim_data(:,3,3)*1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    axis([1 1.8 0 30])
    xlabel('Circ. Stretch'); ylabel('Force (mN)')
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

    figure(2)
    subplot(2,1,1)
    hold on
    plot(exp_data(:,4,1), exp_data(:,6,1)/1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
    plot(exp_data(:,4,2), exp_data(:,6,2)/1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    plot(exp_data(:,4,3), exp_data(:,6,3)/1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
    plot(sim_data(:,4,1), sim_data(:,6,1)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,2), sim_data(:,6,2)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,3), sim_data(:,6,3)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    axis([0.5 2.5 0 160])
    xlabel('Circ. Stretch'); ylabel('Circ. Stress (kPa)')
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
    subplot(2,1,2)
    hold on
    plot(exp_data(:,4,1), exp_data(:,7,1)/1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
    plot(exp_data(:,4,2), exp_data(:,7,2)/1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    plot(exp_data(:,4,3), exp_data(:,7,3)/1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
    plot(sim_data(:,4,1), sim_data(:,7,1)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,2), sim_data(:,7,2)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    plot(sim_data(:,4,3), sim_data(:,7,3)/1D3, 'LineWidth', 2.0, 'Color','r', 'LineStyle', ':')
    axis([0.5 2.5 0 140])
    xlabel('Circ. Stretch'); ylabel('Axial Stress (kPa)')
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

end

%--------------------------------------------------------------------------
%Functions
function J = par_est_obj(par, mat, test_ind, exp_data, Ri, H, V, L, F_iv)

    mat.ce = par(1);
    mat.cm = par(2:3);
    mat.ck = [par(4:5) par(4:5) par(4:5) par(4:5)];
    mat.beta(3) = par(6);
    mat.beta(4) = 2 * pi - par(6);
    mat.ge = [1 / (par(7) * par(8)); par(7); par(8)];
    mat.gm = par(9);
    mat.gc = par(10);
    mat.phik = [par(11) par(12) (1 - par(11) - par(12))/2 (1 - par(11) - par(12))/2];
    %mat.phik = [0 par(11) (1 - par(11))/2 (1 - par(11))/2];

    
    ndata = length(squeeze(exp_data(:,1,1))) * length(test_ind);
    all_P_err = zeros(ndata,1);
    all_L_err = zeros(ndata,1);
    
    datan = 1;
    for testn = test_ind
        
        %loading inputs
        P_range = squeeze(exp_data(:,1,testn));
        lz_range = squeeze(exp_data(:,5,testn));
        mean_P = mean(P_range);
        
        %experimental comparisons
        ro_range = squeeze(exp_data(:,2,testn));
        L_range = squeeze(exp_data(:,3,testn));
        mean_L = mean(L_range);
        
        %number of points to simulate
        nsim = length(P_range);
       
        %simulated comparisons
        P_range_sim = zeros(nsim,1);
        L_range_sim = zeros(nsim,1);
        
        for i = 1:nsim
            
            ro = ro_range(i);
            lz = lz_range(i);
            
            ri = sqrt(ro^2 - V/(pi * lz * L));
            t = calc_stress(ro, mat, Ri, H, V, L, lz, F_iv);
            P_range_sim(i) = t(2) * (ro - ri) / ri;
            L_range_sim(i) = calc_axial_force(ro, ri, t(3));
            
            all_P_err(datan) = (P_range_sim(i) - P_range(i)) / mean_P;
            all_L_err(datan) = (L_range_sim(i) - L_range(i)) / mean_L;
            datan = datan + 1;
            
        end
        
    end
    
    J = [all_P_err all_L_err];

end

function J = active_obj(T_act, ro, mat, P, Ri, H, V, L, lz, F_iv)
    
    mat.T_act = T_act;

    objective = @(ro)loaded_obj(ro, mat, P, Ri, H, V, L, lz, F_iv);
    ro_g = ro;
    ro_act = fzero(objective, ro_g);
    
    J = ro_act - ro;

end


function L = calc_axial_force(ro, ri, tz)

    L = pi*(ro^2 - ri^2) * tz;

end


