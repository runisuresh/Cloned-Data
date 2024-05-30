function write_Native_in(mat, ord, ri_iv, h_iv, P_iv, Q_iv, output_flags)
%Writes input file for growth and remodeling program

    vessel_name = strcat('PAgen_', num2str(ord));
    
    %some default parameter values than need to be written for G&R input
    lambda_z_iv = 1.0; %default in vivo axial stretch
    alpha_mechinfl = zeros(1, 6, 'uint32'); %infl for homeostatic constituent type
    beta_e = -1; %eln alignment
    beta_m = 90; %SMC alignment
    ge_aniso = 1.0; %eln anisotropic pre-stretch
    rhoR_h = 1050; %true mass density
    k_e_exp = 0.0;
    
    %default mechanobiological parameters
    eta_K = 4; %4 deft
    eta_ups = 1.0;
    K_c_tauw = 1.00; %1 deft
    K_c_sigma = eta_K * K_c_tauw;
    K_m_tauw = eta_ups * K_c_tauw;
    K_m_sigma = eta_ups * K_c_sigma;
    
    K_sigma_p = [1 K_m_sigma K_c_sigma K_c_sigma K_c_sigma K_c_sigma];
    K_sigma_d = ones(1,6);
    K_tauw_p = [1 K_m_tauw K_c_tauw K_c_tauw K_c_tauw K_c_tauw];
    K_tauw_d = zeros(1,6);
    
    %default active stress parameters
    lambda_m = 1.10;
    lambda_0 = 0.40;
    k_act = 1.0 / 7.0;
    CB = 0.8326;
    CS = 0.5 * 0.8326;
    
    %default elastin degradation parameters (i.e. none)
    s_edeg_off = 0;
    epsilonR_e_min = 1.00;
    k_e_h = 1 / 14;
    
    %default inflammation parameters (i.e. none)
    Ki_trans = 0.0;
    beta_i = 0.0;
    Ki_deg = 0.0;
    
    Ki_steady = 0.0;
    delta_i = 1.0;
    gamma_i = 0.00;
    Ki_Tact = 1.0;
    phi_Tact0_min = 0.00;
    K_infl_eff = 0.00;
    K_mech_eff = 1.00;
    s_int_infl = 0;
    delta_m = 1.0;
    s_int_mech = 0;
    
    %set constituents to experience infl
    alpha_mechinfl = [0 1 0 0 0 0];
    
%     %Varying infl
%     %linear: * (1 -  (ord - 1) / (10 - 1));
%     %exponential: * (1 - exp(ord - 10));
%     Ki = 0.0;
%     gammai = 16.0;
%     Kme = 1.00;
%     phiTactmin = 1.0;
%     elnmin = 1.0;
%     
%     istim = 0.00 * Ki * (1 - exp(ord - 10)) + 1.00 * Ki; %
%     gammastim = gammai; %0.8 * gammai * (1 - exp(ord - 10)) + 0.2 * gammai;
%     mechstim = Kme; %(1 - Kme) * exp(ord - 10) + Kme;
%     %elnstim = (1 - elnmin) * exp(ord - 10) + elnmin;
%     %elnstim = 1 - (1 - elnmin) * exp(ord - 10);
%     %phiTactminstim = phiTactmin * (1 - exp(ord - 10));
%     
%     Ki_steady = istim;
%     gamma_i = gammastim;
%     K_mech_eff = mechstim;
%     phi_Tact0_min = phiTactmin;
%     Ki_Tact = 1 - phiTactmin;
%     epsilonR_e_min = elnmin;
%     
%     delta_i = 0.5;
%     K_infl_eff = 1.00;
%     s_int_infl = 0;
% 
%     delta_m = 1.0;
%     s_int_mech = 0;

%     %Uniform infl
    Ki_Tact = 0.0;
    phi_Tact0_min = 1.00;
    gamma_i = 32;
    Ki_steady = 1.0; %0.80;
    delta_i = 1.00;
    K_infl_eff = 1.00;
    s_int_infl = 0;

    delta_m = 1.00;
    K_mech_eff = 0.25;
    s_int_mech = 0;

    epsilonR_e_min = 1.00;

%     WSS infl
%     Ki_Tact = 0.0;
%     phi_Tact0_min = 1.0;
%     gamma_i = 8;
%     Ki_steady = 0.00; %0.80;
%     delta_i = 1.0;
%     K_infl_eff = 0.08;
%     s_int_infl = 0;
% 
%     delta_m = 1.0;
%     K_mech_eff = 0.25;
%     s_int_mech = 0;
% 
%     epsilonR_e_min = 0.5;

    %toggle single vessel vs vessel order outputs
    if output_flags == 2
        filename = "Native_in_ord" + ord;
    else
        filename = 'Native_in_';
    end
    
    %toggle mass change for evaluation of initial behavior
    if (output_flags == 3)
        kdeg = 1.0 / 1000000;    
    else
    	kdeg = 1.0 / 7;
    end
    
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', vessel_name);
    fprintf(fid, '%4.10f %4.10f\n', ri_iv, h_iv);
    fprintf(fid, '%4.10f\n', lambda_z_iv);
    fprintf(fid, '%d %d %d %d %d %d\n', alpha_mechinfl);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ...
                 mat.ce, 0, mat.cm, mat.ck);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', beta_e, beta_m, mat.beta*180/pi);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', mat.ge);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', ge_aniso, mat.gm, mat.gc, mat.gc, mat.gc, mat.gc);
    fprintf(fid, '%4.10f\n', rhoR_h);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', mat.phie, mat.phim, mat.phic * mat.phik(1), mat.phic * mat.phik(2), mat.phic * mat.phik(3), mat.phic * mat.phik(4));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', k_e_exp, kdeg * ones(1,5));
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', K_sigma_p);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', K_sigma_d);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', K_tauw_p);
    fprintf(fid, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n', K_tauw_d);
    fprintf(fid, '%4.10f %4.10f\n', P_iv, Q_iv);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', k_act, lambda_0, lambda_m);
    fprintf(fid, '%4.10f %4.10f\n', CB, CS);
    fprintf(fid, '%4.10f\n', mat.T_act);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', s_edeg_off, epsilonR_e_min, k_e_h);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', Ki_trans, Ki_steady, Ki_deg);
    fprintf(fid, '%4.10f %4.10f\n', delta_i, beta_i);
    fprintf(fid, '%4.10f %4.10f\n', K_infl_eff, s_int_infl);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', gamma_i, Ki_Tact, phi_Tact0_min);
    fprintf(fid, '%4.10f %4.10f %4.10f\n', delta_m, K_mech_eff, s_int_mech);
    
    fclose(fid);
    
end

