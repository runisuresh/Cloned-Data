
days_run = 40;

params_uniform_infl.Ki_Tact = 0.0;
params_uniform_infl.phi_Tact0_min = 1.00;
params_uniform_infl.gamma_i = 32;
params_uniform_infl.Ki_steady = 1.0; %0.80;
params_uniform_infl.delta_i = 1.00;
params_uniform_infl.K_infl_eff = 1.00;
params_uniform_infl.s_int_infl = 0;
params_uniform_infl.delta_m = 1.00;
params_uniform_infl.K_mech_eff = 0.25;
params_uniform_infl.s_int_mech = 0;
params_uniform_infl.epsilonR_e_min = 1.00;

simulation_settings.tree_diamrat_prev_flag = 0; %Use diameter ratio from a previous solve 1
simulation_settings.tree_solve_prev_flag = 0; %Use down the tree material parameters from a previous solve 1
simulation_settings.tree_gen_prev_flag = 0; %Use a morphometric tree from a previous solve 1
simulation_settings.single_vessels_flag = 0; %Disable hemodynamic feedback 1
simulation_settings.mech_infl_flag = " --mech_infl_flag 0"; %Use inflammation driven by WSS 1
simulation_settings.dQ = 0.00; %Pertubation in inflow
simulation_settings.dP = 0.00; %Pertubation in outlet pressure
simulation_settings.k_ramp = 1/4; %days to ramp hemodynamic changes over
simulation_settings.step_size = 1.0; %size of step in days
simulation_settings.max_days = days_run; %maximum days to run the G&R
simulation_settings.save_steps = 1; %number of steps between hemo->G&R feedback, min 1
simulation_settings.other = "test_run2";
simulation_settings.loop_iteration = 1; % Initialize as needed

ts = 0: simulation_settings.step_size: simulation_settings.max_days / simulation_settings.step_size - 1;
simulation_settings.dQ_s = simulation_settings.dQ * (1 - exp( -simulation_settings.k_ramp * ts));
%simulation_settings.dQ_s = dQ_inputs;

[~, ~] = fun_run_tree_GnR(simulation_settings, params_uniform_infl);