clear; clc;

% path = "/Users/vaaru/OneDrive/Desktop/CGRL/Extracted_MCT_Data.xlsx";
path = fullfile('/mnt/c/Users/vaaru/OneDrive/Desktop/CGRL', 'Extracted_MCT_Data.xlsx');
sheets = sheetnames(path);
days_run = 40;
inputs = {'Cardiac output (ml/min)','Cardiac index'};
metric_names = "RVSP (mmHG)";
num_mets = numel(metric_names);
exp_data = cell(1,num_mets);
co_data = [];

sheetIndex = find(strcmp(sheets, 'LiuEtAl2016'));
curr_sheet = sheets(sheetIndex);
    
%Keep headers intact
opts = detectImportOptions( ...
            path, ...
            "VariableNamingRule","preserve", ...
            "Sheet",curr_sheet ...
            );
column_headers = opts.VariableNames;

sheet_metric = readtable(path,opts);
disp(sheet_metric)

%Normalize CI data by time-matched control (e.g. "Control Day 8")
if strcmp(curr_sheet,'Cardiac index')

%Convert to cardiac output (CO)
sheet_metric.CO = sheet_metric{:,2} .* sheet_metric{:,3};
sheet_metric(:,[2,3]) = [];
numrows = height(sheet_metric);
time_pts = 0.5 * numrows;
treatment_gp = sheet_metric(1:time_pts,'CO');
ctrl_gp = sheet_metric((time_pts+1):numrows,'CO');
norm_co = treatment_gp ./ ctrl_gp;
sheet_metric(1:time_pts,"CO") = norm_co;
end

%Normalize CO data to only the "Day 0" control
if strcmp(curr_sheet,'Cardiac output (ml/min)')
     sheet_metric = renamevars(sheet_metric,"Cardiac output (ml/min)","CO");
     ctrl = sheet_metric{1,2};
     numrows = height(sheet_metric);
    for row = 1:numrows
        sheet_metric{row,2} = sheet_metric{row,2} / ctrl;
    end
end

%Remove NaN values
sheet_metric = rmmissing(sheet_metric);

%If time is in weeks, change units to days
if strcmp(sheet_metric.Properties.VariableNames{1}, 'Time (wks)')
     sheet_metric{:,1} = 7 * sheet_metric{:,1};
     sheet_metric = renamevars(sheet_metric,"Time (wks)","Time (days)");
end  

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% opt_parameter = 32;
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

simulation_settings.tree_diamrat_prev_flag = 1; %Use diameter ratio from a previous solve 1
simulation_settings.tree_solve_prev_flag = 1; %Use down the tree material parameters from a previous solve 1
simulation_settings.tree_gen_prev_flag = 1; %Use a morphometric tree from a previous solve 1
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
%--------------------------------------------------------------------------

ts = 0: simulation_settings.step_size: simulation_settings.max_days / simulation_settings.step_size - 1;
simulation_settings.dQ_s = simulation_settings.dQ * (1 - exp( -simulation_settings.k_ramp * ts));
%simulation_settings.dQ_s = dQ_inputs;

% Store simulated data
sim_data = cell(1, 10);
total_days = days_run + 1;

selected_metric = table2cell(sheet_metric(:, 2));

metric_sim_vals = cell(1, 5); 

opt_function = @(opt_parameter) objective_combined(sim_data, selected_metric, params_uniform_infl, opt_parameter, simulation_settings);

opt_parameter0 = 1; 

% Run optimization
opt_gamma = lsqnonlin(opt_function, opt_parameter0);

J = objective_combined(sim_data, selected_metric, params_uniform_infl, opt_parameter0, simulation_settings);  

% put simulation_settings into objective function 
function J = objective_combined(sim_data, selected_metric, params, opt_parameter, simulation_settings)

    params.gamma_i = opt_parameter(1);
    
    % Run the gnr code with opt_parameters as the input
    [~, ~] = fun_run_tree_GnR(simulation_settings, params);

    pressure_fold = zeros(1, total_days); 
    
    for order = 1:10
        file = strcat('GnR_out_ord', num2str(order));
        sim_data{order} = load(file);
        
        for day = 1:total_days
            
            order_info = sim_data{order}(day, :);
            
            r_inner = order_info(1) * 10^6;
            h = order_info(2) * 10^6;
            r_outer = r_inner + h;
            od = 2 * r_outer;
            
            % new_col = numel(sim_data{order}(day, :)) + 1;
            % sim_data{order}(day, new_col) = order_info(21) + order_info(22) + order_info(23) + order_info(24);
            % sim_data{order}(day, new_col + 1) = percent_thick;

            if order == 10
                if day == 1
                    pressure_fold(1) = order_info(9);
                end
                pressure_fold(day) = order_info(9) / pressure_fold(1);
            end
        end
    end

    % Compute metrics
    [metric_sim_vals(1), metric_sim_vals(2)] = sim_results(total_days, sim_data, 26, od);
    [metric_sim_vals(3), metric_sim_vals(4)] = sim_results(total_days, sim_data, 15, od);
    metric_sim_vals(5) = pressure_fold;

    % Compute additional metrics as needed
    [~, metric_sim_vals(6)] = sim_results(total_days, sim_data, 16, od);
    [~, metric_sim_vals(7)] = sim_results(total_days, sim_data, 18, od);
    [~, metric_sim_vals(8)] = sim_results(total_days, sim_data, 25, od);

    % Calculate objective function values
    J = zeros(1, 5);
    
    for i = 1:length(selected_metric)
        J_curr = 0;
        if selected_metrics(i) > 0
            J_curr = (metric_sim_vals(i) - selected_metric(i)) / selected_metric(i);
        end
        
        % Store current objective value
        J(i) = J_curr;
    end
end
