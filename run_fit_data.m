clear all
clc

path = '/Users/brettcardenas/Desktop/Stanford_CVI/Extracted_MCT_Data.xlsx';
% path = '/Users/brettcardenas/Desktop/Stanford_CVI/Extracted_SuHx_Data.xlsx';
sheets = sheetnames(path);
days_run = 40;
inputs = {'Cardiac output (ml/min)','Cardiac index'};
wt_vars = {'WT (%), <100 um OD','WT (%), >100 um OD'};
stiff_vars = {'Shear mod (kPa), <100 um OD','Shear mod (kPa), >100 um OD'};
pressure_vars = {'mPAP (mmHg)', 'RVSP (mmHg)'};
ecs_vars = {'Elastase activity ratio to control','Colagenase activity ratio to control','COL1A1 mRNA'};
metrics = [wt_vars,stiff_vars,'mPAP (mmHg)',ecs_vars];
num_mets = numel(metrics);
excel_data = cell(1,num_mets);
metrics = [metrics inputs];
co_data = [];

%Loop for each metric
for metric = 1:numel(metrics)
    curr_metric = metrics{metric};
    total_metric = [];

    %Loop over each sheet
    for sheet = 1:numel(sheets)
        curr_sheet = sheets(sheet);
        time = 'Time (wks)';
    
        %Keep headers intact
        opts = detectImportOptions( ...
            path, ...
            "VariableNamingRule","preserve", ...
            "Sheet",curr_sheet ...
            );

        column_headers = opts.VariableNames;

        %Check if metric is in current sheet, else skip
        if ~any(strcmp(column_headers,curr_metric))
            continue
        end

        if ~any(strcmp(column_headers, time))
            time = 'Time (days)';
        end

        if strcmp(curr_metric,'Cardiac index')
            opts.SelectedVariableNames = {time, curr_metric, 'BW (g)'};
        else
            opts.SelectedVariableNames = {time, curr_metric};
        end

        sheet_metric = readtable(path,opts);

        %Remove -1 values (no data); if error, fix the metrics cell array
        sheet_metric(ismember(sheet_metric{:,curr_metric},-1),:) = [];

        %Normalize CI data by time-matched control (e.g. "Control Day 8")
        if strcmp(curr_metric,'Cardiac index')
            sheet_metric.CO = sheet_metric{:,2} .* sheet_metric{:,3};
            sheet_metric(:,[2,3]) = [];
            numrows = height(sheet_metric);
            time_pts = 0.5 * numrows;
            treatment_gp = sheet_metric(1:time_pts,'CO');
            ctrl_gp = sheet_metric((time_pts+1):numrows,'CO');
            norm_co = treatment_gp ./ ctrl_gp;
            sheet_metric(1:time_pts,"CO") = norm_co;
        end

        %Normalize CO data to only one "Day 0" control
        if strcmp(curr_metric,'Cardiac output (ml/min)')
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

        %Add table to growing matrix for specific metric
        if any(strcmp(curr_metric,inputs))
            co_data = [co_data; sheet_metric];
        else
            total_metric = [total_metric; sheet_metric];
        end

    end

    %Remove 35 wk data point (only applicable for SuHx data)
    %total_metric(ismember(total_metric{:,"Time (days)"},245),:) = [];

    if ~any(strcmp(curr_metric,inputs))
        [temp,~,idx] = unique(total_metric{:,1});
        excel_data{metric} = [temp, accumarray(idx,total_metric{:,2},[],@mean)];
    end

end

[temp,~,idx] = unique(co_data{:,1});
input_data = [temp, accumarray(idx,co_data{:,2},[],@mean)];

Q_exp_time = input_data(:,1);
dQ_exp = input_data(:,2) - input_data(1,2);

dQ_interpolated_times = 0:5:84;
dQ_interpolated = interp1(Q_exp_time,dQ_exp,dQ_interpolated_times,'pchip');
dQ_smooth = smooth(dQ_interpolated, 5);
dQ_inputs_times = 0:1:84;
dQ_inputs = interp1(dQ_interpolated_times,dQ_smooth,dQ_inputs_times,'pchip');

colors = sky(2) .* 0.8;
figure(16)
hold on

plot(Q_exp_time,dQ_exp,'ko', ...
    'MarkerSize', 8, ...
    'LineWidth', 2)
plot(dQ_interpolated_times,dQ_interpolated,'--', ...
    'Color', colors(1,:), ...
    'LineWidth', 2)
plot(dQ_inputs_times,dQ_inputs,'-', ...
    'Color', colors(2,:), ...
    'LineWidth', 2)

axis([0 87 -0.5 0.5])
xlabel('Time (days)'); ylabel('dQ');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', 0:7:84, ...
    'Ytick', -0.5:0.1:0.5);

hold off



%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%2. Run the coupled model script and read in the simulated outputs
%--------------------------------------------------------------------------

selected_metrics = [1 0 0 0 0 0 0 0]; %these are 1 if it's included and 0 if not
opt_function = @(opt_parameter) objective_combined(selected_metrics, metrics, params, opt_parameter);
%function J = objective_combined(selected_metrics, metrics, params, opt_parameter)
%J = 0
%params.Ki_steady = opt_parameter(1)

%run the gnr code with opt_parameters as the input
%read in the simulation results

%for i = 1:length(selected_metrics)
%J_curr = 0
%if selected_metrirics(i) > 0
%J_curr = ((fold change simulation - fold change exp) / fold change exp)^2
%end
%J = J + J_curr
%end
opt_parameter0 = [1, 1];
opt_gamma = lsqnonlin(opt_function, opt_parameter0);

J = J_combined(objective_combined(selected_metrics, metrics, opt_parameter));

params_uniform_infl.Ki_Tact = 0.0;
params_uniform_infl.phi_Tact0_min = 1.00;
params_uniform_infl.gamma_i = opt_parameter;
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
%ts = 0: simulation_settings.step_size: simulation_settings.max_days / simulation_settings.step_size - 1;
%simulation_settings.dQ_s = simulation_settings.dQ * (1 - exp( -simulation_settings.k_ramp * ts));
simulation_settings.dQ_s = dQ_inputs;


curr_param = [0.1, 0.2, 0.4, 0.8, 1.6]; % Ki_steady params
num_params = numel(curr_param);
sim_data = cell(num_params,10);

colors = sky(num_params + 1);

total_days = days_run + 1;


% 
for val = 1:num_params
    params_uniform_infl.Ki_steady = curr_param(val);
    simulation_settings.loop_iteration = val;
    
    [hemo_save, geom_save] = fun_run_tree_GnR(simulation_settings,params_uniform_infl);

    for order = 1:10
        file = strcat('GnR_out_ord',num2str(order));
        sim_data{val,order} = load(file);
        
        for day = 1:total_days
            order_cols = numel(sim_data{val, order}(day,:)) + 1;
            sim_data{val,order}(day,order_cols) = ;
        end
    end

    wt_over_mean = zeros([1 total_days]);
    wt_under_mean = zeros([1 total_days]);
    stiff_over_mean = zeros([1 total_days]);
    stiff_over_fold = zeros([1 total_days]);
    stiff_under_mean = zeros([1 total_days]);
    stiff_under_fold = zeros([1 total_days]);
    pressure_fold = zeros([1 total_days]);
    elastin_degr_mean = zeros([1 total_days]);
    elastin_degr_fold = zeros([1 total_days]);
    collagen_degr_mean = zeros([1 total_days]);
    collagen_degr_fold = zeros([1 total_days]);
    collagen_prod_mean = zeros([1 total_days]);
    collagen_prod_fold = zeros([1 total_days]);
    
    for day = 1:total_days

        wt_over100 = [];
        wt_under100 = [];
        stiff_over100 = [];
        stiff_under100 = [];
        pressure_0 = sim_data{val,10}(1,9);
        elastin_degr_under100 = [];
        collagen_degr_under100 = [];
        collagen_prod_under100 = [];
    
        for order = 1:numel(sim_data(val,:))
    
            % Gets data for current order
            order_info = sim_data{val,order}(day,:);
    
            r_inner = order_info(1) * 10^6;
            h = order_info(2) * 10^6;
            r_outer = r_inner + h;
            od = 2 * r_outer;
    
            outer = (r_outer)^2;
            inner = (r_inner)^2;
            percent_thick = 100 * (outer - inner)/inner;

            if od > 100
                wt_over100 = [wt_over100 percent_thick];
                stiff_over100 = [stiff_over100 order_info(15)];
            else
                wt_under100 = [wt_under100 percent_thick];
                stiff_under100 = [stiff_under100 order_info(15)];
                elastin_degr_under100 = [elastin_degr_under100 order_info(16)];
                collagen_degr_under100 = [collagen_degr_under100 order_info(18)];
                
                collagen_prod_total = order_info(21) + order_info(22) + order_info(23) + order_info(24);
                collagen_prod_under100 = [collagen_prod_under100 collagen_prod_total];
            end

            if order == 10
                pressure = order_info(9);
                fold = pressure / pressure_0;
                pressure_fold(day) = fold;
            end
    
        end
    
        if any(wt_over100)
            wt_over_mean(day) = mean(wt_over100);
        end
        if any(wt_under100)
            wt_under_mean(day) = mean(wt_under100);
        end
        if any(stiff_over100)
            stiff_over_mean(day) = mean(stiff_over100);
            stiff_over_fold(day) = stiff_over_mean(day) / stiff_over_mean(1);
        end
        if any(stiff_under100)
            stiff_under_mean(day) = mean(stiff_under100);
            stiff_under_fold(day) = stiff_under_mean(day) / stiff_under_mean(1);
        end
        if any(elastin_degr_under100)
            elastin_degr_mean(day) = mean(elastin_degr_under100);
            elastin_degr_fold(day) = elastin_degr_mean(day) / elastin_degr_mean(1);
        end
        if any(collagen_degr_under100)
            collagen_degr_mean(day) = mean(collagen_degr_under100);
            collagen_degr_fold(day) = collagen_degr_mean(day) / collagen_degr_mean(1);
        end
        if any(collagen_prod_under100)
            collagen_prod_mean(day) = mean(collagen_prod_under100);
            collagen_prod_fold(day) = collagen_prod_mean(day) / collagen_prod_mean(1);
        end

    
    end
    
    wt_over_mean = wt_over_mean';
    wt_under_mean = wt_under_mean';
    stiff_over_fold = stiff_over_fold';
    stiff_under_fold = stiff_under_fold';
    pressure_fold = pressure_fold';
    elastin_degr_fold = elastin_degr_fold';
    collagen_degr_fold = collagen_degr_fold';
    collagen_prod_fold = collagen_prod_fold';

    x = 0:days_run;
    
    % y(row, column, set)
    y = NaN(length(x),numel(metrics));
    
    y(:,1) = wt_under_mean(:,1);
    y(:,2) = wt_over_mean(:,1);
    y(:,3) = stiff_under_fold(:,1);
    y(:,4) = stiff_over_fold(:,1);
    y(:,5) = pressure_fold;
    y(:,6) = elastin_degr_fold(:,1);
    y(:,7) = collagen_degr_fold(:,1);
    y(:,8) = collagen_prod_fold(:,1);

    % UPDATE HERE
    for metric_number = 1:numel(metrics)
        figure (7 + metric_number)

        if metric_number == numel(metrics)
            semilogy(x, y(:,metric_number), ...
            'LineWidth', 2, ...
            'Color', colors(val,:), ...
            'LineStyle','-')
            hold on
        else
            hold on
            plot(x, y(:,metric_number), ...
                'LineWidth', 2, ...
                'Color', colors(val,:), ...
                'LineStyle','-')
        end
    end

end

wt_under_exp_days = excel_data{1}(:,1);
wt_under_exp_data = excel_data{1}(:,2);
wt_over_exp_days = excel_data{2}(:,1);
wt_over_exp_data = excel_data{2}(:,2);
stiff_under_exp_days = excel_data{3}(:,1);
stiff_under_exp_data = excel_data{3}(:,2);
stiff_under_exp_0 = stiff_under_exp_data(1);
stiff_over_exp_days = excel_data{4}(:,1);
stiff_over_exp_data = excel_data{4}(:,2);
stiff_over_exp_0 = stiff_over_exp_data(1);
pressure_exp_days = excel_data{5}(:,1);
pressure_exp_data = excel_data{5}(:,2);
pressure_exp_0 = pressure_exp_data(1);
elastase_exp_days = excel_data{6}(:,1);
elastase_ratio_exp_data = excel_data{6}(:,2);
colagenase_exp_days = excel_data{7}(:,1);
colagenase_ratio_exp_data = excel_data{7}(:,2);
col1a1_mrna_exp_days = excel_data{8}(:,1);
col1a1_mrna_exp_data = excel_data{8}(:,2);
col1a1_mrna_exp_0 = col1a1_mrna_exp_data(1);

y_exp_wt_under = NaN([1 numel(x)]);
y_exp_wt_over = NaN([1 numel(x)]);
y_exp_stiff_under_fold = NaN([1 numel(x)]);
y_exp_stiff_over_fold = NaN([1 numel(x)]);
y_exp_pressure_fold = NaN([1 numel(x)]);
y_exp_elastase_ratio = NaN([1 numel(x)]);
y_exp_colagenase_ratio = NaN([1 numel(x)]);
y_exp_col1a1_fold = NaN([1 numel(x)]);


min_x = 0;
min_y = 0;

% Plot for WT under 100

for unit = 1:numel(wt_under_exp_days)
    curr_day = wt_under_exp_days(unit) + 1;
    curr_val = wt_under_exp_data(unit);
    y_exp_wt_under(curr_day) = curr_val;
end

figure(8)
hold on
plot(x, y_exp_wt_under, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(wt_under_exp_days)/10)*10;
max_y = 100;
x_step = 5;
y_step = 10;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Wall Thickness (%), OD < 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for WT over 100

for unit = 1:numel(wt_over_exp_days)
    curr_day = wt_over_exp_days(unit) + 1;
    curr_val = wt_over_exp_data(unit);
    y_exp_wt_over(curr_day) = curr_val;
end

figure(9)
hold on
plot(x, y_exp_wt_over, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(wt_over_exp_days)/10)*10;
max_y = 175;
x_step = 5;
y_step = 25;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Wall Thickness (%), OD > 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for stiffness under 100 fold change

for unit = 1:numel(stiff_under_exp_days)
    curr_day = stiff_under_exp_days(unit) + 1;
    curr_val = stiff_under_exp_data(unit);
    stiff_under_fold_exp = curr_val/stiff_under_exp_0;
    y_exp_stiff_under_fold(curr_day) = stiff_under_fold_exp;
end

figure(10)
hold on
plot(x, y_exp_stiff_under_fold, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(stiff_under_exp_days)/10)*10;
max_y = 9;
x_step = 5;
y_step = 1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Stiffness Fold Change, OD < 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for stiffness over 100 fold change

for unit = 1:numel(stiff_over_exp_days)
    curr_day = stiff_over_exp_days(unit) + 1;
    curr_val = stiff_over_exp_data(unit);
    stiff_over_fold_exp = curr_val/stiff_over_exp_0;
    y_exp_stiff_over_fold(curr_day) = stiff_over_fold_exp;
end

figure(11)
hold on
plot(x, y_exp_stiff_over_fold, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(stiff_over_exp_days)/10)*10;
max_y = 11;
x_step = 5;
y_step = 1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Stiffness Fold Change, OD > 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for pressure fold change

for unit = 1:numel(pressure_exp_days)
    curr_day = pressure_exp_days(unit) + 1;
    curr_val = pressure_exp_data(unit);
    pressure_fold_exp = curr_val/pressure_exp_0;
    y_exp_pressure_fold(curr_day) = pressure_fold_exp;
end

figure(12)
hold on
plot(x, y_exp_pressure_fold, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(pressure_exp_days)/10)*10;
max_y = 6;
x_step = 5;
y_step = 1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('mPAP Fold Change, Order 10');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for Elastase Ratio to Control

for unit = 1:numel(elastase_exp_days)
    curr_day = elastase_exp_days(unit) + 1;
    curr_val = elastase_ratio_exp_data(unit);
    y_exp_elastase_ratio(curr_day) = curr_val;
end

figure(13)
hold on
plot(x, y_exp_elastase_ratio, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(elastase_exp_days)/10)*10;
max_y = 1.5;
x_step = 5;
y_step = 0.1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Elastin Degradation Fold Change, OD < 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for Colagenase Ratio to Control

for unit = 1:numel(colagenase_exp_days)
    curr_day = colagenase_exp_days(unit) + 1;
    curr_val = colagenase_ratio_exp_data(unit);
    y_exp_colagenase_ratio(curr_day) = curr_val;
end

figure(14)
hold on
plot(x, y_exp_colagenase_ratio, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(colagenase_exp_days)/10)*10;
max_y = 1.4;
x_step = 5;
y_step = 0.1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Collagen Degradation Fold Change, OD < 100 um');
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', min_y:y_step:max_y);



% Plot for COL1A1 mRNA fold change

for unit = 1:numel(col1a1_mrna_exp_days)
    curr_day = col1a1_mrna_exp_days(unit) + 1;
    curr_val = col1a1_mrna_exp_data(unit);
    col1a1_exp = curr_val / col1a1_mrna_exp_0;
    y_exp_col1a1_fold(curr_day) = col1a1_exp;
end

figure(15)
hold on
semilogy(x, y_exp_col1a1_fold, ...
    'LineWidth', 2, ...
    'Color', [0 0 0], ...
    'LineStyle','none', ...
    'Marker', 'o', ...
    'MarkerSize',8)

max_x = ceil(max(col1a1_mrna_exp_days)/10)*10;
min_y = 10^-3;
max_y = 10;
x_step = 5;
%y_step = 1;

axis([min_x (max_x + x_step/2) min_y max_y])
xlabel('Time (days)'); ylabel('Collagen Production Fold Change, OD < 100 um');
format shortEng
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'FontWeight', 'bold', ...
    'Xtick', min_x:x_step:max_x, ...
    'Ytick', [10.^(-3:0), 5, 10]);



save_dir = "max" + simulation_settings.max_days + "save" + simulation_settings.save_steps + "dQ" + simulation_settings.dQ + "dP" + simulation_settings.dP + simulation_settings.other + "figures";

mkdir(save_dir)

save_fig8 = "wt_under.png";
saveas(figure(8), save_fig8)
copyfile(save_fig8, save_dir);

save_fig9 = "wt_over.png";
saveas(figure(9), save_fig9)
copyfile(save_fig9, save_dir);

save_fig10 = "stiff_under.png";
saveas(figure(10), save_fig10)
copyfile(save_fig10, save_dir);

save_fig11 = "stiff_over.png";
saveas(figure(11), save_fig11)
copyfile(save_fig11, save_dir);

save_fig12 = "pressure.png";
saveas(figure(12), save_fig12)
copyfile(save_fig12, save_dir);

save_fig13 = "elastin_degradation.png";
saveas(figure(13), save_fig13)
copyfile(save_fig13, save_dir);

save_fig14 = "collagen_degradation.png";
saveas(figure(14), save_fig14)
copyfile(save_fig14, save_dir);

save_fig15 = "collagen_production.png";
saveas(figure(15), save_fig15)
copyfile(save_fig15, save_dir);

J_pressure = (pressure_fold -  y_exp_pressure_fold) ./ y_exp_pressure_fold;
J_pressure = J_pressure(J_Pressure == J_Pressure);

J_WT = (wt_fold -  wt_over_exp_data) ./ wt_over_exp_data;
J_WT = J_WT(J_WT == J_WT);
%Pad with zeros to get same length of arrays for each metric

J_combined = [J_pressure; J_WT];

    % over100_no_nan = rmmissing(y(:,:,1));
    % under100_no_nan = rmmissing(y(:,:,2));
    % mpap_no_nan = rmmissing(y2);
    % 
    % %J = ((exp - sim) / sim)^2
    % over100_sim = over100_no_nan(:,1);
    % over100_exp = over100_no_nan(:,2);
    % under100_sim = under100_no_nan(:,1);
    % under100_exp = under100_no_nan(:,2);
    % mpap_sim = mpap_no_nan(:,1);
    % mpap_exp = mpap_no_nan(:,2);
    % 
    % J_over100 = ((over100_exp - over100_sim) ./ over100_sim).^2;
    % J_under100 = ((under100_exp - under100_sim) / under100_sim)^2;
    % J_mpap = ((mpap_exp - mpap_sim) / mpap_sim)^2;
    % 
    % J_over100 = sum(J_over100)/numel(J_over100);
    % J_under100 = sum(J_under100)/numel(J_under100);
    % J_mpap = sum(J_mpap)/numel(J_mpap);


%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%3. Compare the experimental data to the simulated outputs
%--------------------------------------------------------------------------

%Now you need to extract stuff from the simulation results that is
%comparable to the experimental data. The 1st column in GnR_out_ordn is the
%inner radius and the 2nd column is the thickness (in meters). You should 
%be able to calculate the same thickness metric (which is thickness % of 
%outer radius iirc). Each line of GnR_out_ordn is a time step, which you
%specify as part of the simulations settings, default to 1 day. So
%calculate that metric from the simulations results and then calculate:
%J = ((exp - sim) / sim)^2


% and plot the experimental and simulated results (Matlab has good built in
% plotting

%--------------------------------------------------------------------------

%Steps 2 and 3 will eventually need to be in their own function so that we
%can use an optimization method to iteratively input parameters and see if
%we can make the coupled model better match the data.