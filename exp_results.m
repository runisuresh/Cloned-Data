function fold_change = exp_results(exp_means)

days = exp_means(:, 1);
data = exp_means(:, 2);

min_days = min(days);
max_days = max(days);
fold_change = NaN(max_days,1);

for unit = 1:numel(days)
    curr_day = days(unit) - min_days + 1;
    curr_val = data(unit);
    fold_change(curr_day) = curr_val / data(1);

end


% 
% % Plot for Elastase Ratio to Control
% 
% for unit = 1:numel(elastase_exp_days)
%     curr_day = elastase_exp_days(unit) + 1;
%     curr_val = elastase_ratio_exp_data(unit);
%     y_exp_elastase_ratio(curr_day) = curr_val;
% end
% 
% 
% 
% 
% 
% % Plot for Colagenase Ratio to Control
% 
% for unit = 1:numel(colagenase_exp_days)
%     curr_day = colagenase_exp_days(unit) + 1;
%     curr_val = colagenase_ratio_exp_data(unit);
%     y_exp_colagenase_ratio(curr_day) = curr_val;
% end
% 
% 
% 
% 
% 
% % Plot for COL1A1 mRNA fold change
% 
% for unit = 1:numel(col1a1_mrna_exp_days)
%     curr_day = col1a1_mrna_exp_days(unit) + 1;
%     curr_val = col1a1_mrna_exp_data(unit);
%     col1a1_exp = curr_val / col1a1_mrna_exp_0;
%     y_exp_col1a1_fold(curr_day) = col1a1_exp;
% end

end