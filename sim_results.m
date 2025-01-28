function [under_fold, over_fold] = sim_results(total_days, sim_data, col_num, od)

over_fold = zeros(1,total_days);
under_fold = zeros(1,total_days);

for day = 1:total_days

    over = NaN(1,numel(sim_data));
    under = NaN(1,numel(sim_data));

    for order = 1:numel(sim_data)
        
        if od >= 100
            over(order) = sim_data{order}(day,col_num);
        else
            under(order) = sim_data{order}(day,col_num);
        end

    end

    over_fold = mean(over) / over_fold(1);
    under_fold = mean(under) / under_fold(1);

end

end