function J = J_combined(metrics_on_certain_day)
metric = 1:height(metrics_on_certain_day);
sim = metrics_on_certain_day(:,1);
exp = metrics_on_certain_day(:,2);
J = sqrt((sim(metric) - exp(metric)) ./ exp(metric));
end