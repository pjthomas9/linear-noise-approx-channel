% Script to produce MI vs p plot for ACh fully observed AND
% partially observed

dp = 1e-2; % specify the increments of p (probability of high input) to sweep through
p_min = dp/2; % minumum value of p to plot
p_max = 1; % maximum value of p to plot
p_list = p_min:dp:p_max; % creates a vector of the values p over the specified range

MI_full_vec = zeros(1, length(p_list));
MI_part_vec = zeros(1, length(p_list));


for j = 1:length(p_list)
    [MI_full, MI_part] = MI_ACh_colored(p_list(j)); %Computes MI for Fully Observed Chain
    MI_full_vec(j) = MI_full;
    MI_part_vec(j) = MI_part;
end

%% Plots

figure
plot(p_list, MI_full_vec, 'LineWidth', 3)
hold on
plot(p_list, MI_part_vec, 'LineWidth', 3)
title('ACh: MI vs p','FontSize',14)
xlabel('p','FontSize',18)
ylabel('MI','FontSize',18)
legend({'Fully Observed', 'Partially Observed'},'FontSize',12)
%shg
