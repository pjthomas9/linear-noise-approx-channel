%NOTE: Comment/Uncomment lines 20-22 to specify the channel

dp = 3e-2; % specify the increments of p (probability of high input) to sweep through
p_min = dp/2; % minumum value of p to plot
p_max = 1; % maximum value of p to plot
p_list = p_min:dp:p_max; % creates a vector of the values p over the specified range

% Creates a vector for the mutual information (MI) for fully and partially
% observed systems
SE_FOS_vec = zeros(1, length(p_list));
SE_FOS_exact_vec = zeros(1, length(p_list));

SE_POS_low_vec = zeros(1, length(p_list));
SE_POS_high_vec = zeros(1, length(p_list));

gap_ratio_low = zeros(1, length(p_list));
gap_ratio_high = zeros(1, length(p_list));

for j = 1:length(p_list)
    [SE_FOS, SE_POS_low, SE_POS_high, SE_FOS_exact] = SE_3state_chain(p_list(j)); %Computes MI for Fully and Partially Observed Chain
    %[SE_FOS, SE_POS_low, SE_POS_high] = SE_3state_ring(p_list(j)); %Computes MI for Fully and Partially Observed Ring
    %[SE_FOS, SE_POS_low, SE_POS_high] = SE_ACh(p_list(j)); %Computes MI for Fully and Partially Observed ACh model
    
  
    
    SE_FOS_vec(j) = real(SE_FOS);
    SE_POS_low_vec(j) = SE_POS_low;
    SE_POS_high_vec(j) = SE_POS_high;
    SE_FOS_exact_vec(j) = SE_FOS_exact;
    
    gap_ratio_low(j) = real(SE_FOS)/SE_POS_low;
    gap_ratio_high(j) = real(SE_FOS)/SE_POS_high;
end

%% Plots

figure
plot(p_list, SE_FOS_vec, 'k', 'LineWidth', 3)
hold on
plot(p_list, SE_POS_low_vec, '--', 'LineWidth', 5)
plot(p_list, SE_POS_high_vec, '--', 'LineWidth', 5)
%title('Three-State Chain: SE vs p; 2 \rightarrow 3 sensitive','FontSize',14)
xlabel('p','FontSize',18)
%ylabel('SE','FontSize',18)
legend({'Fully Observed', 'Partially Observed (low \omega)', 'Partially Observed (high \omega)'},'FontSize',12)
hold off
%shg
% figure
% plot(p_list, gap_ratio_low, 'LineWidth', 2)
% hold on
% plot(p_list, gap_ratio_high, 'LineWidth', 2)
% title('Three-State: Ratio of Fully to Partially Observed vs p','FontSize',14)
% xlabel('p','FontSize',14)
% ylabel('Ratio','FontSize',14)
% legend({'Low Frequency Limit', 'High Frequency Limit'},'FontSize',12)
% hold off



% filename = 'mi_vs_p_c001_0000';
% savefig(filename)
% saveas(gcf,filename,'jpg')
