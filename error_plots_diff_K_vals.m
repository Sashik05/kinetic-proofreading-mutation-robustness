% --- Publication-ready Plot: Error Ratio vs m for multiple kx values ---
clear; clc;

% --- Parameters ---
alpha = 1e-4;         % Error discrimination parameter
kxp = 250;            % Phosphate hydrolysis rate constant
lx = 0.2;             % Dissociation rate for correct substrate
lxp = 0.01;             % Error pathway contribution
w = 0;                % No perturbation

% --- Range of m values ---
m = logspace(-5, 6, 5000);  % Driving rate constants

% --- kx values to compare ---
kx_values = [1, 15, 90];  % Different correct association rates
colors = lines(length(kx_values));  % Distinct colors
min_points = zeros(length(kx_values), 2);  % Store [m_opt, min_val]

% --- Create figure ---
figure('Color', 'w', 'Position', [100 100 760 540]); hold on;

for i = 1:length(kx_values)
    kx = kx_values(i);
    ky = kx / alpha;
    kyp = kxp;

    ly = lx / alpha;
    lyp = lxp;

    % Compute correct and error formation terms
    fx = ((m .* kxp) ./ ((kx + m) .* (lx + w))) + (lxp / (lx + w));
    fy = ((m .* kyp) ./ ((ky + m) .* (ly + w))) + (lyp / (ly + w));
    error_ratio = fy ./ fx;

    % Find minimum and annotate
    [min_val, min_idx] = min(error_ratio);
    m_opt = m(min_idx);
    min_points(i, :) = [m_opt, min_val];

    % Plot curve
    loglog(m, error_ratio, 'Color', colors(i,:), 'LineWidth', 2.5);
    
    % Mark minimum
    loglog(m_opt, min_val, 'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 8);
    
    % Annotate m_opt
    text(m_opt * 0.2, min_val * 0.5, ...
        sprintf('$m_0=%.2g s^{-1}$', m_opt), ...
        'Interpreter', 'latex', 'FontSize', 13, ...
        'Color', colors(i,:), 'FontWeight', 'bold');
end

% --- Axes and Labels ---
grid on;
set(gca, ...
    'FontSize', 14, ...
    'LineWidth', 1.5, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XScale', 'log', ...
    'YScale', 'log');

xlabel('Driving rate constant ($m$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Error Ratio ($f$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

legend(arrayfun(@(k) sprintf('$k_X = %d$', k), kx_values, ...
    'UniformOutput', false), 'Interpreter', 'latex', ...
    'FontSize', 13, 'Location', 'northeast');

% title('Error Ratio vs Driving Constant $m$', 'Interpreter', 'latex', 'FontSize', 16);

% Optional: save high-res vector output
% print(gcf, 'error_ratio_comparison', '-dpdf', '-r600');
