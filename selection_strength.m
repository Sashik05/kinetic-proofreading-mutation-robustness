clear all; clc;
s = linspace(0, 1, 100);
a = [0.001, 1, 5, 10, 1000];

% Create figure with proper settings for publication
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Professional color palette (RGB triplets)
colors = [0, 0.4470, 0.7410;    % Blue
          0.8500, 0.3250, 0.0980; % Orange-red
          0.9290, 0.6940, 0.1250; % Yellow
          0.4940, 0.1840, 0.5560; % Purple
          0.4660, 0.6740, 0.1880]; % Green

% Alternative: MATLAB's parula colormap for better distinction
% colors = parula(length(a));

hold on;

% Plot each curve
for i = 1:length(a)
    p_off = exp(a(i) .* s) ./ (1 + exp(a(i) .* s));
    plot(s, p_off, 'Color', colors(i,:), 'LineWidth', 3, ...
         'DisplayName', sprintf('a = %g', a(i)));
end

% Enhanced formatting for publication
xlabel('s', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('p_{off}', 'FontSize', 14, 'FontWeight', 'bold');

% Mathematical expression in title using LaTeX interpreter
% title('p_{off} = e^{as}/(1 + e^{as}) for Different a Values', ...
%       'FontSize', 16, 'FontWeight', 'bold');

% Enhanced legend
legend('Location', 'southeast', 'FontSize', 12, 'Box', 'on');

% Grid with subtle appearance
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--');

% Axis limits and ticks
xlim([0, 1]);
ylim([0.4, 1.02]); % Adjusted to avoid cutting off data

% Enhanced axis properties
set(gca, 'FontSize', 12, 'LineWidth', 1.5, ...
         'TickDir', 'out', 'Box', 'on');

% Set aspect ratio for better visualization
pbaspect([1.5 1 1]);

hold off;

% Optional: Save as high-quality PNG/EPS for publication
% print('sigmoid_plot', '-dpng', '-r300');
% print('sigmoid_plot', '-depsc', '-r300');