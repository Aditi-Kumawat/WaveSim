T = 1;
A = 1;
t_vect = 0:1e-3:T;
rect_vect = rectangularPulse(0, T, t_vect);
cos_vect = A * sin(pi/T * t_vect);
g_t_vect = cos_vect .* rect_vect;

figure
plot(t_vect, g_t_vect, 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% Define the ticks
xTicks = [0, T/2, T];
yLim = get(gca, 'YLim');
yTicks = [yLim(1), mean(yLim), yLim(2)];

% Set ticks
set(gca, 'XTick', xTicks, 'YTick', yTicks);
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% Set figure size
set(gcf,'Units','inches', 'Position', [18 3 3 2.5],...
    'PaperUnits', 'Inches', 'PaperSize', [3 2.5]);

cd SAVE_FIGS
saveas(gcf, 'HSP.svg');
cd ..
