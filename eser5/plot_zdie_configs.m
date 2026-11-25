%% Impedance at the die port for different decoupling configurations
% Reads the LTspice text exports of the AC analysis (1 A injected at the die, so that
% V(vdie) = Z_die in ohm), exported in (real, imaginary) form, and plots |Z_die| of all
% the configurations in the same figure. read_ltspice_ac.m must be in the current folder
% (or on the MATLAB path).
clear; close all; clc

%% 1) Configurations to plot
% One row per curve: {file name, legend label, color [R G B], line width}
% To exclude a curve, comment its row with %  (or delete it).
% >>> EDIT THE FILE NAMES (first column) <<<
cfg = {
    'zdie_to_plot/zdie_nocaps.txt',     'no decaps',          [0.85 0.10 0.10], 2.2;
    'zdie_to_plot/zdie_best.txt',       'best configuration', [0.00 0.00 0.00], 2.4;
    'zdie_to_plot/zdie_all_1p5nF.txt',  '5 x 1.5 nF',         [0.00 0.45 0.74], 1.2;
    'zdie_to_plot/zdie_all_2p2nF.txt',  '5 x 2.2 nF',         [0.93 0.69 0.13], 1.2;
    'zdie_to_plot/zdie_all_3p3nF.txt',  '5 x 3.3 nF',         [0.47 0.67 0.19], 1.2;
    'zdie_to_plot/zdie_all_4p7nF.txt',  '5 x 4.7 nF',         [0.49 0.18 0.56], 1.2;
    'zdie_to_plot/zdie_all_6p8nF.txt',  '5 x 6.8 nF',         [0.30 0.75 0.93], 1.2;
    'zdie_to_plot/zdie_all_10nF.txt',   '5 x 10 nF',          [0.64 0.08 0.18], 1.2;
};
trace = 'V(vdie)';                % name of the exported trace (column header of the files)
f_mark = 100e6;                   % frequency at which |Z| is also reported [Hz]
xl = [0.1 5000];                  % x limits [MHz]   ([] = automatic)
yl = [1e-3 10];                   % y limits [ohm]   ([] = automatic)

%% 2) Read, print the summary and plot
figure('Color', 'w');
ax = axes;  hold(ax, 'on');
n = size(cfg, 1);
fprintf('%-22s %16s %12s %20s\n', 'Configuration', 'peak |Z| [ohm]', 'at f [MHz]', ...
    sprintf('|Z| at %g MHz [ohm]', f_mark/1e6));
for m = 1:n
% for m = [1, 2, 8]
    [f, Zall, names] = read_ltspice_ac(cfg{m,1});
    col = find(strcmpi(names, trace), 1);
    if isempty(col)
        if size(Zall, 2) == 1
            col = 1;                                  % only one trace in the file: use it
        else
            error('Trace %s not found in %s. Traces found: %s', trace, cfg{m,1}, strjoin(names, ', '));
        end
    end
    z = Zall(:, col);                                 % complex impedance (1 A injected) [ohm]
    plot(ax, f/1e6, abs(z), 'Color', cfg{m,3}, 'LineWidth', cfg{m,4}, 'DisplayName', cfg{m,2});

    [zpk, ipk] = max(abs(z));
    z_mark = exp(interp1(log(f), log(abs(z)), log(f_mark)));      % |Z| at f_mark (log-log interpolation)
    fprintf('%-22s %16.4f %12.2f %20.4f\n', cfg{m,2}, zpk, f(ipk)/1e6, z_mark);
end

% log scales set explicitly: with hold on, plot() keeps the scale you set
set(ax, 'XScale', 'log', 'YScale', 'log');
grid(ax, 'on');  box(ax, 'on');
if ~isempty(xl), xlim(ax, xl); end
if ~isempty(yl), ylim(ax, yl); end
xlabel(ax, 'Frequency [MHz]');
ylabel(ax, '|Z_{die}| [\Omega]');
title(ax, 'Die-port impedance for different decoupling configurations');
legend(ax, 'Location', 'northeast');

% optional: save the figure (MATLAB R2020a or later)
% exportgraphics(gcf, 'zdie_configs.png', 'Resolution', 300);
