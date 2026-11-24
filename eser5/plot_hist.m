% plot_hist.m  -  histogram of the peak-to-peak ripple of all the decap configurations
% Reads all_vpp.txt (written by best_decaps.py): first line is a header, then one row per
% configuration: "Vpp_V C01 C02 C03 C04 C05"

fid = fopen('all_vpp.txt', 'r');
fgetl(fid);                                           % skip the header line
D = textscan(fid, '%f %s %s %s %s %s');
fclose(fid);

vpp = D{1} * 1e3;                                     % ripple [mV]
cfg = [D{2} D{3} D{4} D{5} D{6}];                     % N x 5 cell array of capacitor labels

[vbest, ibest] = min(vpp);                            % best configuration
iuni10 = find(all(strcmp(cfg, '10nF'), 2), 1);        % five 10 nF capacitors

fprintf('%d configurations: min %.1f mV, median %.1f mV, max %.1f mV\n', ...
        numel(vpp), min(vpp), median(vpp), max(vpp));
fprintf('best: %s  %.2f mV\n', strjoin(cfg(ibest, :), ' '), vbest);

figure('Color', 'w', 'Position', [100 100 700 360]);
histogram(vpp, 80, 'FaceColor', [0.27 0.47 0.67], 'DisplayName', 'configurations');
hold on; grid on;
xline(vbest,        '--', 'Color', [0.80 0.20 0.07], 'LineWidth', 1.5, ...
      'DisplayName', sprintf('best: %.1f mV', vbest));
xline(vpp(iuni10),  '--', 'Color', [0.13 0.53 0.20], 'LineWidth', 1.5, ...
      'DisplayName', sprintf('5 x 10 nF: %.1f mV', vpp(iuni10)));
xlabel('Peak-to-peak ripple at the die [mV]');
ylabel('Number of configurations');
legend('Location', 'northeast');

exportgraphics(gcf, 'ripple_hist.png', 'Resolution', 200);