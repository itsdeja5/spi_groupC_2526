%% Stand-alone impedance of the decoupling capacitors from vendor data
% Source: IdEM workspace capacitors_workspace.idem.mat (vendor S-parameters)
% Run from the SW2_STM root folder (the one containing the "capacitors" folder).
clear; close all; clc

%% 1) Load the workspace
% The file is a MATLAB v7.3 file containing one struct called "idemws".
% idemws.DUM is a 6x1 cell: one entry per capacitor, each a plain struct
% with fields: Name, f (1x401, Hz), S (2x2 cell of 1x401 complex), R0 (50 ohm).
% The IdEM toolbox is NOT needed just to read the data.
wsfile = fullfile('../input_files/capacitors', 'capacitors_workspace.idem.mat');
load(wsfile, 'idemws');
N = numel(idemws.DUM);

%% 2) Convert S-parameters to impedance
% The vendor measured the parts in SERIES MODE: the capacitor is placed
% between port 1 and port 2. For a series impedance Z between two ports of
% reference impedance R0:   S21 = 2*R0/(Z + 2*R0)   =>   Z = 2*R0*(1-S21)/S21
f    = cell(N,1);
Z    = cell(N,1);
name = cell(N,1);
C_nF = zeros(N,1);
for k = 1:N
    D       = idemws.DUM{k};
    name{k} = D.Name;
    f{k}    = D.f(:);
    S21     = D.S{2,1}(:);
    Z{k}    = 2*D.R0*(1 - S21) ./ S21;
    % capacitance code in the Murata part number (characters 11-13):
    % two significant digits followed by the number of zeros, in pF
    C_nF(k) = str2double(D.Name(11:12)) * 10^str2double(D.Name(13)) * 1e-3;
end

% order from the largest to the smallest capacitance (as in the slide)
[~, idx] = sort(C_nF, 'descend');

%% 3) Plot |Z| vs frequency (same axes as the slide: MHz, ohm)
figure('Color', 'w'); hold on
colors = lines(N);
for m = 1:N
    k = idx(m);
    plot(f{k}/1e6, abs(Z{k}), 'LineWidth', 1.5, 'Color', colors(m,:), ...
        'DisplayName', sprintf('%s (%.4g nF)', name{k}, C_nF(k)));
end
% log-log axes (set explicitly: with 'hold on', loglog() would NOT change the scale)
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on; box on
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
xlim([0.1 1e4]); ylim([1e-2 1e4]);
xlabel('Frequency [MHz]'); ylabel('|Z| [\Omega]');
title('Stand-alone capacitor impedance (vendor S-parameters)');
legend('Location', 'southwest', 'Interpreter', 'none');

%% 4) Self-resonance frequency and minimum impedance
% Note: the vendor grid has only 401 log-spaced points (about 4.6 % per step),
% so the minimum found on the grid is accurate only to that resolution.
fprintf('%-20s %8s %12s %12s\n', 'Part', 'C [nF]', 'fSRF [MHz]', '|Z|min [ohm]');
for m = 1:N
    k = idx(m);
    [zmin, i] = min(abs(Z{k}));
    fprintf('%-20s %8.3g %12.2f %12.4f\n', name{k}, C_nF(k), f{k}(i)/1e6, zmin);
end

%% 5) Save the results for the later comparison with the SPICE models
save('cap_Z_vendor_data.mat', 'f', 'Z', 'name', 'C_nF');
