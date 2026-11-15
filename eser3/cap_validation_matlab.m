%% Validation of the vendor SPICE capacitor models against the vendor data
% Vendor data : cap_Z_vendor_data.mat  (saved by cap_impedance_matlab.m: f, Z, name, C_nF)
% SPICE data  : LTspice text export of V(z103) ... V(z152) from the capacitor netlist
%               (caps_zgrid.cir or caps_spice_ac.cir); read by read_ltspice_txt.m
% Both files and read_ltspice_txt.m must be in the current folder (or edit the names below).
clear; close all; clc

matfile   = 'cap_Z_vendor_data.mat';
spicefile = 'caps_spice_ac.txt';                % <-- name of your exported LTspice file

%% 1) Vendor data, from the .mat file
V    = load(matfile);                        % fields: f, Z, name, C_nF
f    = V.f;      Zv = V.Z;      name = V.name;      C_nF = V.C_nF(:);
N    = numel(name);
code = cell(N,1);
for k = 1:N
    code{k} = name{k}(11:13);                % capacitance code in the part number, e.g. '103'
end
[~, idx] = sort(C_nF, 'descend');            % from the largest to the smallest capacitance

%% 2) SPICE data, from the LTspice text export
[fs, Zs] = read_ltspice_txt(spicefile);      % fs [Hz]; Zs.z103, Zs.z682, ... (complex or |Z|)
for k = 1:N
    if ~isfield(Zs, ['z' code{k}])
        error('Trace V(z%s) (%s) not found in %s. Traces found: %s', code{k}, name{k}, ...
            spicefile, strjoin(fieldnames(Zs).', ' '));
    end
end

%% 3) Comparison and plots
figure('Color', 'w');
ax1 = subplot(2,1,1);  hold(ax1, 'on');
ax2 = subplot(2,1,2);  hold(ax2, 'on');
colors = lines(N);
fprintf('\n%-20s %18s %14s %20s\n', 'Part', 'max |Z| error [%]', 'at f [MHz]', 'max phase err [deg]');
for m = 1:N
    k  = idx(m);
    zs = Zs.(['z' code{k}]);                 % SPICE data on the SPICE frequency grid fs
    % compare only where the vendor frequencies fall inside the SPICE frequency range
    ok = f{k}(:) >= fs(1) & f{k}(:) <= fs(end);
    fv = f{k}(ok);   zv = Zv{k}(ok);
    % SPICE |Z| at the vendor frequencies (interpolation in log-log)
    mag_i = exp(interp1(log(fs), log(abs(zs)), log(fv), 'linear'));
    err   = abs(mag_i - abs(zv)) ./ abs(zv) * 100;                   % |Z| error [%]
    [emax, i] = max(err);
    if isreal(zs)
        phmax = NaN;                         % magnitude-only export: no phase to compare
    else
        ph_i  = interp1(log(fs), unwrap(angle(zs)), log(fv), 'linear');
        phmax = max(abs(ph_i - unwrap(angle(zv)))) * 180/pi;         % phase error [deg]
    end
    fprintf('%-20s %18.3f %14.2f %20.3f\n', name{k}, emax, fv(i)/1e6, phmax);

    plot(ax1, f{k}/1e6, abs(Zv{k}), '-', 'Color', colors(m,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (%.4g nF)', name{k}, C_nF(k)));
    plot(ax1, fs/1e6, abs(zs), 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
    plot(ax2, fv/1e6, err, 'Color', colors(m,:), 'LineWidth', 1.2);
end
plot(ax1, NaN, NaN, 'k--', 'DisplayName', 'SPICE models (dashed)');

% log scales set explicitly: with hold on, plot() keeps the scale you set
set(ax1, 'XScale', 'log', 'YScale', 'log');
set(ax2, 'XScale', 'log');
grid(ax1, 'on');  box(ax1, 'on');  grid(ax2, 'on');  box(ax2, 'on');
xlim(ax1, [0.1 1e4]);  ylim(ax1, [1e-2 1e4]);
xlim(ax2, [0.1 1e4]);
ylabel(ax1, '|Z| [\Omega]');
title(ax1, 'Vendor S-parameter data (solid) vs vendor SPICE models (dashed)');
legend(ax1, 'Location', 'southwest');
xlabel(ax2, 'Frequency [MHz]');  ylabel(ax2, '|Z| error [%]');
