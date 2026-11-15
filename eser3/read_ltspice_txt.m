function [f, S] = read_ltspice_txt(file)
% READ_LTSPICE_TXT  Reads a text file written by LTspice "File > Export data as Text" (AC analysis).
%
%   [f, S] = read_ltspice_txt(file)
%
%   f : frequency vector [Hz] (first column of the file)
%   S : struct with one field for every trace whose name contains zXXX (XXX = 3 digits),
%       for example S.z103 for the trace V(z103). The field is complex when the export
%       contains complex pairs, real when it contains plain magnitudes.
%
%   The layout is detected automatically:
%     - (real, imaginary) pairs    3.9e-002,-2.5e-001
%     - (dB, phase) pairs          (-1.18e+001dB,-8.11e+001<degree sign>)
%     - plain real numbers         2.3e+000         (e.g. mag(V(z103)) traces)
%   Columns that are not zXXX traces (for instance an extra "frequency" column) are ignored.

txt = fileread(file);
if isempty(txt)
    error('read_ltspice_txt:empty', 'The file %s is empty.', file);
end
nl   = find(txt == char(10), 1);                 % end of the header line
hdr  = strtrim(txt(1:nl-1));
body = txt(nl+1:end);
names = strtrim(strsplit(hdr, sprintf('\t')));
nc    = numel(names) - 1;                        % number of trace columns after the frequency

% layout detection on the first data line
nl2   = find(body == char(10), 1);
first = body(1:nl2-1);
if ~isempty(strfind(first, 'dB'))
    layout = 'dB';
elseif ~isempty(strfind(first, ','))
    layout = 'reim';
else
    layout = 'real';
end

% keep only numbers: every other character (parentheses, commas, 'dB', degree sign, ...) becomes a blank
body = regexprep(body, '[^0-9eE+\-.\s]', ' ');
if strcmp(layout, 'real'), per = 1; else, per = 2; end
vals = sscanf(body, '%f');
if mod(numel(vals), 1 + per*nc) ~= 0
    error('read_ltspice_txt:format', ...
        ['Unexpected number of values in %s (layout "%s", %d trace columns). ' ...
         'Check that every trace was exported as a number or as a complex pair.'], file, layout, nc);
end
A = reshape(vals, 1 + per*nc, []).';             % one row per frequency
f = A(:,1);

S = struct();
for c = 1:nc
    switch layout
        case 'real'
            x = A(:, 1 + c);
        case 'reim'
            x = A(:, 2*c) + 1i*A(:, 2*c + 1);
        case 'dB'
            x = 10.^(A(:, 2*c)/20) .* exp(1i * A(:, 2*c + 1) * pi/180);
    end
    tok = regexp(names{c+1}, 'z(\d{3})', 'tokens', 'once');
    if ~isempty(tok)
        S.(['z' tok{1}]) = x;
    end
end
fprintf('%s: %d frequencies (%.4g Hz to %.4g Hz), layout "%s", traces: %s\n', ...
    file, numel(f), f(1), f(end), layout, strjoin(fieldnames(S).', ' '));
end
