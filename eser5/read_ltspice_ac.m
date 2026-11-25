function [f, Z, names] = read_ltspice_ac(file)
% READ_LTSPICE_AC  Reads ALL the traces of a text file written by LTspice
% "File > Export data as Text" for an AC analysis.
%
%   [f, Z, names] = read_ltspice_ac(file)
%
%   f     : frequency vector [Hz] (first column of the file), Nf x 1
%   Z     : Nf x Ntraces matrix, one column per exported trace; complex when the export
%           contains complex pairs, real when it contains plain magnitudes
%   names : cell array with the trace names taken from the header, e.g. {'V(vdie)'}
%
%   The layout is detected automatically:
%     - (real, imaginary) pairs    3.9e-002,-2.5e-001
%     - (dB, phase) pairs          (-1.18e+001dB,-8.11e+001<degree sign>)
%     - plain real numbers         2.3e+000

txt = fileread(file);
if isempty(txt)
    error('read_ltspice_ac:empty', 'The file %s is empty.', file);
end
nl   = find(txt == char(10), 1);                 % end of the header line
hdr  = strtrim(txt(1:nl-1));
body = txt(nl+1:end);
names = strtrim(strsplit(hdr, sprintf('\t')));
names = names(~cellfun(@isempty, names));        % ignore empty header fields
nc    = numel(names) - 1;                        % number of traces after the frequency column

% layout detection on the first data line
nl2   = find(body == char(10), 1);
if isempty(nl2), nl2 = numel(body) + 1; end
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
    error('read_ltspice_ac:format', ...
        ['Unexpected number of values in %s (layout "%s", %d traces in the header: %s). ' ...
         'Check that every trace was exported as a number or as a complex pair.'], ...
        file, layout, nc, strjoin(names(2:end), ', '));
end
A = reshape(vals, 1 + per*nc, []).';             % one row per frequency
f = A(:,1);

Z = zeros(size(A,1), nc);
for c = 1:nc
    switch layout
        case 'real'
            Z(:,c) = A(:, 1 + c);
        case 'reim'
            Z(:,c) = A(:, 2*c) + 1i*A(:, 2*c + 1);
        case 'dB'
            Z(:,c) = 10.^(A(:, 2*c)/20) .* exp(1i * A(:, 2*c + 1) * pi/180);
    end
end
names = names(2:end);                            % trace names only (without the frequency column)
end
