% script_uq.m
if nargin < 3
    error('Not enough args. (params, density, seed_num)');
end

% Read args
params = str2num(argv(){2});  % Convert string to number
density = str2double(argv(){3});
seed_num = round(str2double(argv(){4}));  % Ensure seed_num is an integer

% Call main function
generateOnePatternComposition(params, density, seed_num, 0.005);
