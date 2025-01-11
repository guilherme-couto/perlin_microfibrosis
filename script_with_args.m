% script_with_args.m
if nargin < 5
    error('Not enough args. (params, param_name, density, N_patterns, fibro_typename)');
end

% Read args
params = str2num(argv(){2});  % Convert string to array
param_name = argv(){3};
density = str2double(argv(){4});
N_patterns = argv(){5};
fibro_typename = argv(){6};

% disp(argv());

% Call main function
generatePatternsModified(params, param_name, density, N_patterns, fibro_typename);
