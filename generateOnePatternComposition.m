function pattern = generateOnePatternComposition(params, density, seed_num, tolerance, mesh)
% This function takes the provided set of parameter values, and creates the
% representation of a fibrosis pattern with the desired density.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the pattern to be generated
% seed_num - seed for the random number generator
% tolerance - tolerance for density matching
% (mesh) - optionally provided mesh to specify size of pattern
%
% PARAMETER INFORMATION:
%
% 1 - FIBRENESS: The extent to which patterns exhibit long, thin fibres 
%     aligned in consistent directions
%     ::: If set to NaN, a pattern without fibres will be created :::
% 2 - FIBRE SEPARATION: The average spacing between fibres (in units
%     matching those used in input mesh
% 3 - PATCHINESS: The extent of inconsistency in pattern density (high
%     patchiness will produce distinct regions of higher and lower density)
% 4 - FEATURE SIZE: The overall size of obstacle features in the mesh (in
%     units matching the mesh)
% 5 - ROUGHNESS: The roughness of feature edges (values from [0,1], may
%     cause issues for values of precisely 0 or 1)
% 6 - PATCH SIZE: The size of regions over which density varies (with 
%     extent according to PATCHINESS)
% 7 - FIBRE ALIGNMENT: The extent to which non-fibre features are aligned
%     to the fibre direction (i.e. extent of feature anisotropy)
% 8 - DIRECTION: An angle (in radians) defining the orientation of fibres
%     and/or feature anisotropy

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 5
    mesh = buildMesh(250, 400, 1/136);
end

% Initialise counters
restart_counter = 0;
tries_counter = 0;
max_tries = 10;
iterations = 0;

% Use the fibre-free generator if NaNs are present in input params
% vector, or if only non-fibre parameters provided, otherwise 
% use the standard generator
threshold = getThreshold(density);
[permute_table, offset_table] = generateTables(seed_num);
if any(isnan(params))
    [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
elseif length(params) == 6
    [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
else
    [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
end

% Get the density of the pattern generated
actual_density = getPatternDensity(presence);

% While the density is not the desired one, keep generating patterns and summing them
while abs(actual_density - density) > tolerance

    iterations = iterations + 1;

    % Generate a new pattern
    abs_diff = abs(actual_density - density);
    threshold = getThreshold(abs_diff);
    [permute_table, offset_table] = generateTables(seed_num + iterations);
    if any(isnan(params))
        [aux_presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
    elseif length(params) == 6
        [aux_presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
    else
        [aux_presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
    end

    % Sum the patterns
    aux_presence = presence + aux_presence;
    aux_presence(aux_presence == 2) = 1; % if element is 2, it means that both patterns have 1 in that position, so it is a 1

    % Get the density of the sum of the patterns
    aux_density = getPatternDensity(aux_presence);
    
    % Update the tries counter
    tries_counter = tries_counter + 1;

    % Only update the presence if the density obtained is less than the desired density 
    if aux_density < density+tolerance
        presence = aux_presence;
        actual_density = aux_density;
        tries_counter = 0;
    end

    % If the number of tries is greater than the maximum number of tries, restart the presence
    if tries_counter > max_tries

        restart_counter = restart_counter + 1;

        % Generate a new pattern
        threshold = getThreshold(density);
        if any(isnan(params))
            [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
        elseif length(params) == 6
            [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
        else
            [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
        end
        actual_density = getPatternDensity(presence);
        tries_counter = 0;
    end
end

% Create folders to save the images
output_folder = 'uncertainty_quantification';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fibrosis_pattern_folder = fullfile(output_folder, 'fibrosis_patterns');
if ~exist(fibrosis_pattern_folder, 'dir')
    mkdir(fibrosis_pattern_folder);
end

% Save the resulting patterns as images
fig = figure('Visible', 'off');
imagesc(presence);
axis('equal', 'off');
title(['Pattern ', num2str(seed_num), ' - Density: ', num2str(actual_density)]);
colormap(fibroclr);
saveas(gcf, fullfile(fibrosis_pattern_folder, ['fibrosis_pattern_', num2str(seed_num), '.png']));
close(fig);

fprintf('Pattern %d generated successfully! Saved in the %s folder.\n', seed_num, fibrosis_pattern_folder);

end

function [permute_table, offset_table] = generateTables(seed)
% This function generates the permutation and offset tables for the provided seed.
%
% INPUTS:
%
% seed - the seed to be used for the random number generator
%
% OUTPUTS:
%
% permute_table - the permutation table generated
% offset_table - the offset table generated

% Set the seed for the random number generator
rng(seed);

% Assume a decent safe number like eight for the number of offsets
N_freqs = 8;

% Permutation tables for this seed
for j = 1:N_freqs
    permute_table(j,:) = int32(randperm(256) - 1);
end

% Offset table for this seed
offset_table = rand(N_freqs, 2) - 0.5;
    
end

function threshold = getThreshold(density)
% This function returns the threshold to be used for the density of the pattern.
%
% INPUTS:
%
% density - the density of the pattern to be generated
%
% OUTPUTS:
%
% threshold - the threshold to be used for the density of the pattern

if density > 0.1
    threshold = 0.1;
elseif density > 0.05
    threshold = 0.05;
else
    threshold = 0.01;
end

end
