function presence = generatePatternsCompositionModified(params, density, N_patterns, tolerance, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
% (tolerance) - optional tolerance for density matching (default 0.005)
% (mesh) - optionally provided mesh to specify size of patterns
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

% Load in the seed data
load('fibro_seedinfo.mat', 'permute_tables', 'offset_tables');

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Set the tolerance if it wasn't provided
if nargin < 4
    tolerance = 0.005;
end

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 5
    mesh = buildMesh(250, 400, 1/136);
end

% Create the requested number of patterns
for m = 1:N_patterns

    % Initialise the patterns storage
    patterns{m} = {};
    seeds{m} = {};
    
    % Use the fibre-free generator if NaNs are present in input params
    % vector, or if only non-fibre parameters provided, otherwise 
    % use the standard generator
    if density > 0.1
        threshold = 0.1;
    elseif density > 0.05
        threshold = 0.05;
    else
        threshold = 0.01;
    end
    seed = m;
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibresModified(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    end

    % Get the density of the pattern generated
    actual_density = getPatternDensity(presence);

    % While the density is not the desired one, keep generating patterns and summing them
    max_iterations = 10;
    iteration = 0;
    while abs(actual_density - density) > tolerance

        % Use the fibre-free generator if NaNs are present in input params
        % vector, or if only non-fibre parameters provided, otherwise 
        % use the standard generator
        abs_diff = abs(actual_density - density);
        if abs_diff > 0.1
            threshold = 0.1;
        elseif abs_diff > 0.05
            threshold = 0.05;
        else
            threshold = 0.01;
        end
        seed = randi([1, 250]);
        if any(isnan(params))
            [aux_presence, ~, ~] = createFibroPatternNoFibresModified(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
        elseif  length(params) == 6
            [aux_presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
        else
            [aux_presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
        end

        % Sum the patterns
        aux_presence = presence + aux_presence;
        aux_presence(aux_presence == 2) = 1; % if element is 2, it means that both patterns have 1 in that position, so it is a 1

        % Get the density of the sum of the patterns
        aux_density = getPatternDensity(aux_presence);

        iteration = iteration + 1;

        % Only update the presence if the density obtained is less than the desired density 
        if aux_density < density+tolerance
            presence = aux_presence;
            actual_density = aux_density;
            iteration = 0;
        end

        % If the number of iterations is greater than the maximum number of iterations, restart the presence
        if iteration > max_iterations

            % Use the fibre-free generator if NaNs are present in input params
            if density > 0.1
                threshold = 0.1;
            elseif density > 0.05
                threshold = 0.05;
            else
                threshold = 0.01;
            end
            if any(isnan(params))
                [presence, ~, ~] = createFibroPatternNoFibresModified(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
            elseif  length(params) == 6
                [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
            else
                [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
            end
            actual_density = getPatternDensity(presence);
            iteration = 0;

        end
    end
end

end