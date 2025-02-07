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

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 5
    mesh = buildMesh(250, 400, 1/136);
end

% Create the requested number of patterns
for m = 1:N_patterns

    % Initialise counters
    initial_seed = m;
    restart_counter = 0;
    tries_counter = 0;
    max_tries = 10;
    iterations = 0;
    
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
    seed = getFeasibleSeed(initial_seed);
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    end

    % Get the density of the pattern generated
    actual_density = getPatternDensity(presence);

    % While the density is not the desired one, keep generating patterns and summing them
    while abs(actual_density - density) > tolerance

        iterations = iterations + 1;

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
        seed = getFeasibleSeed(initial_seed + iterations);
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
            tries_counter = 0;
        end
    end
end

end

function feasible_seed = getFeasibleSeed(seed)
    % This function returns a feasible seed for the random number generator
    % based on the provided seed. If the seed is greater than 250, it returns
    % the remainder of the division by 250. If the seed is less than 1, it
    % returns 1.
    %
    % INPUTS:
    %
    % seed - the seed to be checked
    %
    % OUTPUTS:
    %
    % feasible_seed - the feasible seed to be used
    
    if seed > 250
        feasible_seed = mod(seed, 250);
        if feasible_seed == 0
            feasible_seed = 1;
        end
    elseif seed < 1
        feasible_seed = 1;
    else
        feasible_seed = seed;
    end
    
    end