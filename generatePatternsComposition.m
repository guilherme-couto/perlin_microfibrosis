function patterns = generatePatternsComposition(params, density, N_patterns, tolerance, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
% tolerance - tolerance for density matching
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

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 5
    mesh = buildMesh(250, 400, 1/136);
end

% Variable to store each pattern used to generate the desired patterns
patterns = cell(N_patterns, 1);
seeds = cell(N_patterns, 1);

% Initialise a figure to plot some of the patterns
figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1], 'Visible', 'off');
num_plots_y = round(sqrt(N_patterns));
num_plots_x = ceil(N_patterns/num_plots_y);

% Create the requested number of patterns
for m = 1:N_patterns

    % Initialise counters
    initial_seed = m * 10; % x10 to avoid seed repetition in composition  
    restart_counter = 0;
    tries_counter = 0;
    max_tries = 10;
    iterations = 0;

    % Initialise the patterns storage
    patterns{m} = {};
    seeds{m} = {};
    
    % Use the fibre-free generator if NaNs are present in input params
    % vector, or if only non-fibre parameters provided, otherwise 
    % use the standard generator
    threshold = getThreshold(density);
    seed = initial_seed;
    [permute_table, offset_table] = generateTables(seed);
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
    end

    % Get the density of the pattern generated
    actual_density = getPatternDensity(presence);

    % Save the pattern's presence and seed
    patterns{m}{end+1} = presence;
    seeds{m}{end+1} = seed;

    % While the density is not the desired one, keep generating patterns and summing them
    while abs(actual_density - density) > tolerance

        iterations = iterations + 1;

        % Generate a new pattern
        abs_diff = abs(actual_density - density);
        threshold = getThreshold(abs_diff);
        seed = initial_seed + iterations;
        [permute_table, offset_table] = generateTables(seed);
        if any(isnan(params))
            [aux_presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
        elseif  length(params) == 6
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

            % Save the pattern's presence and seed
            patterns{m}{end+1} = presence;
            seeds{m}{end+1} = seed;
        end

        % If the number of tries is greater than the maximum number of tries, restart the presence
        if tries_counter > max_tries

            % Restart the pattern storage
            restart_counter = restart_counter + 1;
            patterns{m} = {};
            seeds{m} = {};

            % Generate a new pattern
            threshold = getThreshold(density);
            if any(isnan(params))
                [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
            elseif  length(params) == 6
                [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
            else
                [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
            end
            actual_density = getPatternDensity(presence);
            tries_counter = 0;

            % Save the pattern's presence and seed
            patterns{m}{end+1} = presence;
            seeds{m}{end+1} = seed;
        end

    end
    
    % Plot this pattern
    subplot(num_plots_y, num_plots_x, m);
    imagesc(presence);
    axis('equal', 'off');
    title(['D: ', num2str(actual_density)]);
    colormap(fibroclr);
end

% Create folders to save the images
output_folder = 'generated_patterns_composition';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

samples_folder = fullfile(output_folder, 'composition_samples');
if ~exist(samples_folder, 'dir')
    mkdir(samples_folder);
end

gifs_folder = fullfile(output_folder, 'composition_gifs');
if ~exist(gifs_folder, 'dir')
    mkdir(gifs_folder);
end

seeds_folder = fullfile(output_folder, 'composition_seeds');
if ~exist(seeds_folder, 'dir')
    mkdir(seeds_folder);
end

fibrosis_pattern_folder = fullfile(output_folder, 'fibrosis_patterns');
if ~exist(fibrosis_pattern_folder, 'dir')
    mkdir(fibrosis_pattern_folder);
end

% Save the figure (by Guilherme)
saveas(gcf, fullfile(output_folder, ['patterns_', num2str(N_patterns), '.png']));

% Save the presences and the seeds used to generate the patterns
for m = 1:N_patterns
    % GIF filename
    gif_filename = fullfile(gifs_folder, ['pattern_', num2str(m), '.gif']);

    % Loop para salvar os frames do GIF
    for n = 1:length(patterns{m})
        % Create figure to save the pattern image
        fig = figure('Visible', 'off');
        imagesc(patterns{m}{n});
        axis('equal', 'off');
        title(['Pattern ', num2str(m), ' - Density: ', num2str(getPatternDensity(patterns{m}{n}))]);
        colormap(fibroclr);
        saveas(gcf, fullfile(samples_folder, ['pattern_', num2str(m), '_presence_', num2str(n), '.png']));
        
        % Read the image
        img = imread(fullfile(samples_folder, ['pattern_', num2str(m), '_presence_', num2str(n), '.png']));
        
        % Save to GIF
        if n == 1
            imwrite(img, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
            colormap(fibroclr);
        else
            imwrite(img, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
            colormap(fibroclr);
        end
        
        close(fig);
    end
    close(gif_filename);

    % Save the seeds used to generate the patterns
    % Open file to save seeds used to generate the patterns
    seed_file = fullfile(seeds_folder, ['pattern_', num2str(m), '_seeds.txt']);
    fid = fopen(seed_file, 'w');
    fprintf(fid, 'Pattern Number | Seed Used\n');
    for n = 1:length(seeds{m})
        if n == 1
            fprintf(fid, '%d | %d\n', m, seeds{m}{n});
        else
            fprintf(fid, '+ | %d\n', seeds{m}{n});
        end
        if n == length(seeds{m})
            fprintf(fid, '(Total: %d)\n', length(seeds{m}));
        end
    end
    fclose(fid);

    % Save the resulting patterns as images
    fig = figure('Visible', 'off');
    imagesc(patterns{m}{end});
    axis('equal', 'off');
    title(['Pattern ', num2str(m), ' - Density: ', num2str(getPatternDensity(patterns{m}{end}))]);
    colormap(fibroclr);
    saveas(gcf, fullfile(fibrosis_pattern_folder, ['fibrosis_pattern_', num2str(m), '.png']));
    close(fig);
end

fprintf('Patterns generated successfully! Files saved in the %s folder.\n', output_folder);

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
