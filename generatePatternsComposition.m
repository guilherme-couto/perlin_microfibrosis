function patterns = generatePatternsComposition(params, density, N_patterns, tolerance, mesh)
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

% Create a folder to save the images
output_folder = 'generated_patterns_composition';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Variable to store each pattern used to generate the desired patterns
patterns = cell(N_patterns, 1);
seeds = cell(N_patterns, 1);

% Initialise a figure to plot some of the patterns
figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
num_plots_y = round(sqrt(N_patterns));
num_plots_x = ceil(N_patterns/num_plots_y);

% Create the requested number of patterns
for m = 1:N_patterns

    % Initialise counters
    initial_seed = m;
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
    if density > 0.1
        threshold = 0.1;
    elseif density > 0.05
        threshold = 0.05;
    else
        threshold = 0.01;
    end
    seed = getFeasibleSeed(initial_seed);
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibresModified(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
    end

    % Get the density of the pattern generated
    actual_density = getPatternDensity(presence);

    % Save the pattern's presence and seed
    patterns{m}{end+1} = presence;
    seeds{m}{end+1} = seed;

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
        % seed = randi([1, 250]);
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

            % Use the fibre-free generator if NaNs are present in input params
            if density > 0.1
                threshold = 0.1;
            elseif density > 0.05
                threshold = 0.05;
            else
                threshold = 0.01;
            end
            if any(isnan(params))
                [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_tables{seed}, offset_tables{seed});
            elseif  length(params) == 6
                [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
            else
                [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_tables{seed}, offset_tables{seed});
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

% Save the figure (by Guilherme)
saveas(gcf, fullfile(output_folder, 'patterns.png'));

% Open file to save seeds used to generate the patterns
seed_file = fullfile(output_folder, 'pattern_seeds.txt');
fid = fopen(seed_file, 'w');
fprintf(fid, 'Pattern Number | Seed Used\n');

% Save the presences and the seeds used to generate the patterns
for m = 1:N_patterns
    % GIF filename
    gif_filename = fullfile(output_folder, ['pattern_', num2str(m), '.gif']);

    % Loop para salvar os frames do GIF
    for n = 1:length(patterns{m})
        % Criar a figura do padrão atual
        fig = figure('Visible', 'off');
        imagesc(patterns{m}{n});
        axis('equal', 'off');
        title(['Pattern ', num2str(m), ' - Density: ', num2str(getPatternDensity(patterns{m}{n}))]);
        colormap(fibroclr);
        saveas(gcf, fullfile(output_folder, ['pattern_', num2str(m), '_presence_', num2str(n), '.png']));
        
        % Capturar o frame da imagem
        img = imread(fullfile(output_folder, ['pattern_', num2str(m), '_presence_', num2str(n), '.png']));
        
        % Salvar no GIF (se for o primeiro frame, cria o arquivo, senão adiciona)
        if n == 1
            imwrite(img, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
            colormap(fibroclr);
        else
            imwrite(img, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
            colormap(fibroclr);
        end
        
        % Fechar a figura
        close(fig);
    end
    close(gif_filename);
    % Save the seeds used to generate the patterns
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
end

% Close the file
fclose(fid);

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