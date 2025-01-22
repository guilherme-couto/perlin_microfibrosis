function patterns = generatePatternsTests(params, density, N_patterns, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
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
if nargin < 4
    mesh = buildMesh(250, 400, 1/136);
end

% Initialise a figure to plot some of the patterns
figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Create the requested number of patterns
for m = 1:N_patterns
    
    % Use the fibre-free generator if NaNs are present in input params
    % vector, or if only non-fibre parameters provided, otherwise 
    % use the standard generator
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibresModified(mesh, density, params(3:8), permute_tables{m}, offset_tables{m});
        [presence2, ~, ~] = createFibroPatternNoFibresModified(mesh, density, params(3:8), permute_tables{randi([1, 250])}, offset_tables{randi([1, 250])});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, params, permute_tables{m}, offset_tables{m});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, density, params, permute_tables{m}, offset_tables{m});
        [presence2, ~, ~, ~] = createFibroPattern(mesh, density, params, permute_tables{randi([1, 250])}, offset_tables{randi([1, 250])});
    end
    
    % Store this pattern
    % patterns{m} = presence;

    % Count how many elements are in the presence, specifying the number of 1 and 0
    total_elements = numel(presence);
    count_1 = sum(presence(:) == 1);
    count_0 = sum(presence(:) == 0);
    density = count_1/total_elements;
    disp(['Pattern ', num2str(m), ' -> Total elements: ', num2str(total_elements), ' | 1s: ', num2str(count_1), ' | 0s: ', num2str(count_0), ' | Density: ', num2str(density)]);

    % Count how many elements are in the presence2, specifying the number of 1 and 0
    total_elements2 = numel(presence2);
    count_1_2 = sum(presence2(:) == 1);
    count_0_2 = sum(presence2(:) == 0);
    density2 = count_1_2/total_elements2;
    disp(['Pattern2 ', num2str(m), ' -> Total elements: ', num2str(total_elements2), ' | 1s: ', num2str(count_1_2), ' | 0s: ', num2str(count_0_2), ' | Density: ', num2str(density2)]);

    % Sum and analyze the somatory of the two patterns
    sum_presences = presence + presence2;
    sum_presences(sum_presences == 2) = 1; % if element is 2, it means that both patterns have 1 in that position, so it is a 1
    total_elements_sum = numel(sum_presences);
    count_1_sum = sum(sum_presences(:) == 1);
    count_0_sum = sum(sum_presences(:) == 0);
    density_sum = count_1_sum/total_elements_sum;
    disp(['SumPattern ', num2str(m), ' -> Total elements: ', num2str(total_elements_sum), ' | 1s: ', num2str(count_1_sum), ' | 0s: ', num2str(count_0_sum), ' | Density: ', num2str(density_sum)]);

    % If m is less than 10, plot the patterns
    if m <= 10
        % Plot this pattern
        subplot(3,3,m);
        imagesc(presence);
        axis('equal', 'off');
        title(['Density: ', num2str(density)]);
        colormap(fibroclr);
        
        % Plot this pattern
        subplot(3,3,m+3);
        imagesc(presence2);
        axis('equal', 'off');
        title(['Density: ', num2str(density2)]);
        colormap(fibroclr);

        % Plot this pattern
        subplot(3,3,m+6);
        imagesc(sum_presences);
        axis('equal', 'off');
        title(['Density: ', num2str(density_sum)]);
        colormap(fibroclr);
    end

end

% % Save the figure (by Guilherme)
saveas(gcf, 'pattern_examples_somatory.png');

end