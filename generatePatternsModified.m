function patterns = generatePatternsModified(params, param_name, density, N_patterns, fibro_typename, mesh)
% This function takes the provided set of parameter values, and creates the
% requested number of representative patterns.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% param_name - the name of the parameter being varied
% density - density of fibrosis in the patterns to be generated
% N_patterns - the number of patterns to generate
% fibro_typename - the name of the fibrosis type being generated
% (mesh) - optionally provided mesh to specify size of patterns
%
% PARAMETER INFORMATION:
%
% 1 - FIBRENESS (f): The extent to which patterns exhibit long, thin fibres
%     aligned in consistent directions
%     ::: If set to NaN, a pattern without fibres will be created :::
% 2 - FIBRE SEPARATION (L): The average spacing between fibres (in units
%     matching those used in input mesh
% 3 - PATCHINESS (d): The extent of inconsistency in pattern density (high
%     patchiness will produce distinct regions of higher and lower density)
% 4 - FEATURE SIZE (lb): The overall size of obstacle features in the mesh (in
%     units matching the mesh)
% 5 - ROUGHNESS (gamma): The roughness of feature edges (values from [0,1], may
%     cause issues for values of precisely 0 or 1)
% 6 - PATCH SIZE (ld): The size of regions over which density varies (with
%     extent according to PATCHINESS)
% 7 - FIBRE ALIGNMENT (R): The extent to which non-fibre features are aligned
%     to the fibre direction (i.e. extent of feature anisotropy)
% 8 - DIRECTION (phi): An angle (in radians) defining the orientation of fibres
%     and/or feature anisotropy

% Load in the seed data
load('fibro_seedinfo.mat', 'permute_tables', 'offset_tables');

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 6
    mesh = buildMesh(250, 400, 1/136);
end

% Get the value of the parameter being varied
param_val = params(strcmp(param_name, {'fibreness', 'fibre_separation', 'patchiness', 'feature_size', 'roughness', 'patch_size', 'fibre_alignment', 'direction'}));
if param_name == 'threshold'
    param_val = density;
end

% Create the requested number of patterns
for m = 1:N_patterns
    
    % Use the fibre-free generator if NaNs are present in input params
    % vector, or if only non-fibre parameters provided, otherwise 
    % use the standard generator
    if any(isnan(params))
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, params(3:8), permute_tables{m}, offset_tables{m});
    elseif  length(params) == 6
        [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, params, permute_tables{m}, offset_tables{m});
    else
        [presence, ~, ~, ~] = createFibroPattern(mesh, density, params, permute_tables{m}, offset_tables{m});
    end
    
    % Store this pattern
    patterns{m} = presence;

    % Save the pattern as a PNG image silently
    fig = figure('visible', 'off');
    imagesc(presence);
    axis('off', 'image');
    colormap(fibroclr);
    title(sprintf('%s = %0.2f m = %d', param_name, param_val, m), 'Interpreter', 'none');
    % Adjust the figure size to remove white space
    set(gcf, 'Position', [100, 100, size(presence, 2), size(presence, 1)]);
    filename = sprintf('./patterns/%s/%s/samples/%0.2f_%d.png', fibro_typename, param_name, param_val, m);
    check_and_create_dirs(filename);
    print(filename, '-dpng', '-r300');
    close(fig);

end

end