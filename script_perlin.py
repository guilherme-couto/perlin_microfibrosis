import subprocess
import sys
import numpy as np
import imageio.v2

# This script generates fibrosis patterns using the Perlin noise generator
def execute_fibrosis_generator(params, param_name, density, N_patterns, fibro_typename):
    command = f"octave --no-gui -q -f ./script_with_args.m -- '{params}' {param_name} {density} {N_patterns} {fibro_typename}"
    print(f'Running command {command}')
    subprocess.run(command, shell=True)

# PARAMETER INFORMATION:
# 1 - FIBRENESS (f): The extent to which patterns exhibit long, thin fibres
#     aligned in consistent directions
#     ::: If set to NaN, a pattern without fibres will be created :::
# 2 - FIBRE SEPARATION (L): The average spacing between fibres (in units
#     matching those used in input mesh)
# 3 - PATCHINESS (d): The extent of inconsistency in pattern density (high
#     patchiness will produce distinct regions of higher and lower density)
# 4 - FEATURE SIZE (lb): The overall size of obstacle features in the mesh (in
#     units matching the mesh)
# 5 - ROUGHNESS (gamma): The roughness of feature edges (values from [0,1], may
#     cause issues for values of precisely 0 or 1)
# 6 - PATCH SIZE (ld): The size of regions over which density varies (with
#     extent according to PATCHINESS)
# 7 - FIBRE ALIGNMENT (R): The extent to which non-fibre features are aligned
#     to the fibre direction (i.e. extent of feature anisotropy)
# 8 - DIRECTION (phi): An angle (in radians) defining the orientation of fibres
#     and/or feature anisotropy

params_names_ranges = { 'fibreness': [0.0, 0.4],
                        'fibre_separation': [0.3, 2.0],
                        'patchiness': [0.0, 0.5],
                        'feature_size': [0.01, 2.0],
                        'roughness': [0.0, 0.99],
                        'patch_size': [1, 8],
                        'fibre_alignment': [0.5, 50.0],
                        'direction': [-np.pi/2, np.pi/2]}

f_paper =       {'interstitial': 0.3,     'compact': 'NaN',    'diffuse': 'NaN',   'patchy': 0.38}
L_paper =       {'interstitial': 0.31,    'compact': 'NaN',    'diffuse': 'NaN',   'patchy': 0.31}
d_paper =       {'interstitial': 0.32,    'compact': 0.44,     'diffuse': 0.49,    'patchy': 0.42}
lb_paper =      {'interstitial': 0.24,    'compact': 0.96,     'diffuse': 0.07,    'patchy': 0.32}
gamma_paper =   {'interstitial': 0.96,    'compact': 0.59,     'diffuse': 0.44,    'patchy': 0.78}
ld_paper =      {'interstitial': 4.67,    'compact': 2.03,     'diffuse': 1.22,    'patchy': 2.1}
R_paper =       {'interstitial': 1.89,    'compact': 2.47,     'diffuse': 2.17,    'patchy': 2.5}
phi_paper =     {'interstitial': 0.59341, 'compact': -0.15708, 'diffuse': 0.19199, 'patchy': 1.1868}
density_paper = {'interstitial': 0.096,   'compact': 0.472,    'diffuse': 0.220,   'patchy': 0.269}

num_param_samples = 10

num_patterns = 10

fibros_typenames = ['interstitial', 'compact', 'diffuse', 'patchy']
for fibro_typename in fibros_typenames:
    for param_name, param_range in params_names_ranges.items():
        param_samples = np.linspace(param_range[0], param_range[1], num_param_samples)
        for m, param_value in enumerate(param_samples):
            
            d = d_paper[fibro_typename]
            L = L_paper[fibro_typename]
            f = f_paper[fibro_typename]
            lb = lb_paper[fibro_typename]
            gamma = gamma_paper[fibro_typename]
            ld = ld_paper[fibro_typename]
            R = R_paper[fibro_typename]
            phi = phi_paper[fibro_typename]
            density = density_paper[fibro_typename]

            if param_name == 'fibreness':
                f = param_value
            elif param_name == 'fibre_separation':
                L = param_value
            elif param_name == 'patchiness':
                d = param_value
            elif param_name == 'feature_size':
                lb = param_value
            elif param_name == 'roughness':
                gamma = param_value
            elif param_name == 'patch_size':
                ld = param_value
            elif param_name == 'fibre_alignment':
                R = param_value
            elif param_name == 'direction':
                phi = param_value
            else:
                print(f'Error: parameter unknown {param_name}')
                sys.exit(1)

            params = f'[{f}, {L}, {d}, {lb}, {gamma}, {ld}, {R}, {phi}]'
            execute_fibrosis_generator(params, param_name, density, num_patterns, fibro_typename)

        # Generate a gif fixing the pattern number and changing the param value 
        for m in range(num_patterns):
            images = []
            for param_value in param_samples:
                images.append(imageio.v2.imread(f'./patterns/{fibro_typename}/{param_name}/samples/{(param_value):.2f}_{m+1}.png'))
            imageio.v2.mimsave(f'./patterns/{fibro_typename}/{param_name}/pattern_{m+1}.gif', images, duration=0.5)

        # Generate a gif fixing the param value and changing the pattern number
        for param_value in param_samples:
            images = []
            for m in range(num_patterns):
                images.append(imageio.v2.imread(f'./patterns/{fibro_typename}/{param_name}/samples/{(param_value):.2f}_{m+1}.png'))
            imageio.v2.mimsave(f'./patterns/{fibro_typename}/{param_name}/param_{(param_value):.2f}.gif', images, duration=0.5)
        
         
