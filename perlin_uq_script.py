import cv2, os, csv, subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
from skimage.measure import regionprops, label
from skimage.filters import gaussian

# Definir o colormap 'fibrosis' fora das funções
fibroclr = [[0.8, 0.2, 0.2], [0.95, 0.85, 0.55]]
fibro_cmap = LinearSegmentedColormap.from_list("fibrosis", fibroclr)

def remove_white_border(image):
    """Remove borda branca da imagem e ignora o título da imagem."""
    _, binary = cv2.threshold(image, 240, 255, cv2.THRESH_BINARY)
    binary = cv2.bitwise_not(binary)
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        contour = max(contours, key=cv2.contourArea)
        x, y, w, h = cv2.boundingRect(contour)
        return image[y:y+h, x:x+w]
    return image  

def plot_comparison(image, spectrum, spectrum_title, save_path):
    """Plota a imagem original e seu espectro de Fourier e salva a imagem gerada."""
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].imshow(image, cmap=fibro_cmap)
    axs[0].set_title("Imagem Binária Sem Borda")
    axs[0].axis('off')
    
    # Plotando o espectro de Fourier
    im1 = axs[1].imshow(np.log1p(spectrum), cmap='inferno')  # Escala log para melhor visualização
    axs[1].set_title(spectrum_title)
    axs[1].axis('off')
    
    # Adicionando o colorbar ao segundo subplot
    fig.colorbar(im1, ax=axs[1], label='Log Intensity')
    
    plt.savefig(save_path)
    plt.close()

def plot_comparison_with_ellipses(image, spectrum, ellipses, spectrum_title, save_path):
    """Plota a imagem original com as elipses e seu espectro de Fourier e salva a imagem gerada."""
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].imshow(image, cmap=fibro_cmap, alpha=0.5)
    axs[0].set_title("Imagem Binária Sem Borda")
    axs[0].axis('off')

    # Desenhar elipses com cores diferentes
    num_ellipses = len(ellipses.keys())
    colors = plt.cm.tab10.colors[:num_ellipses]
    for i, (threshold, (major, minor, angle)) in enumerate(ellipses.items()):
        ellipse = Ellipse(
            xy=(image.shape[1] // 2, image.shape[0] // 2),  # Centro da imagem
            width=major,
            height=minor,
            angle=np.degrees(angle),  # Convertendo radianos para graus
            edgecolor=colors[i],
            facecolor='none',
            linewidth=1.5,
            label=f'Elipse {threshold:.0f}%',
        )
        axs[0].add_patch(ellipse)

    axs[0].legend(loc='upper right')

    # Plotando o espectro de Fourier
    im1 = axs[1].imshow(np.log1p(spectrum), cmap='inferno')  # Escala log para melhor visualização
    axs[1].set_title(spectrum_title)
    axs[1].axis('off')
    
    # Adicionando o colorbar ao segundo subplot
    fig.colorbar(im1, ax=axs[1], label='Log Intensity')
    
    plt.savefig(save_path)
    plt.close()
    
def plot_convergence(stats_history, power_thresholds, save_path):
    """Plota os gráficos de convergência para cada métrica."""
    metric_labels = ['major_axis', 'minor_axis']
    num_metrics = len(metric_labels)
    fig, axs = plt.subplots(num_metrics, 1, figsize=(12, 8))
    
    for row, metric in enumerate(metric_labels):
        for threshold in power_thresholds:
            values = [stats[threshold][f'mean_{metric}'] for stats in stats_history]
            axs[row].plot(range(1, len(values) + 1), values, marker='.', label=f'{threshold}%')
            axs[row].set_title(f'{metric}')
            axs[row].set_xlabel('Number of Seeds')
            axs[row].set_ylabel('Mean Value')
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    
def save_comparison_plots(cropped_image, power_spectrum, smoothed_spectrum, ellipses, seed, base_folder):
    # Criar diretórios para salvar imagens
    os.makedirs(f"{base_folder}/spectra", exist_ok=True)
    os.makedirs(f"{base_folder}/ellipses", exist_ok=True)
    
    # Salvar espectros originais e suavizados
    plot_comparison(cropped_image, power_spectrum, "Espectro de Fourier (Potência)", f"{base_folder}/spectra/spectrum{seed}.png")
    plot_comparison_with_ellipses(cropped_image, smoothed_spectrum, ellipses, "Espectro Suavizado", f"{base_folder}/ellipses/ellipse{seed}.png")

def save_ellipses_metrics(seed, density, ellipses, csv_path):
    # Check if the file already exists, if not, create it and write the header, otherwise append the new data
    if not os.path.exists(csv_path):
        with open(csv_path, "w") as f:
            f.write("seed,density,power_threshold,major_axis,minor_axis,orientation\n")

    with open(csv_path, "a") as f:
        for threshold, (major, minor, angle) in ellipses.items():
            f.write(f"{seed},{density},{threshold},{major},{minor},{angle}\n")

def execute_fibrosis_generator(params, density, seed_num):
    command = f"octave --no-gui --no-window-system --silent ./script_uq.m -- '{params}' {density} {seed_num}"
    print(f'Running command {command}')
    subprocess.run(command, shell=True)

def compute_fft2(image):
    """Aplica a FFT2 e retorna o espectro de potência."""
    image = image - np.mean(image)  # Subtrai a média global da imagem
    f_transform = np.fft.fft2(image)
    f_transform = np.fft.fftshift(f_transform)
    power_spectrum = np.abs(f_transform) ** 2  # Intensidade do espectro
    return power_spectrum

def smooth_spectrum(spectrum, sigma=4):
    """Suaviza o espectro de potência com um filtro gaussiano."""
    return gaussian(spectrum, sigma=sigma)

def extract_ellipse_metrics(smoothed_spectrum, power_thresholds):
    """Calcula métricas das elipses a partir do espectro de potência suavizado."""
    
    # Calcular a potência total
    total_power = np.sum(smoothed_spectrum)

    # Ordenar valores de potência para determinar os limiares
    sorted_power = np.sort(smoothed_spectrum.ravel())[::-1]
    cumulative_power = np.cumsum(sorted_power)

    ellipses = {}

    for threshold in power_thresholds:
        # Determinar o valor de corte para cada limiar de potência
        cutoff_value = sorted_power[np.searchsorted(cumulative_power, threshold * total_power / 100)]
        
        # Criar máscara binária
        mask = smoothed_spectrum > cutoff_value
        
        # Identificar propriedades da maior região conectada
        labeled_mask = label(mask)
        props = regionprops(labeled_mask)
        
        if props:
            largest_region = max(props, key=lambda r: r.area)
            metrics = (
                largest_region.major_axis_length,
                largest_region.minor_axis_length,
                largest_region.orientation
            )
            ellipses[threshold] = metrics

    return ellipses

def is_evaluation_done(csv_path, seed, density, power_thresholds):
    """Verifica se as métricas para a semente já foram extraídas."""
    # If the file does not exist or just has the header, return False
    if not os.path.exists(csv_path) or os.stat(csv_path).st_size == 0:
        return False, power_thresholds

    missing_thresholds = power_thresholds
    with open(csv_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if (int(row['seed']) == seed and float(row['density']) == density):
                if all(float(row['power_threshold']) == threshold for threshold in power_thresholds):
                    return True, None
                else:
                    # Remove thresholds that have already been evaluated from the missing_thresholds list
                    missing_thresholds = [threshold for threshold in missing_thresholds if threshold != float(row['power_threshold'])]
    
    return False, missing_thresholds

def compute_statistics(csv_path, power_thresholds):
    # Calculate mean and standard deviation for major and minor axis for each threshold
    data = {threshold: {'major_axis': [], 'minor_axis': []} for threshold in power_thresholds}
    
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            threshold = float(row['power_threshold'])
            if threshold in power_thresholds:
                data[threshold]['major_axis'].append(float(row['major_axis']))
                data[threshold]['minor_axis'].append(float(row['minor_axis']))
    
    # Compute statistics
    stats = {}
    for threshold, metrics in data.items():
        if threshold not in stats:
            stats[threshold] = {}
        for metric in metrics.keys():
            if metric not in stats[threshold]:
                stats[threshold][metric] = {}
            values = metrics[metric]
            stats[threshold][metric]['mean'] = np.mean(values)
            stats[threshold][metric]['std'] = np.std(values, ddof=1)
            stats[threshold][metric]['cv'] = np.std(values, ddof=1) / np.mean(values)
            stats[threshold][metric]['min'] = np.min(values)
            stats[threshold][metric]['max'] = np.max(values)
            stats[threshold][metric]['median'] = np.median(values)
            stats[threshold][metric]['q1'] = np.percentile(values, 25)
            stats[threshold][metric]['q3'] = np.percentile(values, 75)
            stats[threshold][metric]['iqr'] = np.percentile(values, 75) - np.percentile(values, 25)
            
    return stats

def check_convergence(stats_history, tolerance=0.01):
    # Check if there are at least two sets of statistics to compare
    if len(stats_history) < 2:
        return False

    prev_stats = stats_history[-2]
    curr_stats = stats_history[-1]
    
    # Check if the mean and standard deviation for each metric are within the tolerance
    for threshold in prev_stats.keys():
        for metric in prev_stats[threshold].keys():
            prev_mean = prev_stats[threshold][metric]['mean']
            curr_mean = curr_stats[threshold][metric]['mean']
            prev_std = prev_stats[threshold][metric]['std']
            curr_std = curr_stats[threshold][metric]['std']
            
            # While mean or standard deviation above the tolerance, return False
            if abs(prev_mean - curr_mean) / (prev_mean + 1e-8) > tolerance:
                return False
            if abs(prev_std - curr_std) / (prev_mean + 1e-8) > tolerance:
                return False
    return True

def main():
    # Base folder
    base_folder = 'uncertainty_quantification'
    os.makedirs(base_folder, exist_ok=True)

    images_folder = f"{base_folder}/fibrosis_patterns"
    
    os.makedirs(f"{base_folder}/ellipses_metrics", exist_ok=True)
    os.makedirs(f"{base_folder}/uncertainty_stats", exist_ok=True)

    params = '[0.38, 0.31, 0.42, 0.32, 0.78, 2.1, 2.5, 1.18680]'
    density = 0.269
    power_thresholds = [80]
    save_plots = False

    total_seeds = 10000
    seeds_interval = np.arange(100, total_seeds + 1, 100)
    stats_history = []
    
    for num_seeds in seeds_interval:
        # Create file to store the ellipses metrics
        csv_path = f"{base_folder}/ellipses_metrics/ellipses_metrics_{num_seeds}.csv"
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as f:
                f.write("seed,density,power_threshold,major_axis,minor_axis,orientation\n")
        
        # Else, load the file and check if the number of seeds (lines) is equal to the current number of seeds multiplied by the number of thresholds
        else:
            with open(csv_path, "r") as f:
                num_lines = sum(1 for line in f)
            if num_lines == num_seeds * len(power_thresholds) + 1:
                print(f"File {csv_path} already contains {num_seeds} seeds evaluated for {len(power_thresholds)} thresholds. Skipping...")
                
                # Compute statistics and check convergence
                stats = compute_statistics(csv_path, power_thresholds)
                stats_history.append(stats)
                if check_convergence(stats_history):
                    print(f'Convergence reached after evaluating {num_seeds} seeds.')
                    break
                continue
            else:
                with open(csv_path, "w") as f:
                    f.write("seed,density,power_threshold,major_axis,minor_axis,orientation\n")
                
        stats_file = f"{base_folder}/uncertainty_stats/uncertainty_stats_{num_seeds}.txt"
        
        # Generate random seeds
        seeds = np.random.randint(1, 100000, num_seeds)
        
        for seed in seeds:
            # Check if the evaluation has already been done for the current seed and thresholds
            evaluation_done, missing_thresholds = is_evaluation_done(csv_path, seed, density, power_thresholds)
            if evaluation_done:
                print(f"Seed {seed} has already been evaluated for density {density}. See {csv_path}. Skipping...")
                continue
            
            # Execute fibrosis generator
            execute_fibrosis_generator(params, density, seed)
            
            # Load image
            image_path = f"{images_folder}/fibrosis_pattern_{seed}.png"
            image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
            
            # Remove pattern image to save memory
            os.remove(image_path)
            
            # Remove white border from image, compute FFT2 and smooth the spectrum
            cropped_image = remove_white_border(image)
            power_spectrum = compute_fft2(cropped_image)
            smoothed_spectrum = smooth_spectrum(power_spectrum)
            
            # Extract ellipse metrics and save them to a CSV file
            ellipses = extract_ellipse_metrics(smoothed_spectrum, missing_thresholds)
            save_ellipses_metrics(seed, density, ellipses, csv_path)

            if save_plots:
                save_comparison_plots(cropped_image, power_spectrum, smoothed_spectrum, ellipses, seed, base_folder)
                
        # Compute statistics and check convergence
        stats = compute_statistics(csv_path, power_thresholds)
        with open(stats_file, "w") as f:
            f.write(f"Uncertainty Quantification Statistics with {seed} seeds\n\n")
            for threshold, metrics in stats.items():
                f.write(f"Threshold {threshold}%:\n")
                for key, value in metrics.items():
                    f.write(f"  {key}: {value}\n")
                f.write("\n")
                
        stats_history.append(stats)
        if check_convergence(stats_history):
            print(f'Convergence reached after evaluating {num_seeds} seeds.')
            break

    # Plot convergence graphs and save statistics to a file
    plot_convergence(stats_history, power_thresholds, f"{base_folder}/convergence_plot.png")
    print(f"Statistics and convergence plot saved to {base_folder}.")
    
if __name__ == "__main__":
    main()
