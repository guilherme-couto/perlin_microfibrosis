import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

def analyze_csvs(directory, output_txt="general_stats_output.txt"):
    file_paths = glob.glob(os.path.join(directory, "*.csv"))
    
    # Ler todos os arquivos e concatená-los em um único DataFrame
    dfs = [pd.read_csv(f) for f in file_paths if os.stat(f).st_size > 0]
    if dfs:
        df = pd.concat(dfs, ignore_index=True)
    else:
        df = pd.DataFrame()
    
    metrics = ['major_axis', 'minor_axis', 'orientation']
    stats = {}
    
    # Calcular estatísticas
    for metric in metrics:
        values = df[metric]
        stats[metric] = {
            'mean': np.mean(values),
            'std': np.std(values, ddof=1),
            'cv': np.std(values, ddof=1) / np.mean(values),
            'min': np.min(values),
            'max': np.max(values),
            'median': np.median(values),
            'q1': np.percentile(values, 25),
            'q3': np.percentile(values, 75),
            'iqr': np.percentile(values, 75) - np.percentile(values, 25)
        }
    
    with open(output_txt, 'w') as f:
        for metric in metrics:
            f.write(f"Statistics for {metric}:\n")
            f.write(f"Mean: {stats[metric]['mean']}\n")
            f.write(f"Standard Deviation: {stats[metric]['std']}\n")
            f.write(f"Coefficient of Variation: {stats[metric]['cv']}\n")
            f.write(f"Minimum: {stats[metric]['min']}\n")
            f.write(f"Maximum: {stats[metric]['max']}\n")
            f.write(f"Median: {stats[metric]['median']}\n")
            f.write(f"Q1 (25%): {stats[metric]['q1']}\n")
            f.write(f"Q3 (75%): {stats[metric]['q3']}\n")
            f.write(f"IQR: {stats[metric]['iqr']}\n\n")
    
    # Gerar gráficos
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    for i, metric in enumerate(metrics):
        values = df[metric]
        
        # Histograma
        sns.histplot(values, bins=10, kde=True, ax=axes[i, 0], color='blue', edgecolor='black', alpha=0.7)
        axes[i, 0].axvline(stats[metric]['mean'], color='red', linestyle='dashed', linewidth=1, label='Média')
        axes[i, 0].axvline(stats[metric]['median'], color='blue', linestyle='dashed', linewidth=1, label='Mediana')
        axes[i, 0].axvline(stats[metric]['q1'], color='green', linestyle='dashed', linewidth=1, label='Q1 (25%)')
        axes[i, 0].axvline(stats[metric]['q3'], color='purple', linestyle='dashed', linewidth=1, label='Q3 (75%)')
        axes[i, 0].set_title(f"Distribuição de {metric}")
        axes[i, 0].set_xlabel(metric)
        axes[i, 0].set_ylabel("Frequência")
        axes[i, 0].legend()
        
        # Boxplot
        sns.boxplot(x=values, ax=axes[i, 1], color='lightblue')
        axes[i, 1].set_title(f"Boxplot de {metric} - CV: {stats[metric]['cv']:.2f}")
        axes[i, 1].set_xlabel(metric)
    
    plt.tight_layout()
    plt.savefig("metrics_analysis.png")

def read_uncertainty_stats(directory):
    file_paths = glob.glob(os.path.join(directory, "uncertainty_stats_*.txt"))    
    file_paths = sorted(file_paths, key=lambda x: int(x.split('_')[-1].split('.')[0]))
    stats_history = []
    power_thresholds = set()
    samples_steps = []
    
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        num_samples = int(file_path.split('_')[-1].split('.')[0])
        samples_steps.append(num_samples)
        
        stats = {}
        current_threshold = None
        current_metric = None
        
        for line in lines:
            line = line.strip()
            if line.startswith("Threshold"):
                current_threshold = int((line.split()[1][:-1]).split('%')[0])
                power_thresholds.add(current_threshold)
                stats[current_threshold] = {}
            elif line.startswith("Statistics for"):
                current_metric = line.split()[2][:-1]
                stats[current_threshold][current_metric] = {}
            elif current_metric and ':' in line:
                key, value = line.split(':')
                stats[current_threshold][current_metric][key.strip()] = float(value.strip())
        
        stats_history.append(stats)
    
    return stats_history, sorted(power_thresholds), sorted(samples_steps)

def plot_convergence(stats_history, power_thresholds, samples_steps, stat):
    metric_labels = ['major_axis', 'minor_axis', 'orientation']
    num_metrics = len(metric_labels)
    fig, axs = plt.subplots(num_metrics, 1, figsize=(12, 4*num_metrics))
    
    for row, metric in enumerate(metric_labels):
        for threshold in power_thresholds:
            values = [stats[threshold][metric][stat] for stats in stats_history]
            axs[row].plot(samples_steps, values, marker='.', label=f'Ellipse {threshold}%')
            axs[row].set_title(f'{metric}', fontweight='bold')
            axs[row].set_xlabel('Number of Samples')
            axs[row].set_ylabel(f'{stat.capitalize()} Value')
            axs[row].legend()
    
    plt.suptitle(f'{stat.capitalize()} Convergence for Ellipse Metrics', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{stat}_convergence_plot.png')
    plt.close()

def main():
    analyze_csvs("ellipses_metrics")
    stats_history, power_thresholds, samples_steps = read_uncertainty_stats("uncertainty_stats")
    plot_convergence(stats_history, power_thresholds, samples_steps, 'mean')
    plot_convergence(stats_history, power_thresholds, samples_steps, 'std')

    
if __name__ == "__main__":
    main()
    