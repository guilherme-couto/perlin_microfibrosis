import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_csv(file_path, output_txt="stats_output.txt"):
    # Ler o arquivo CSV
    df = pd.read_csv(file_path)
    
    # Selecionar as colunas de interesse
    metrics = ['major_axis', 'minor_axis']
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
    
    # Escrever estatísticas em um arquivo TXT
    with open(output_txt, "w") as f:
        for metric, values in stats.items():
            f.write(f"Statistics for {metric}:\n")
            for stat_name, stat_value in values.items():
                f.write(f"  {stat_name}: {stat_value:.4f}\n")
            f.write("\n")
    
    print(f"Estatísticas salvas em {output_txt}")
    
    # Gerar gráficos
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
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
        axes[i, 1].set_title(f"Boxplot de {metric}")
        axes[i, 1].set_xlabel(metric)
    
    plt.tight_layout()
    plt.savefig("metrics_analysis.png")

def main():
    analyze_csv("ellipses_metrics.csv")
    
if __name__ == "__main__":
    main()
