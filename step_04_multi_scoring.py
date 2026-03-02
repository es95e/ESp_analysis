import pandas as pd
import numpy as np
import os
from scipy.stats import zscore
from sklearn.cluster import AgglomerativeClustering

OUTPUT_DIR = "output_folder"
BASIC_MEAN_FILE = os.path.join(OUTPUT_DIR, "BasicMeanOnly.csv")
SCORED_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "Biological_Scoring_Comparison.csv")

# A) Paramethers biological weight
WEIGHTS_A = {
    "dG_binding": 3,  
    "MGW": 2,         
    "HelT": 2,        
    "Shift": 1.5,     
}

# B) Paramethers redundancy reduction 
CLUSTER_GROUPS_B = {
    'Torsione_e_Rotazione': ["Twist", "Propeller", "HelT"],
    'Translazione_e_Flessione': ["Slide", "Rise", "Roll", "Tilt"],
    'Energetici_Stabilità': ["dG_binding", "dH", "dS", "dG_stacked"],
    'Solco_Minore': ["MGW", "ProT"],
    'Locali_e_Apertura': ["Shift", "Buckle", "Opening", "Stretch_BP", "Stagger", "Shear"],
}

# C) Paramethers biological map (LOW/HIGH/BOTH)
PARAM_DIRECTION = {
    "Shift": 'LOW',     "Slide": 'BOTH',    "Rise": 'LOW', 
    "Tilt": 'BOTH',     "Roll": 'BOTH',     "Twist": 'LOW',
    "Propeller": 'BOTH', "Buckle": 'BOTH',  "Opening": 'LOW', 
    "Stretch_BP": 'BOTH', "Stagger": 'BOTH', "Shear": 'BOTH', 
    "dG_binding": 'HIGH', "dH": 'HIGH',      "dS": 'BOTH', 
    "dG_stacked": 'HIGH',
    "MGW": 'LOW',         "ProT": 'BOTH',     "HelT": 'BOTH'
}

def calculate_anomaly_score(df_means, p_low, p_high):
    """Calculation of indexes A, B and C"""
    
    p_low_float = p_low / 100
    p_high_float = p_high / 100
    
    percentile_low_values = df_means.quantile(p_low_float)
    percentile_high_values = df_means.quantile(p_high_float)
    
    param_columns = df_means.columns.tolist()
    
    df_binary = pd.DataFrame(index=df_means.index)

    for param in param_columns:
        direction = PARAM_DIRECTION.get(param, 'BOTH')
        
        is_outlier = pd.Series(False, index=df_means.index)
        
        if direction == 'HIGH':
            is_outlier = (df_means[param] >= percentile_high_values[param])
        elif direction == 'LOW':
            is_outlier = (df_means[param] <= percentile_low_values[param])
        elif direction == 'BOTH':
            is_outlier = (df_means[param] <= percentile_low_values[param]) | \
                         (df_means[param] >= percentile_high_values[param])
        
        df_binary[param] = is_outlier.astype(int)
    
    score_A = pd.Series(0, index=df_means.index, dtype=float)
    for param in param_columns:
        weight = WEIGHTS_A.get(param, 1.0)
        score_A += df_binary[param] * weight
    
    score_B = pd.Series(0, index=df_means.index, dtype=int)
    
    for cluster_name, params in CLUSTER_GROUPS_B.items():

        hit_in_cluster = (df_binary[params].sum(axis=1) >= 1).astype(int)
        score_B += hit_in_cluster

    df_zscore = df_means.apply(zscore)
    
    score_C = pd.Series(0, index=df_means.index, dtype=float)
    for param in param_columns:
        direction = PARAM_DIRECTION.get(param, 'BOTH')

        is_outlier_mask = df_binary[param] == 1
        if direction in ['LOW', 'HIGH']:
            z_score_abs = df_zscore[param].abs()
            score_C += z_score_abs * is_outlier_mask
        else: 
             z_score_abs = df_zscore[param].abs()
             score_C += z_score_abs * is_outlier_mask

    score_C = score_C / score_C.max() * 10 

    df_out = pd.DataFrame({
        'Gene_ID': df_means.index,
        'Count_A_Weighted': score_A,
        'Count_B_Cluster': score_B,
        'Count_C_ZScore': score_C
    })
    
    return df_out.sort_values(by='Count_C_ZScore', ascending=False)

def run_multi_scoring_analysis():
    print("Inizio analisi multi-scoring (A, B, C)...")
    
    try:
        df_means = pd.read_csv(BASIC_MEAN_FILE)
        df_means.set_index("gene_id", inplace=True)

        param_columns = [col for col in df_means.columns if col not in ['gene_id', 'Promoter_Signature_Count']]
        df_means = df_means[param_columns]
    except FileNotFoundError:
        print(f"ERRORE: Missing {BASIC_MEAN_FILE} file. Please follow pipeline order.")
        return

    P_LOW_TEST = 5
    P_HIGH_TEST = 95
    
    print(f"Calculating anomaly scores for the fixed threshold P{P_LOW_TEST}/{P_HIGH_TEST}...")
    
    df_comparison = calculate_anomaly_score(df_means, P_LOW_TEST, P_HIGH_TEST)

    df_comparison.to_csv(SCORED_OUTPUT_FILE, index=False)
    print(f"Indexes A, B and C reported in: {SCORED_OUTPUT_FILE}")

    print(" Top 10 putative candidates ")
    
    # Visualizza i primi 10
    print(df_comparison.head(10).to_string(index=False))

    print("Difference analysis")
    print(f"Range Score A (Pesi): Max {df_comparison['Count_A_Weighted'].max():.2f} (Theorical max: 4*3 + 3*2 + 12*1 = 30)")
    print(f"Range Score B (Cluster): Max {df_comparison['Count_B_Cluster'].max()} (Theorical max: {len(CLUSTER_GROUPS_B)})")
    print(f"Range Score C (Z-Score): Max {df_comparison['Count_C_ZScore'].max():.2f}")

if __name__ == "__main__":
    run_multi_scoring_analysis()
