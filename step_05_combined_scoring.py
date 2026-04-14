# Copyright (c) 2026 Emanuele Selleri - University of Parma
# Licensed under the MIT License
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore

OUTPUT_DIR = "output_folder"
BASIC_MEAN_FILE = os.path.join(OUTPUT_DIR, "BasicMeanOnly.csv")
SCORED_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "Continuous_Combined_Scoring.csv")

WEIGHTS_A = {
    "dG_binding": 3,
    "MGW": 2,
    "HelT": 2,
    "Shift": 1.5,
}

CLUSTER_GROUPS_B = {
    'Torsion and Rotation': ["Twist", "Propeller", "HelT"],
    'Translation and Flexion': ["Slide", "Rise", "Roll", "Tilt"],
    'Energy and Stability': ["dG_binding", "dH", "dS", "dG_stacked"],
    'Minor Groove': ["MGW", "ProT"],
    'Premises and Opening': ["Shift", "Buckle", "Opening", "Stretch_BP", "Stagger", "Shear"],
}

PARAM_DIRECTION = {
    "Shift": 'LOW', "Slide": 'BOTH', "Rise": 'LOW',
    "Tilt": 'BOTH', "Roll": 'BOTH', "Twist": 'LOW',
    "Propeller": 'BOTH', "Buckle": 'BOTH', "Opening": 'LOW',
    "Stretch_BP": 'BOTH', "Stagger": 'BOTH', "Shear": 'BOTH',
    "dG_binding": 'HIGH', "dH": 'HIGH', "dS": 'BOTH',
    "dG_stacked": 'HIGH',
    "MGW": 'LOW', "ProT": 'BOTH', "HelT": 'BOTH'
}

def calculate_continuous_scores(df_means):
    """
    Calculates the continuous score based on the directional Z-score.
    """
    param_columns = df_means.columns.tolist()
    
    df_zscore = df_means.apply(zscore)
    
    df_directed_score = pd.DataFrame(index=df_means.index)
    
    for param in param_columns:
        direction = PARAM_DIRECTION.get(param, 'BOTH')
        z_scores = df_zscore[param]
        
        directed_score = pd.Series(0.0, index=df_means.index)

        if direction == 'HIGH':
            directed_score = z_scores.clip(lower=0) 
            
        elif direction == 'LOW':
            directed_score = (-z_scores).clip(lower=0)
            
        elif direction == 'BOTH':
            directed_score = z_scores.abs()
        
        df_directed_score[param] = directed_score

    df_directed_score_norm = df_directed_score.div(df_directed_score.max(axis=0), axis=1).fillna(0)

    score_C_base = df_directed_score_norm.sum(axis=1)

    score_A_weighted = pd.Series(0.0, index=df_means.index)
    for param in param_columns:
        weight = WEIGHTS_A.get(param, 1.0)
        score_A_weighted += df_directed_score_norm[param] * weight
        
    score_B_cluster = pd.Series(0.0, index=df_means.index)
    
    for cluster_name, params in CLUSTER_GROUPS_B.items():
        existing_params = [p for p in params if p in df_directed_score_norm.columns]
        if existing_params:
            cluster_contribution = df_directed_score_norm[existing_params].mean(axis=1).fillna(0)
            score_B_cluster += cluster_contribution

    scores = {
        'Count_A_Weighted': score_A_weighted,
        'Count_B_Cluster': score_B_cluster,
        'Count_C_ZScore': score_C_base
    }
    
    df_out = pd.DataFrame({
        'Gene_ID': df_means.index,
        'Score_A_Weighted': (scores['Count_A_Weighted'] / scores['Count_A_Weighted'].max()) * 10,
        'Score_B_Cluster_Compressed': (scores['Count_B_Cluster'] / scores['Count_B_Cluster'].max()) * 10,
        'Score_C_ZScore_Base': (scores['Count_C_ZScore'] / scores['Count_C_ZScore'].max()) * 10,
    })
    
    df_out['Final_Combined_Score'] = (df_out['Score_A_Weighted'] + df_out['Score_B_Cluster_Compressed'] + df_out['Score_C_ZScore_Base']) / 3

    return df_out.sort_values(by='Final_Combined_Score', ascending=False)

def run_combined_continuous_scoring():
    print("Analyzing the continuous score (Indexes A, B and C)")
    
    try:
        df_means = pd.read_csv(BASIC_MEAN_FILE)
        df_means.set_index("gene_id", inplace=True)
        param_columns = [col for col in df_means.columns if col not in ['gene_id']]
        df_means = df_means[param_columns]
    except FileNotFoundError:
        print(f"missing {BASIC_MEAN_FILE} file. Please follow pipeline order.")
        return

    print("Calculation of continuous scores based on directional Z-score")
    
    df_comparison = calculate_continuous_scores(df_means)
    
    df_comparison.to_csv(SCORED_OUTPUT_FILE, index=False)
    print(f"Index confrontation A, B and C reported in: {SCORED_OUTPUT_FILE}")
    
    print(df_comparison.head(10).to_string(index=False))

    print("Analysis end")

if __name__ == "__main__":
    run_combined_continuous_scoring()
