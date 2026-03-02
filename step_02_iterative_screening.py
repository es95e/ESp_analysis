import pandas as pd
import numpy as np
import os

OUTPUT_DIR = "output_folder"
BASIC_MEAN_FILE = os.path.join(OUTPUT_DIR, "BasicMeanOnly.csv")

PARAM_DIRECTION = {
    "Shift": 'LOW',     "Slide": 'BOTH',    "Rise": 'LOW',  
    "Tilt": 'BOTH',     "Roll": 'BOTH',     "Twist": 'LOW',
    "Propeller": 'BOTH', "Buckle": 'BOTH',  "Opening": 'LOW',  
    "Stretch_BP": 'BOTH', "Stagger": 'BOTH', "Shear": 'BOTH',  

    "dG_binding": 'HIGH', "dH": 'HIGH',      "dS": 'BOTH',  
    "dG_stacked": 'HIGH',
    
    "MGW": 'LOW',         "ProT": 'BOTH',     "HelT": 'BOTH'
}

def run_iterative_screening():
    
    print("--- Setting screening paramethers ---")
    
    while True:
        try:
            p_range_input = input("Percentile threshold (example 70-95): ")
            p_start, p_end = map(int, p_range_input.split('-'))
            if not (1 <= p_start < p_end <= 99):
                raise ValueError
            break
        except:
            print("Invalid or missing format (correct example 70-95).")

    while True:
        try:
            p_step = int(input("Increment step (example. 5): "))
            if p_step <= 0:
                raise ValueError
            break
        except:
            print("Invalid number (expected a integer, positive number.")

    high_percentiles = list(range(p_end, p_start - p_step, -p_step))
    
    print(f"Starting the screening with threshold: {high_percentiles}")
    
    try:
        df = pd.read_csv(BASIC_MEAN_FILE)
        df.set_index("gene_id", inplace=True)
    except FileNotFoundError:
        print(f"Missing {BASIC_MEAN_FILE}.")
        print("The scirpt must be execute inside the folder containing the subfolder 'output_folder'.")
        return

    param_columns = df.columns.tolist()

    for p_high in high_percentiles:
        p_low = 100 - p_high 
        
        print(f"Execution for LOW={p_low}% / HIGH={p_high}% ***")

        percentile_low_values = df.quantile(p_low / 100)
        percentile_high_values = df.quantile(p_high / 100)
        

        df_general_binary = pd.DataFrame(index=df.index)

        for param in param_columns:

            col_low_name = f"{param}_LOW"
            df_general_binary[col_low_name] = (df[param] <= percentile_low_values[param]).astype(int)
            
            col_high_name = f"{param}_HIGH"
            df_general_binary[col_high_name] = (df[param] >= percentile_high_values[param]).astype(int)

        df_general_binary['Total_General_Outliers'] = df_general_binary.sum(axis=1)
        df_general_binary.sort_values(by='Total_General_Outliers', ascending=False, inplace=True)

        df_general_binary.reset_index(inplace=True)
        cols = df_general_binary.columns.tolist()
        cols.insert(1, cols.pop(cols.index('Total_General_Outliers'))) 
        df_general_binary = df_general_binary[cols]

        general_output_name = os.path.join(OUTPUT_DIR, f"General_Profile_P{p_low}_{p_high}.csv")
        df_general_binary.to_csv(general_output_name, index=False)
        print(f"General profiles saved in : {general_output_name}")

        df_biological_binary = pd.DataFrame(index=df.index)

        for param in param_columns:
            direction = PARAM_DIRECTION.get(param, 'BOTH')
            is_outlier = pd.Series(False, index=df.index)
            
            if direction == 'HIGH':
                is_outlier = (df[param] >= percentile_high_values[param])
            elif direction == 'LOW':
                is_outlier = (df[param] <= percentile_low_values[param])
            elif direction == 'BOTH':
                is_outlier = (df[param] <= percentile_low_values[param]) | \
                             (df[param] >= percentile_high_values[param])
            
            df_biological_binary[param] = is_outlier.astype(int)

        df_biological_binary['Promoter_Signature_Count'] = df_biological_binary.sum(axis=1)

        df_biological_binary.sort_values(by='Promoter_Signature_Count', ascending=False, inplace=True)
        df_biological_binary.reset_index(inplace=True)
        cols = df_biological_binary.columns.tolist()
        cols.insert(1, cols.pop(cols.index('Promoter_Signature_Count')))
        df_biological_binary = df_biological_binary[cols]

        biological_output_name = os.path.join(OUTPUT_DIR, f"Biological_Profile_P{p_low}_{p_high}.csv")
        df_biological_binary.to_csv(biological_output_name, index=False)
        print(f" Biological profiles saved in : {biological_output_name}")

        MIN_SIGNATURE_COUNT = 3 # Manteniamo la soglia di 3 per la stampa a video
        strong_candidates = df_biological_binary[df_biological_binary['Promoter_Signature_Count'] >= MIN_SIGNATURE_COUNT]
        
        print(f" {len(strong_candidates)} putative candidates found")

if __name__ == "__main__":
    run_iterative_screening()
