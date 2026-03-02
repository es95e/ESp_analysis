import pandas as pd
import os
import glob
import matplotlib.pyplot as plt

INPUT_DIR = "output_folder"

AGGREGATE_OUTPUT_FILE = "Biological_Anomaly_Counts_Trends.csv"
GRAPH_OUTPUT_FILE = "Biological_Anomaly_Trends.png"

def run_anomaly_frequency_analysis_and_plot():
    print("Anomaly frequency analysis")
    
    search_pattern = os.path.join(INPUT_DIR, "Biological_Profile_P*_*.csv")
    
    biological_files = sorted(glob.glob(search_pattern), 
                              key=lambda x: int(os.path.basename(x).split('_P')[-1].split('_')[1].split('.')[0]), 
                              reverse=False)
    
    if not biological_files:
        print("Missing proper file 'Biological_Profile_PXX_YY.csv' in 'output_folder'.")
        print("Please respect pipeline order.")
        return

    results = []
    
    try:
        df_example = pd.read_csv(biological_files[0])
        param_columns = [col for col in df_example.columns if col not in ['gene_id', 'Promoter_Signature_Count']]
        total_genes = len(df_example)
        print(f"Genes within the database: {total_genes}")
    except Exception as e:
        print(f"Error in the uploading the data: {e}")
        return

    for file_path in biological_files:
        file_name = os.path.basename(file_path)
        try:
            p_low, p_high = file_name.split('_P')[-1].split('.')[0].split('_')
            p_low, p_high = int(p_low), int(p_high)
        except:
            continue

        df_binary = pd.read_csv(file_path)
        
        anomaly_counts = df_binary[param_columns].sum().to_dict()
        total_signatures = df_binary['Promoter_Signature_Count'].sum()
        
        result_row = {
            'P_LOW': p_low,
            'P_HIGH': p_high,
            'Total_Genes': total_genes,
            'Total_Signatures_All_Geni': total_signatures,
            **anomaly_counts
        }
        results.append(result_row)
        
        print(f"[{p_low}-{p_high}]: Total Signatures: {total_signatures}")

    if not results:
        print("No significant data found")
        return

    df_trends = pd.DataFrame(results)
    
    df_trends.sort_values(by='P_HIGH', ascending=True, inplace=True) 
    
    df_trends.to_csv(AGGREGATE_OUTPUT_FILE, index=False)
    print(f"Final table in: {AGGREGATE_OUTPUT_FILE}")

    x_labels = [f"P{p}" for p in df_trends['P_HIGH']]
    
    plt.figure(figsize=(12, 8))
    
    for param in param_columns:
        plt.plot(x_labels, df_trends[param], marker='o', label=param)
        
    plt.title('Trend of Absolute Anomaly Count by Percentile Threshold')
    plt.xlabel(f"Selectivity threshold (P_HIGH, P_LOW=100-P_HIGH)")
    plt.ylabel("Number of anomalous genes ")
    plt.xticks(rotation=45)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', title="Parametri")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    plt.savefig(GRAPH_OUTPUT_FILE)
    
    print(f"Trend graph: {GRAPH_OUTPUT_FILE}")


if __name__ == "__main__":

    try:
        import matplotlib.pyplot as plt
        run_anomaly_frequency_analysis_and_plot()
    except ImportError:
        print("missing conda libraby (matplotlib). Please follow the .yml file.")
