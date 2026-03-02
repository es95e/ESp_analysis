#!/usr/bin/env python3
import csv
from Bio import SeqIO
import numpy as np
import os

UPSTREAM = 200
DOWNSTREAM = 70 

PARAMETERS = [
    # Structural parameters from Olson2001, reported in tab1.csv)
    "Shift", "Slide", "Rise", "Tilt",
    "Roll", "Twist",
    "Propeller", "Buckle", "Opening", "Stretch_BP",
    "Stagger", "Shear",
    
    # Energetic / Thermodynamic parameters from SantaLucia1998, reported in tab2.csv
    "dG_binding", "dH", "dS", "dG_stacked",
    
    # DNA shape feature models from Young2022 and DNAshapeR, reported in tab3.csv)
    "MGW", "ProT", "HelT"
]

GENOME_FASTA = "your_fasta_file.fasta"
GENE_GFF = "your_gff_file.gff"
TAB1_FILE = "tab1.csv"
TAB2_FILE = "tab2.csv"
TAB3_FILE = "tab3.csv"
OUTPUT_DIR = "output_folder"

def read_genome(fasta_file):
    """Upload genome"""
    genome_seq = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            genome_seq[record.id] = str(record.seq).upper()
        return genome_seq
    except FileNotFoundError:
        print(f"missing file FASTA {fasta_file}")
        exit()

def parse_gff(gff_file):
    """
    Upload gff files.
    """
    genes = []
    acceptable_features = ["gene", "CDS", "tRNA", "rRNA", "ncRNA"]
    
    try:
        with open(gff_file, "r") as f:
            for line in f:
                if line.startswith("#") or line.strip() == "":
                    continue
                
                parts = line.strip().split("\t")
                
                if len(parts) < 9 or parts[2] not in acceptable_features:
                    continue
                
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                locus_tag = None
                
                for attr in attributes.split(";"):
                    if attr.startswith("locus_tag="):
                        locus_tag = attr.replace("locus_tag=", "")
                    elif attr.startswith("ID="):
                         if not locus_tag: 
                            locus_tag = attr.replace("ID=", "")
                            
                gene_id = locus_tag if locus_tag else attributes.replace(";", "_")
                
                genes.append({
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "id": gene_id
                })
        return genes
    except FileNotFoundError:
        print(f"missing GFF file {gff_file}")
        exit()

def extract_uptake_region(genome_seq, gene, upstream=200, downstream=70):
    """Upstream regions extraction"""

    if gene["seqid"] not in genome_seq:
        return "" 
        
    seq = genome_seq[gene["seqid"]]
    
    if gene["strand"] == "+":
        start = max(0, gene["start"] - upstream - 1)
        end = min(len(seq), gene["start"] + downstream - 1)
        region = seq[start:end]
    else:
        start = max(0, gene["end"] - downstream)
        end = min(len(seq), gene["end"] + upstream)
        
        genomic_region = seq[start:end]
        
        region = genomic_region[::-1].translate(str.maketrans("ATGC", "TACG"))
        
    return region

def load_table(file, param_list):
    """Upload reference tables for ESp"""
    table = {}
    try:
        with open(file, 'r') as f:
            reader = csv.DictReader(f)

            if not all(p in reader.fieldnames for p in param_list):
                missing = [p for p in param_list if p not in reader.fieldnames]
                print(f"Missing proprieties in the reference files: {missing}.")
            
            dinucleotide_key = reader.fieldnames[0]
            
            for row in reader:
                key = row[dinucleotide_key].upper()
                table[key] = {}
                for param in param_list:
                    if param in row and row[param]:
                        try:
                            table[key][param] = float(row[param])
                        except ValueError:
                            table[key][param] = np.nan  
                    else:
                        table[key][param] = np.nan 
        return table
    except FileNotFoundError:
        print(f"Missing reference file {file}")
        exit()
    except Exception as e:
        print(f"Error in uploading reference file {file}: {e}")
        exit()

def calculate_parameter_means(tab1, tab2, tab3):
    """Average paramthers calculation for each upstream region."""

    param_map = {}
    for p in PARAMETERS[:12]: param_map[p] = tab1   
    for p in PARAMETERS[12:16]: param_map[p] = tab2  
    for p in PARAMETERS[16:]: param_map[p] = tab3   
    
    global_means = {}
    for p in PARAMETERS:
        values = []
        current_table = param_map[p]
        for step in current_table:
            if p in current_table[step] and not np.isnan(current_table[step][p]):
                values.append(current_table[step][p])
         
        global_means[p] = np.mean(values) if values else 0.0
        if np.isnan(global_means[p]):
            global_means[p] = 0.0 
            
    return global_means

def calculate_parameters_for_sequence(seq, tab1, tab2, tab3, global_means):
    param_dict = {p: [] for p in PARAMETERS}
    
    param_table_map = {}
    for p in PARAMETERS[:12]: param_table_map[p] = tab1
    for p in PARAMETERS[12:16]: param_table_map[p] = tab2
    for p in PARAMETERS[16:]: param_table_map[p] = tab3

    for i in range(len(seq)-1):
        step = seq[i:i+2]
        
        for p in PARAMETERS:
            current_table = param_table_map[p]
            step_value = None

            if step in current_table and p in current_table[step] and not np.isnan(current_table[step][p]):
                step_value = current_table[step][p]
            
            if step_value is None:
                fallback_value = global_means.get(p, 0.0)
                param_dict[p].append(fallback_value)
            else:
                param_dict[p].append(step_value)
    
    for p in PARAMETERS:
        if len(param_dict[p]) < len(seq):
            param_dict[p].append(param_dict[p][-1] if param_dict[p] else global_means.get(p, 0.0)) 
            
    return param_dict

def generate_basic_mean(param_dict):
    mean_dict = {}
    for param, values in param_dict.items():
        cleaned_values = [v for v in values if not np.isnan(v)] 
        mean_dict[param] = float(np.mean(cleaned_values)) if cleaned_values else 0.0
    return mean_dict

def generate_statistical(param_dict):
    stat_dict = {}
    for param, values in param_dict.items():
        arr = np.array([v for v in values if not np.isnan(v)])
        if arr.size > 0:
            stat_dict[param + "_mean"] = float(np.mean(arr))
            stat_dict[param + "_std"] = float(np.std(arr))
            stat_dict[param + "_min"] = float(np.min(arr))
            stat_dict[param + "_max"] = float(np.max(arr))
        else:
            stat_dict[param + "_mean"] = stat_dict[param + "_std"] = stat_dict[param + "_min"] = stat_dict[param + "_max"] = 0.0
    return stat_dict

def main():
    print("Start analysis...")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    genome_seq = read_genome(GENOME_FASTA)
    genes = parse_gff(GENE_GFF)

    print(f"Upload references files ({TAB1_FILE}, {TAB2_FILE}, {TAB3_FILE})...")

    tab1 = load_table(TAB1_FILE, PARAMETERS[:12])
    tab2 = load_table(TAB2_FILE, PARAMETERS[12:16]) 
    tab3 = load_table(TAB3_FILE, PARAMETERS[16:])
    
    global_means = calculate_parameter_means(tab1, tab2, tab3)

    basic_mean_rows = []
    statistical_rows = []
    full_vector_lines = []
    
    for gene in genes:
        region_seq = extract_uptake_region(genome_seq, gene, UPSTREAM, DOWNSTREAM)
        if not region_seq:
            print(f"Upstream region for {gene['id']} empty. Skip {gene['id']}.")
            continue
            
        params = calculate_parameters_for_sequence(region_seq, tab1, tab2, tab3, global_means)

        basic_mean = generate_basic_mean(params)
        basic_mean_row = {"gene_id": gene["id"]}
        basic_mean_row.update(basic_mean)
        basic_mean_rows.append(basic_mean_row)

        statistical = generate_statistical(params)
        statistical_row = {"gene_id": gene["id"]}
        statistical_row.update(statistical)
        statistical_rows.append(statistical_row)

        for param, values in params.items():
            line = f"{gene['id']}_{param}: {','.join([f'{v:.4f}' for v in values])}"
            full_vector_lines.append(line)

    print("Generating output files")
    
    basic_csv_file = os.path.join(OUTPUT_DIR, "BasicMeanOnly.csv")
    with open(basic_csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene_id"] + PARAMETERS)
        writer.writeheader()
        writer.writerows(basic_mean_rows)

    stat_fields = ["gene_id"] + [p + suffix for p in PARAMETERS for suffix in ["_mean","_std","_min","_max"]]
    stat_csv_file = os.path.join(OUTPUT_DIR, "Statistical.csv")
    with open(stat_csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=stat_fields)
        writer.writeheader()
        writer.writerows(statistical_rows)

    full_vector_file = os.path.join(OUTPUT_DIR, "FullVector.txt")
    with open(full_vector_file, "w") as f:
        for line in full_vector_lines:
            f.write(line + "\n")

    desc_file = os.path.join(OUTPUT_DIR, "Parameter_Description.txt")
    with open(desc_file, "w") as f:
        f.write("Parameter\tDescription\tSource\n")
        f.write("--- Olson2001 Parameters (Structural) ---\n")
        for p in PARAMETERS[:12]:
            f.write(f"{p}\tStructural parameter (Base-pair step)\tfrom {TAB1_FILE}\n")
        f.write("--- SantaLucia1998 Parameters (Energetic) ---\n")
        for p in PARAMETERS[12:16]:
            f.write(f"{p}\tEnergetic/Thermodynamic parameter\tfrom {TAB2_FILE}\n")
        f.write("--- DNAshapeR Parameters (Shape) ---\n")
        for p in PARAMETERS[16:]:
            f.write(f"{p}\tDNA Shape Feature\tfrom {TAB3_FILE}\n")


    print("--- End of the analysis ---")
    print(f"Files {OUTPUT_DIR}:\n- BasicMeanOnly.csv\n- Statistical.csv\n- FullVector.txt\n- Parameter_Description.txt")

if __name__ == "__main__":
    main()
