import os
import numpy as np
import pandas as pd

################# Pick variables #################
name_of_test = "test11"
bebic = "BEBIC_B"
# bebic = "BEBIC_C"
exp_list = ["B3", "B4", "B5", "B6"]
# exp_list = ["C3", "C5", "C6"]
# modBase = "C"
modBase = "A"
significance_threshold = 0.05

motif = "A"

###################################################

################# Pipeline paths #################

wd = os.getcwd()
print(f"working directory: {wd}")
ad = "/external/analyses/lucy/bebic_pipeline"
print(f"analyses directory: {ad}")
output_dir = os.path.join(ad, name_of_test, "intermediate_data_files")
os.makedirs(output_dir, exist_ok=True)
plotting_dir = os.path.join(wd, "plots")

path2results = f"{ad}/{name_of_test}/modkit"
path2_dmrs = f"{path2results}/all_test"

if significance_threshold == 0.01:
    write_sig = "sig01"
    path2write_dmrdata = os.path.join(output_dir, f'sig_DMPs.{bebic}.modkit.{write_sig}.csv')
else:
    path2write_dmrdata = os.path.join(output_dir, f'sig_DMPs.{bebic}.modkit.csv')

###################################################

################# Set up stuff #################
chr_names = [
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM"
]
colour_mapping = {
    "B3_Control": "#b9647a",
    "B3_Test": "#E3C1CA",
    "B4_Control": "#c843c9",
    "B4_Test": "#E9B5E9",
    "B5_Control": "#d34143",
    "B5_Test": "#EDB3B4",
    "B6_Control": "#c059a3",
    "B6_Test": "#E6BCDA",
    "C3_Control": "#3b6ae1",
    "C3_Test": "#BEE4E1",
    "C5_Control": "#5cbab4",
    "C5_Test": "#B1C4F3",
    "C6_Control": "#6c8db3",
    "C6_Test": "#C4D1E0"
}

sample_names_list = list()
for exp in exp_list:
    for var in ["Control", "Test"]:
        sample_name = exp + "_" + var
        sample_names_list.append(sample_name)
        
dmr_ss_header=["chrom", "start position", "end position", "name", "score", "sample A counts", "sample A total", "sample B counts", "sample B total", "sample A fractions", "sample B fractions", "sample A percent modified", "sample B percent modified", "MAP-based p-value", "effect size", "balanced MAP-based p-value", "balanced effect size", "pct_a_samples", "pct_b_samples", "per-replicate p-values", "per-replicate effect sizes"]

###################################################

################# Functions #################

def top_num(dmr_df, num):
    # Sort by column: 'score' (descending)
    dmr_df = dmr_df.sort_values(['score'], ascending=[False])
    dmr_df = dmr_df.head(num)
    return dmr_df

def filter_to_significance(df, pval=0.05, effect=0.1):
    neg_effect = effect * -1
    df = df[(df['pct_a_samples'] >= 75) & (df['pct_b_samples'] >= 75)]
    df = df[df['balanced_map_pvalue'] <= pval]
    df = df[(df['balanced_effect_size'] >= effect) | (df['balanced_effect_size'] <= neg_effect)]
    return df

def clean_data(plot_data, cpg=False):
    # Drop columns: 'location', 'name' and 2 other columns
    plot_data = plot_data.drop(columns=['location', 'name', 'replicate_map_pvalues', 'replicate_effect_sizes', 'a_mod_percentages', 'b_mod_percentages'])
    
    # Derive column 'a_base' from column: 'a_counts'
    def base(a_counts):
        """
        Transform based on the following examples:
           a_counts      Output
        1: "a:33"     => "A"
        2: "h:2,m:42" => "C"
        3: "h:0,m:25" => "C"
        4: "a:18"     => "A"
        """
        if len(a_counts) - len(a_counts.replace(":", "")) == 2:
            return "C"  # has h in it
        if len(a_counts) - len(a_counts.replace(":", "")) == 1:
            return a_counts.split(":")[0].upper()
        return None
    
    def base_CpG(a_counts):
        return "CpG"
    
    if not cpg:
        plot_data.insert(4, "base", plot_data.apply(lambda row : base(row["a_counts"]), axis=1))
    else:
        plot_data.insert(4, "base", plot_data.apply(lambda row : base_CpG(row["a_counts"]), axis=1))
    
    # Derive column 'a_hydroxymethylated' from column: 'a_counts'
    def hydroxymethylated(counts):
        """
        Transform based on the following examples:
           a_counts       Output
        1: "h:1,m:104" => "1"
        2: "a:38"      => "0"
        """
        if len(counts) - len(counts.replace(":", "")) == 2:
            return int(counts[counts.find(":") + 1:counts.find(",")])
        if len(counts) - len(counts.replace(":", "")) == 1:
            return 0
        return None
    
    if not cpg:
        plot_data.insert(6, "a_hydroxymethylated", plot_data.apply(lambda row : hydroxymethylated(row["a_counts"]), axis=1))
        plot_data.insert(9, "b_hydroxymethylated", plot_data.apply(lambda row : hydroxymethylated(row["b_counts"]), axis=1))
    
    # Derive column 'a_methylated' from column: 'a_counts'
    # Transform based on the following examples:
    #    a_counts      Output
    # 1: "h:2,m:42" => "42"
    # 2: "a:18"     => "18"
    plot_data.insert(6, "a_methylated", plot_data["a_counts"].str.split(":").str[-1].astype(int))
    plot_data.insert(9, "b_methylated", plot_data["b_counts"].str.split(":").str[-1].astype(int))
    
    # Convert to number of samples
    plot_data.insert(19, "num_a_samples", (plot_data["pct_a_samples"]/25).astype(int))
    plot_data.insert(20, "num_b_samples", (plot_data["pct_b_samples"]/25).astype(int))
    
    plot_data = plot_data.drop(columns=['a_counts', 'b_counts', 'pct_a_samples', 'pct_b_samples'])
    plot_data.reset_index(drop=True, inplace=True)
    
    return plot_data


###################################################

dmr_all_chrs_dfs = dict()  
dmr_top20_dfs = dict()
dmr_top50_dfs = dict()
dmr_top100_dfs = dict()
dmr_top150_dfs = dict()
# meth_vals_df = pd.DataFrame()

print(f'{bebic}')
for chr_name in chr_names:
    # dmr_name = control + '_' + test
    path2dmr = os.path.join(path2_dmrs, f'raw_DMPs.{bebic}.{chr_name}.bed')
    
    if not os.path.exists(path2dmr):
        print(f'{path2dmr} does not exist')
        continue

    with open(path2dmr, 'r') as fd:
        try:
            dmr_df = pd.read_csv(fd, sep='\t')#, names=dmr_ss_header)
        except Exception as e:
            print(f"File {path2dmr} did not load => {e}")
            continue

    # Add column for genome browser coords and set as index
    dmr_df.insert(3, "location", dmr_df["chrom"] + ":" + dmr_df["start"].astype(str) + "-" + dmr_df["end"].astype(str))
    dmr_df.set_index(dmr_df['location'], inplace=True, drop=True)
    
    if bebic in dmr_all_chrs_dfs.keys():
        dmr_all_chrs_dfs[bebic] = pd.concat(
            [dmr_all_chrs_dfs[bebic], dmr_df])
    else:
        dmr_all_chrs_dfs[bebic] = dmr_df

# print(dmr_all_chrs_dfs[bebic].shape)
# display(dmr_all_chrs_dfs[bebic].head())

dmr_all_chrs_dfs[bebic] = filter_to_significance(dmr_all_chrs_dfs[bebic].copy(), pval = significance_threshold)
dmr_top20_dfs[bebic] = top_num(dmr_all_chrs_dfs[bebic], 20)
dmr_top50_dfs[bebic] = top_num(dmr_all_chrs_dfs[bebic], 50)
dmr_top100_dfs[bebic] = top_num(dmr_all_chrs_dfs[bebic], 100)

# Print summary data to console
df = dmr_all_chrs_dfs[bebic]
print(f'{bebic} number of significant positions (total)')
print(len(df))
print(f'{bebic} number of significant positions (positive)')
print(len(df[df['balanced_effect_size'] > 0]))
print(f'{bebic} number of significant positions (negative)')
print(len(df[df['balanced_effect_size'] < 0]))

# clean up dataframe to info I want
df = clean_data(df, cpg=True)
# Export to csv
df.to_csv(path2write_dmrdata, sep=',')
