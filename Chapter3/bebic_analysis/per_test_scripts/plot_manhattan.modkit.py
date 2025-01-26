import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors

# sns.set_theme(style='white')
sns.set_theme()

################# Pick variables #################
name_of_test = "test11"
bebic = "BEBIC_B"
# bebic = "BEBIC_C"
exp_list = ["B3", "B4", "B5", "B6"]
# exp_list = ["C3", "C5", "C6"]

# modBase = "C"
modBase = "A"
motif = "A"
both_bases = False
# both_bases = True

# source is modkit or minfi?
data_source = "modkit"
# data_source = "minfi"
# score_threshold = 50

###################################################

################# Pipeline paths #################
wd = os.getcwd()
print(f"working directory: {wd}")
ad = "/external/analyses/lucy/bebic_pipeline"
print(f"analyses directory: {ad}")
output_dir = os.path.join(ad, name_of_test, "intermediate_data_files")
os.makedirs(output_dir, exist_ok=True)
plotting_dir = os.path.join(wd, "plots")
os.makedirs(plotting_dir, exist_ok=True)

sigDMP_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.{data_source}.csv")
sigDMPnames_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.{data_source}.named.csv")

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


################# Functions #################

# input dataframe:
# 	chrom	start	end	score	base	a_methylated	a_hydroxymethylated	a_total	b_methylated	b_hydroxymethylated	b_total	a_pct_modified	b_pct_modified	map_pvalue	effect_size	balanced_map_pvalue	balanced_effect_size	num_a_samples	num_b_samples	name
# 0	chr1	121695871	121695872	4.099420	C	15	1	141	34	0	146	0.113475	0.232877	0.327975	-0.119402	0.015769	-0.151005	3	4	1.1
# 1	chr1	124153984	124153985	4.901008	A	3	0	88	20	0	119	0.034091	0.168067	0.002442	-0.133976	0.005758	-0.119363	3	3	1.2
# 2	chr1	124177350	124177351	2.108251	C	42	2	140	28	1	149	0.314286	0.194631	0.033023	0.119655	0.034888	0.109184	4	3	1.3
# 3	chr1	124203828	124203829	2.788887	C	25	0	102	17	0	146	0.245098	0.116438	0.084494	0.128660	0.001078	0.177288	3	3	1.4
# 4	chr1	124204400	124204401	2.835215	C	49	1	105	104	1	165	0.476190	0.636364	0.041863	-0.160173	0.006696	-0.191870	3	4	1.5

# Important columns
# chrom	start	score	base    balanced_effect_size    name

def add_names_col(df):
    # get list for names
    df['name'] = None
    chr_count = 1
    chrom = 'chr1'
    for i, row in df.iterrows():
        # print(i)
        # print(row)
        if row['chrom'] == chrom:
            df.at[i, 'name'] = f'{str(row["chrom"])[3:]}.{str(chr_count)}'
            chr_count += 1
        else:
            chrom = row['chrom']
            chr_count = 1
            df.at[i, 'name'] = f'{str(row["chrom"])[3:]}.{str(chr_count)}'
            chr_count += 1
    return df


# def adjust_x_for_chromosomes(df, spacer_value=-1):
#     df = df.copy()
#     new_rows = []
#     count = 1

#     chromosomes = df['chrom'].unique()
#     num_to_add = 0
#     spacer_rows = []

#     for chrom in chromosomes:
#         chrom_df = df[df['chrom'] == chrom].copy()
#         chrom_df['x_pos'] = range(count, count + len(chrom_df))
#         count += len(chrom_df)

#         # Add the chrom_df to the new_rows list
#         new_rows.append(chrom_df)

#         # Append a spacer bar with NaN effect size and transparent color
#         spacer_row = {
#             'chrom': chrom,
#             'start': spacer_value,
#             'effect_size': 0,
#             'balanced_effect_size': 0,
#             'score': 0,
#             'x_pos': count,
#             'name': ""
#         }

#         spacer_df = pd.DataFrame([spacer_row])
#         new_rows.append(spacer_df)
#         spacer_rows.append(count - 1)

#         count += 1  # Increment the count for the next chromosome
#     new_rows.pop()
#     # Concatenate the list of DataFrames into a single DataFrame
#     df_adjusted = pd.concat(new_rows, ignore_index=True)

#     return df_adjusted, chromosomes, spacer_rows


# def run_plot(base, plot_data):
#     if data_source == 'modkit':
#         # switch sides so that increased methylation in Test is up, and vice versa
#         plot_data["balanced_effect_size"] = plot_data["balanced_effect_size"] * -1
    
#     if score_threshold:
#         plot_data = plot_data.loc[plot_data["score"] >= score_threshold]
    
#     # Apply the function to adjust x values
#     # plot_data_adjusted, chromosomes = adjust_x_for_chromosomes(plot_data)
#     plot_data_adjusted, chromosomes, spacer_rows = adjust_x_for_chromosomes(plot_data)
#     # display(plot_data_adjusted.head())

#     # Create a single plot with adjusted x values
#     fig = plt.figure(figsize=(24, 6))

#     # Define a color map based on the "score" column
#     norm = mcolors.Normalize(
#         vmin=0, vmax=plot_data_adjusted['score'].max())
#         # vmin=plot_data_adjusted['score'].min(), vmax=plot_data_adjusted['score'].max())
#     cmap = plt.get_cmap('RdPu')
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     plot_data_adjusted['color'] = plot_data_adjusted['score'].apply(
#         lambda x: sm.to_rgba(x))

#     # Mark spacer bars with a gray color
#     plot_data_adjusted['color'] = np.where(
#         plot_data_adjusted['effect_size'] == 0, 'gray', plot_data_adjusted['color'])

#     # Plot the data with graduated shading
#     bar_plot = sns.barplot(
#         data=plot_data_adjusted, x='x_pos', y='balanced_effect_size',
#         hue='x_pos', legend=False, palette=plot_data_adjusted['color'].tolist(), errorbar=None)

#     # Add labels to bars
#     # for i, bar in enumerate(bar_plot.containers):
#     #     if (i in spacer_rows):
#     #         continue
#     #     names = list(plot_data_adjusted['name'])
#     #     # name = f"S{names[i]} - {plot_data_adjusted['base'][i]}"
#     #     name = f"S{names[i]}"
#     #     bar_plot.bar_label(bar, [name], fontsize=9, padding=5, rotation=90)

#     # Customize plot
#     # bar_plot.set_ylim(0.0,0.27)
#     bar_plot.set_ylim(-0.4, 0.4)
#     if bebic == "BEBIC_B":
#         plt.title(
#             f'SW780: Effect Size by DMP per Chromosome: Motif {motif}, score > {score_threshold}',
#             fontsize=20, pad=22)
#     else:
#         plt.title(
#             f'RT4: Effect Size by DMP per Chromosome: Motif {motif}, score > {score_threshold}',
#             fontsize=20, pad=22)
#     plt.xlabel('Chromosome', fontsize=14, labelpad=10)
#     plt.ylabel('Effect Size', fontsize=14, labelpad=10)

#     # Add vertical lines at the end of each chromosome and store positions for tick labels
#     prev_val = 0
#     chrom_tick_positions = list()
#     for chrom in chromosomes:
#         max_pos = plot_data_adjusted[plot_data_adjusted['chrom'] == chrom]['x_pos'].max(
#         )
#         # Add vertical lines extending beyond the plot limits
#         val = max_pos - 1
#         if (chrom != chromosomes[-1]):
#             plt.gca().axvline(
#                 x=val, color='gray', linestyle='--',
#                 linewidth=1, ymin=-0.05, ymax=1, clip_on=False)
#         chrom_tick_positions.append((prev_val + ((val - prev_val) / 2)))
#         prev_val = val
#     chrom_tick_positions[-1] += 0.75

#     # Set x-axis tick labels with chromosome names
#     chrom_tick_labels = plot_data_adjusted['chrom'].unique()

#     plt.xticks(
#         ticks=chrom_tick_positions,
#         labels=chrom_tick_labels, rotation=90)

#     # Add color bar for score
#     cbar = plt.colorbar(sm, aspect=20, ax=plt.gca(), pad=0.02)
#     cbar.set_label('Score')

#     # Remove all borders
#     bar_plot.margins(0.01)
#     sns.despine(right=True)
#     plt.tight_layout()
#     plt.savefig(os.path.join(plotting_dir, f"effect_size_DMPs.{bebic}.{motif}.{data_source}.png"), dpi=300)

#     return


def plot_manhattan(dmp_all_chrs_df):
    plt.cla()
    sns.set_style('white')

    # for exp_name, dmr_effect_df in dmr_all_chrs_dfs.items():
    # print(f'---- {exp_name} ----')
    # display(dmr_effect_df.head())

    # -log_10(pvalue)
    # df['minuslog10pvalue'] = -np.log10(df.pvalue)
    dmp_all_chrs_df.chrom = dmp_all_chrs_df.chrom.astype('category')
    dmp_all_chrs_df.chrom = dmp_all_chrs_df.chrom.cat.set_categories(chr_names, ordered=True)
    dmr_effect_df = dmp_all_chrs_df.sort_values('chrom')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    dmr_effect_df['ind'] = range(len(dmr_effect_df))
    df_grouped = dmr_effect_df.groupby(['chrom'], observed=True)

    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)
    colors = ['#8960b3','#56ae6c']#, 'yellow']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        # group.plot(kind='scatter', x='ind', y='score',color=colors[num % len(colors)], ax=ax, s=2)
        group.plot(kind='scatter', x='ind', y='score',color=colors[num % len(colors)], ax=ax, s=4)
        if len(group) == 0:
            print(f'{name} has no data')
            continue
        x_labels.append(name[0])
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, rotation=45)
    ax.set_xlim(0, len(dmr_effect_df))
    # ax.set_ylim(0, 100)
    ax.set_xlabel('Chromosome')
    ax.grid(axis='y')
    sns.despine(right=True, top=True)
    ax.spines['left'].set_color('#6b6b6bff')
    ax.spines['bottom'].set_color('#6b6b6bff')
    
    if data_source == "minfi":
        # ax.set_yscale('log')
        p_text = "0.01"
    else:
        p_text = "0.05"
    
    if bebic == "BEBIC_B":
        cell_line = "SW780"
    else:
        cell_line = "RT4"
    plt.title(f"Significant {motif} DMPs on Cell Line {cell_line} (p < {p_text})", size=14)
    plt.savefig(os.path.join(plotting_dir, f"manhattan_DMPs.{bebic}.{motif}.{data_source}.png"), dpi=300) 

    return


# Load plot data
# with open("plotting/mod_data_for_genome_plot_names.BEBIC_B.txt", 'r') as fp:
#     plot_data_B = pd.read_csv(fp, sep='\t', index_col=0)
with open(sigDMP_fp, 'r') as fp:
    if data_source == 'modkit':
        plot_data = pd.read_csv(fp, sep=',', index_col=0)
    else:
        plot_data = pd.read_csv(fp, sep=',', index_col=None)
        

if data_source == "minfi":
    plot_data = plot_data.sort_values(by=['chrom', 'start'])
    
if "names" not in plot_data.columns:
    plot_data = add_names_col(plot_data)

plot_manhattan(plot_data)

plot_data.to_csv(sigDMPnames_fp, sep=',', index=False)

