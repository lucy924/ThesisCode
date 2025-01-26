import os
import numpy as np
import pandas as pd

# Options
name_of_test = "test11"
bebic = "BEBIC_B"

################# Pipeline paths #################
wd = os.getcwd()
print(f"working directory: {wd}")
ad = "/external/analyses/lucy/bebic_pipeline"
print(f"analyses directory: {ad}")
output_dir = os.path.join(ad, name_of_test, "intermediate_data_files")
plotting_dir = os.path.join(wd, "plots")

sigDMP_modkit_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.modkit.named.csv")
sigDMP_minfi_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.minfi.named.csv")
sigDMP_modkit_new_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.modkit.noted_minfi.csv")
sigDMP_minfi_new_fp = os.path.join(output_dir, f"sig_DMPs.{bebic}.minfi.noted_modkit.csv")

###################################################

################# Load in data #################
with open(sigDMP_modkit_fp, 'r') as fp:
    sigDMP_modkit = pd.read_csv(fp, sep=',', index_col=None)
with open(sigDMP_minfi_fp, 'r') as fp:
    sigDMP_minfi = pd.read_csv(fp, sep=',', index_col=None)

################# Do comparison #################
# Add column also_in_other_BEBIC: "yes" or "no"
# Create a set of tuples for quick lookup
set_B = set(zip(sigDMP_modkit['chrom'], sigDMP_modkit['start'], sigDMP_modkit['end']))

# Add the new column to sigDMP_C
sigDMP_minfi['also_in_other'] = sigDMP_minfi.apply(
    lambda row: 'yes' if (row['chrom'], row['start']) in set_B else 'no', axis=1
)

# Create a set of tuples for sigDMP_C
set_C = set(zip(sigDMP_minfi['chrom'], sigDMP_minfi['start']))

# Add the new column to sigDMP_B
sigDMP_modkit['also_in_other'] = sigDMP_modkit.apply(
    lambda row: 'yes' if (row['chrom'], row['start'], row['end']) in set_C else 'no', axis=1
)

# write file
sigDMP_modkit.to_csv(sigDMP_modkit_new_fp)
sigDMP_minfi.to_csv(sigDMP_minfi_new_fp)

print('files here: ')
print(sigDMP_modkit_new_fp)
print(sigDMP_minfi_new_fp)
print('...done')
