import time
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
from core_functions import get_flip_snps, flip_alleles, infer_initial

#write output to csv in wd
output_filename = "out.csv"

# - - - - -  Development vars - - - - - - #

final_report_path = "data/sim_genotype_file.txt"
manifest_path = "data/GDA-8v1-0_D1.csv"
write = False
initial_conv = None
target_conv = "PLUS"

# - - - - -  MAIN PROGRAM - - - - - - #

manif_str = time.time()

#First check if inital_convention is empty
if not(target_conv):
    print(f"\nNo Target convention specified\n")
    sys.exit()

if not(initial_conv):
    print(f"\nInferring initial convention from provided data...")
    
    #call function to infer genotype (filler function fr the time being)
    initial_conv = infer_initial(final_report_path, manifest_path)

#read in manifestfile
flip_list = get_flip_snps(manifest_path, initial_conv, target_conv)

manif_end = time.time()
manif_time = manif_end - manif_str

gen_str = time.time()

#call function to read in final_report.txt
genotype_df, snps_flipped, rows_processed = flip_alleles(final_report_path)

# - - - - -  PRINTING UPDATES - - - - - - #
genend = time.time()
gen_time = genend - gen_str

run_time = genend - manif_str

print("\nTotal run time: {:.2f} seconds".format(run_time))
print("-----------------------")
print("  --Manifest processing time: {:.2f} seconds".format(manif_time))
print("  --Genotype processing time: {:.2f} seconds".format(gen_time))
print(f"\nTotal rows processed: {rows_processed}")
print(f"Total SNPs flipped: {snps_flipped}")

if write:
    genotype_df.to_csv(output_filename, index=False)
    print(f"Data written to {output_filename}\n\n")