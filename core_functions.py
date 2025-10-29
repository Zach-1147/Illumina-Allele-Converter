import pandas as pd
import numpy as np

# - - - - Dev params - - - #
manifest_size = 2000000
final_report_size = 150000
# - - - - - - - - - - - - -

#define a function to infer initial_convention from the file, prior to obtaining SNPs to flip.
def infer_initial(final_report_path, manifest_path):
    
    initial_conv = "FWD"
    
    return initial_conv

#define a function to process the manifest file, depending on the user-inputed initial_convention and target_convention, and return a list of SNPs that need to be flipped to achieve the target convention
def get_flip_snps(final_report_path, initial_convention=None, target_convention=None): 
    
    global flip_list

    #make a list from our two convention vars
    conv_vars = [initial_convention, target_convention]
    
    #flip one is going from FWD to PLUS or PLUS TO FWD
    flip_one = (("FWD" in conv_vars) and ("PLUS" in conv_vars))
    
    #generate combinations for the other possibilities
    flip_two = (("TOP" in conv_vars) and ("PLUS" in conv_vars))
    flip_three = (("TOP" in conv_vars) and ("FWD" in conv_vars))

    #proceed with only columns we need, based on the situation (differs by only one, but we want performance)
    if flip_one or flip_three:
        cols_needed =["Name", "IlmnStrand", "RefStrand", "SourceStrand"]
    
    #we dont need the sourcestrand for this combination
    elif flip_two:
        cols_needed = ["Name", "IlmnStrand", "RefStrand"]
        
    #read in chunks with pandas
    chunk_iter = pd.read_csv(final_report_path, sep=",", chunksize=round(manifest_size/10), dtype=str, skiprows=7, usecols=cols_needed)     
    
    flip_list = []
    rows_read = 0
    max_rows = manifest_size
         
    for chunk in chunk_iter:
        
        # **************************************************************************
        # Main Logic for identifying which SNPs need to be reverse complemented.
        # --------------------------------------------------------------------------
        '''
        We will now identify data to keep within these chunks conditionally, based on the initial_convention. The goal is to keep manifest data only for SNPs that need to be flipped in our final_report.txt. The keep condition will change depending on the initial_convention.
        '''
        
        if flip_one:
            #create boolean masks to identify which SNPs we need to flip, based on criteria specific to each flip situation.
            mask_1 = (chunk['RefStrand'] == '+') & (chunk['IlmnStrand'] != chunk['SourceStrand'])
            mask_2 = (chunk['RefStrand'] == '-') & (chunk['IlmnStrand'] == chunk['SourceStrand'])
            
            mask_conditions = (mask_1 | mask_2)
            
        elif flip_two:
            mask_1 = (chunk['RefStrand'] == '+') & (chunk['IlmnStrand'] == "BOT")
            mask_2 = (chunk['RefStrand'] == '-') & (chunk['IlmnStrand'] == "TOP")
            
            mask_conditions = (mask_1 | mask_2)
            
            
        # --------------------------------------------------------------------------
        # **************************************************************************
            

        #filter the chunk and retain only the Name column, since we just need the list of SNPs
        filtered = chunk.loc[mask_conditions, 'Name']

        #append SNP names to list
        flip_list.extend(filtered.tolist())
        
        rows_read += len(chunk)
        
        if max_rows and rows_read >= max_rows:
            break
        
    #print update
    print(f"\nAlleles for {len(flip_list)} loci will be reverse complimented to convert from [{initial_convention}] to [{target_convention}].\n")
    flip_list = set(name.strip().lower() for name in flip_list)
    
    return flip_list

#Function to flip alleles in provided Illumina genotype final report
def flip_alleles(filename, chunk_size=round(final_report_size/10), max_rows=final_report_size, sep='\t'):
    
    skiprows = 0

    data_chunks = []
    rows_processed = 0
    snps_flipped = 0
    
    #read in chunks with pandas
    chunk_iter = pd.read_csv(filename, sep=sep, chunksize=chunk_size, dtype=str, skiprows=skiprows)
    for chunk in chunk_iter:
        
        #for now, we will handle only this one.
        snp_col = "manifest_name"
        
        #create mask for only SNPs that need to be flipped using manifest_name
        mask = chunk[snp_col].apply(lambda x: x.strip().lower() in flip_list)
        
        #create a transition table for use with numpy's char translate function
        flip_table = str.maketrans('ATCGatcg', 'TAGCtagc')
        
        #obtain list of columns with genotype data (all but the snp col really)
        geno_cols = [c for c in chunk.columns if c != snp_col]    

        #convert subset to numpy array for efficiency
        arrai = chunk.loc[mask, geno_cols].to_numpy(dtype=str)

        #use numpy's char.translate (vectorized string ops)
        flipped = np.char.translate(arrai, flip_table)

        #assign back
        chunk.loc[mask, geno_cols] = flipped
            
        snps_flipped += mask.sum()

        #append all to to data_chunks
        data_chunks.append(chunk)
        
        rows_processed += chunk.shape[0]

        if max_rows and rows_processed >= max_rows:
            break
        
    #return df obeject, along with some variables for printing updates
    data = pd.concat(data_chunks, ignore_index=True)
    return data, snps_flipped, rows_processed