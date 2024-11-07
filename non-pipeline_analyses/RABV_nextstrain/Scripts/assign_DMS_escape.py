# Description:
# Python script that predicts DMS
# escape for natural sequences

# Author:
# Caleb Carr

# Imports
import pandas as pd
from Bio import AlignIO 
from augur.utils import write_json
from collections import defaultdict


# Functions
def process_alignment(alignment, metadata, reference_strain, filter_params=None):
    """
    Process natural sequences and return a dataframe with
    all mutations relative to the reference strain
    """

    # Load alignment and metadata info
    natural_seqs_df = pd.DataFrame(columns=["strain", "sequence"])
    metadata_df = pd.read_csv(metadata, sep="\t")

    # Add alignment sequence to dataframe
    for curr_fasta in AlignIO.read(alignment, "fasta"):
        natural_seqs_df.loc[len(natural_seqs_df.index)] = [
            str(curr_fasta.id),
            str(curr_fasta.seq),
        ] 
    
    # Merge sequences and metadata
    natural_seqs_df = (
        natural_seqs_df.merge(
            metadata_df,
            how="left",
            on=["strain"],
            validate="one_to_one",
        )
    )

    # Filter for only desired sequences
    if filter_params != None:
        for param_key,param_value in filter_params.items():
            natural_seqs_df = (
                natural_seqs_df
                .query(f"`{param_key}` == @param_value")
                .reset_index(drop=True)
            )
    

    def get_muts(seq1, seq2):
        """
        Function that returns a list of all muts
        that are different between two sequences.
        """
        list_of_muts = []
        site = 1
        for s1, s2 in zip(seq1, seq2):
            if s1 != s2 and s1 != "-" and s2 != "-":
                list_of_muts.append(s1 + str(site) + s2)
            site += 1
        return list_of_muts
    
    # Get Reference sequence for comparison
    ref_sequence = (
        natural_seqs_df.loc[natural_seqs_df["strain"] == reference_strain].reset_index(drop=True).at[0,"sequence"]
    )
    
    # Get all mutations relative to reference sequence
    natural_seqs_df["natural_muts"] = (
        natural_seqs_df["sequence"].apply(lambda x: get_muts(ref_sequence, x))
    )

    return natural_seqs_df

def calculate_antibody_escape_and_write_json(
        escape_data, 
        natural_seqs_df, 
        site_column, 
        escape_column, 
        antibody_name, 
        output_json,
        output_log_file,
        site_map_csv
    ):
    """
    Function that calculates predicted escape for each natural sequence
    and then outputs to a json file.
    """

    # Load escape data
    escape_df = pd.read_csv(escape_data)

    # Load site map and merge with escape_df
    site_map = (
        pd.read_csv(site_map_csv)
        .rename(columns={
            "reference_site" : "site",
            "reference_wt" : "wildtype",
        })
        .drop(columns=["sequential_site", "sequential_wt"])
    )
    escape_df = (
        escape_df.merge(
            site_map,
            how="left",
            on=["site", "wildtype"],
            validate="many_to_one",
        )
    )

    # Create natural sequence mutation col
    escape_df["natural_sequence_mutation"] = (
        escape_df["wildtype"] + escape_df["natural_sequence_site"].astype(str) + escape_df["mutant"]
    )

    # Floor escape and extract single antibody data
    escape_df["floored_escape"] = escape_df["escape"].clip(lower=0)
    escape_df = escape_df.query("antibody == @antibody_name")

    # Create escape mutation dictonary
    escape_dict = dict(zip(escape_df[site_column], escape_df[escape_column]))

    def get_escape(strain, muts_list, escape_dict, output_log_file):
        """
        Function to calculate predicted DMS escape for each strain
        """

        # Initialize log file
        log_file = open(output_log_file, "a")

        total_escape = 0
        if len(muts_list) == 0:
            log_file.write(f"WARNING: {strain} does not have any mutations!\n")
            log_file.write("\n")
            log_file.close()

            return total_escape
        else:
            for mut in muts_list:
                escape = escape_dict.get(mut, -1000)
                if escape == -1000:
                    log_file.write(f"WARNING: {mut} in {strain} does not have any DMS data!\n")
                else:
                    total_escape += escape

            log_file.write(f"{strain} has a total escape of {total_escape}\n")
            log_file.write("\n")
            log_file.close()
            return total_escape

    
    # Get all mutations relative to reference sequence
    predicted_DMS_escape_col = f"predicted {antibody_name} escape"
    natural_seqs_df[predicted_DMS_escape_col] = (
        natural_seqs_df.apply(lambda x: get_escape(x["strain"], x["natural_muts"], escape_dict, output_log_file), axis=1)
    )

    # Create return dictionary to output as json
    ret_json = {
        "nodes": defaultdict(dict)
    }
    for idx, row in natural_seqs_df.iterrows():
        ret_json["nodes"][row["strain"]][predicted_DMS_escape_col] = row[predicted_DMS_escape_col]
    
    # Write final json
    write_json(ret_json, output_json)




def main():
    """
    Main method
    """

# Input files
escape_data = str(snakemake.input.escape_data)
alignment = str(snakemake.input.alignment)
metadata = str(snakemake.input.metadata)

# Params
antibody_name = str(snakemake.params.antibody_name)
site_map = str(snakemake.params.site_map)
filter_params = snakemake.params.filter_params
reference_strain = str(snakemake.params.reference_strain)
site_column = str(snakemake.params.site_column)
escape_column = str(snakemake.params.escape_column)

# Output files
output_log_file = str(snakemake.log)
output_json = str(snakemake.output.output_json)

# Process natural sequences and find mutations compared to reference
natural_seqs_df = process_alignment(alignment, metadata, reference_strain, filter_params)
# Calculate predicted escape
calculate_antibody_escape_and_write_json(
    escape_data, 
    natural_seqs_df, 
    site_column, 
    escape_column, 
    antibody_name, 
    output_json,
    output_log_file,
    site_map,
)

if __name__ == "__main__":
    main()