# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr

# Imports
import datetime
import os
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
Entrez.email = "RABVexample@RABVexample.com"     # Always tell NCBI who you are

# Functions
def read_and_process_accession_list(
    input_file_name, 
    output_fasta_file_name, 
    output_metadata_file_name, 
    length_threshold, 
    desired_segment
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    def country_extraction(location):
        """
        Function to extract country
        """
        # Upper case location and extract before colon
        location = location.upper().split(":")[0]
            
        return location

     # Open output files
    output_fasta_file = open(output_fasta_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")
    total_accessions_count = 0
    removed_accessions_count = 0

    # Write header for metadata file
    header = [
        "strain",
        "virus",
        "segment",
        "host",
        "accession",
        "date",
        # For now, just using a single location field
        # "region",
        # "country",
        # "division",
        # "city",
        "location",
        "country",
        "database",
        "authors",
        "url",
        "title",
        "journal",
        "paper_url",
    ]
    header = "\t".join(header) + "\n"
    output_metadata_file.write(header)
    
    # Open file with list of accession numbers
    with open(input_file_name, "r") as input_file:
        for line in input_file:

            # Initialize results for metadata
            strain = "MISSING"
            virus = "?"
            segment = "?"
            host = "?"
            accession = "?"
            date = "?"
            # For now, just using a single location field
            # region = ""
            # country = ""
            # division = ""
            # city = ""
            location = "?"
            country = "?"
            database = "?"
            authors = "?"
            url = "?"
            title = "?"
            journal = "?"
            paper_url = "?"
            features_flag = False
            sequence_flag = False
            nucleotide_sequence = ""
            products = None
            reference_count = 0
            length = 0
            backup_date = ""

            # Extract current accession ID
            accession_from_list = line.split()[0]
            total_accessions_count += 1

            # Check if accession is to be excluded
            if accession_from_list[:-2] in snakemake.params.accesstions_to_exclude:
                print(f"{accession_from_list} excluded based on config file!\n")
                removed_accessions_count += 1
                continue

            # Retrieve genbank file for accession ID
            entrez_genbank = Entrez.efetch(
                db="nucleotide", 
                id=accession_from_list, 
                rettype="genbank", 
                retmode="text"
                )
            print(f"Processing {accession_from_list}")
            # Parse genbank file line by line to retrieve all metadata
            for line in entrez_genbank:

                # Process line by removing spaces
                line = " ".join([ele for ele in line.split(" ") if ele != ""])
                split_line = line.split(" ")

                # Extract CDS products and determine which segment they come from
                if "/product" in line and desired_segment != None:
                    curr_product = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "glycoprotein" in curr_product:
                        products = "G"
                    else:
                        print(f"ERROR: {curr_product} product not known")



                # Extract feature information from genbank file
                if features_flag == True:
                    if "/isolate" in line or "/strain" in line:
                        strain = line.replace("\n", "").replace("\"", "").split("=")[1]
                        # Remove special characters ()[]{}|#>< not desired in nextstrain
                        strain = strain.replace("@", "-")
                        for char in ["(",")","[","]","{","}","|","#",">","<"]:
                            strain = strain.replace(char, "")
                    if "/organism" in line:
                        virus = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/host" in line:
                        host = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/segment" in line:
                        segment = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/country" in line:
                        # For now, just using a single location field
                        location = line.replace("\n", "").replace("\"", "").split("=")[1]
                        country = country_extraction(location)
                    if "/collection_date" in line:
                        unformatted_date = line.replace("\n", "").replace("\"", "").split("=")[1]
                        if len(unformatted_date.split("-")) == 1:
                            if "/" in unformatted_date:
                                date = datetime.datetime.strptime(unformatted_date.split("/")[0], "%Y").strftime("%Y") + "-XX-XX"
                            else:
                                date = datetime.datetime.strptime(unformatted_date, "%Y").strftime("%Y") + "-XX-XX"
                        elif len(unformatted_date.split("-")) == 2:
                            for fmt in ["%b-%Y", "%Y-%m"]:
                                test_res = True
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                except:
                                    print(f"ERROR! {unformatted_date} not recongized as {fmt}")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                    break
                        elif len(unformatted_date.split("-")) == 3:
                            for fmt in ["%d-%b-%Y", "%Y-%m-%d"]:
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                except:
                                    print(f"ERROR! {unformatted_date} not recongized as {fmt}")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                    break
                        else:
                            print("Datetime not between 1 and 3")

                # Extract nucleotide sequence
                if sequence_flag == True and split_line[0] != "//": 
                    nucleotide_sequence += "".join(split_line[1:])
                
                if split_line[0] == "LOCUS":
                    # Create a backup date based on top line of genbank
                    unformatted_backup_date = split_line[-1].replace("\n", "")
                    if len(unformatted_backup_date.split("-")) == 1:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%Y").strftime("%Y") + "-XX-XX"
                    elif len(unformatted_backup_date.split("-")) == 2:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%b-%Y").strftime("%Y-%m") + "-XX"
                    elif len(unformatted_backup_date.split("-")) == 3:
                        for fmt in ["%d-%b-%Y", "%Y-%m-%d"]:
                            try:
                                datetime.datetime.strptime(unformatted_backup_date, fmt).strftime("%Y-%m-%d")
                            except:
                                print(f"ERROR! {unformatted_backup_date} not recongized as {fmt}")
                                continue
                            else:
                                backup_date = datetime.datetime.strptime(unformatted_backup_date, fmt).strftime("%Y-%m-%d")
                                break
                    else:
                        print("Datetime not between 1 and 3")
                    # Get sequence length
                    length = int(split_line[2])
                    continue
                if split_line[0] == "ACCESSION":
                    accession = split_line[1].replace("\n", "")
                    genbank_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
                    url = genbank_base_url + accession
                    continue
                if split_line[0] == "DBLINK":
                    database = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "REFERENCE":
                    reference_count += 1
                    continue
                if split_line[0] == "AUTHORS" and reference_count == 1:
                    authors = split_line[1].split(",")[0] + " et al"
                    continue
                if split_line[0] == "TITLE" and reference_count == 1:
                    title = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "JOURNAL" and reference_count == 1:
                    journal = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "PUBMED" and reference_count == 1:
                    pubmed_base_url = "https://pubmed.ncbi.nlm.nih.gov/"
                    paper_url = pubmed_base_url + split_line[1].replace("\n", "")
                    continue
                if split_line[0] == "FEATURES":
                    features_flag = True # Set features flag to then extract feature information
                    continue 
                if split_line[0] == "ORIGIN":
                    sequence_flag = True # Set sequence flag to extract nucleotide sequence
                    continue

            # Check length of sequence and do not add if below threshold
            if length < length_threshold[0] or length > length_threshold[1]:
                print(f"{accession_from_list} excluded because {length} not within length limits!\n")
                removed_accessions_count += 1
                continue

            # Check if only specific segments should be kept
            if desired_segment != None and segment != desired_segment:
                print(f"{accession_from_list} excluded because {segment} is not the desired segment ({desired_segment})!\n")
                removed_accessions_count += 1
                continue

            # If not collection date is found, add date from top of genbank file
            if date == "?":
                date = backup_date

            # Check if sequence is below ambiguous base threshold
            if nucleotide_sequence.upper().count("N")/length > snakemake.params.max_frac_N:
                print(f"{accession_from_list} excluded because of high N count!\n")
                removed_accessions_count += 1
                continue
            
            # Join strain/isolate name with accession and date to make sure it is unique
            if strain != "MISSING":
                strain = strain + "_" + accession + "_" + date
            else:
                strain = accession + "_" + date

            # Replace slashes, periods, and spaces in name with underscores
            strain = strain.replace("/", "_")
            strain = strain.replace(". ", "-")
            strain = strain.replace(" ", "_")
            strain = strain.replace(".", "-")

            # Make sure every sequence has a fasta strain name
            assert strain != "MISSING", "Virus strain name is missing"

            # Check if segment is labeled and label if products are known
            if segment == "?" and products != None:
                segment = products

            # Create new metadata line
            new_metadata_line = "\t".join([
                strain,
                virus,
                segment,
                host,
                accession,
                date,
                # For now, just using a single location field
                # region,
                # country,
                # division,
                # city,
                location,
                country,
                database,
                authors,
                url,
                title,
                journal,
                paper_url,
            ])
            new_metadata_line += "\n"

            # Write new metadata line
            output_metadata_file.write(new_metadata_line)

            # Write current fasta sequence to output file
            output_fasta_file.write(f">{strain}\n")
            output_fasta_file.write(f"{nucleotide_sequence}\n")

    print(f"A total of {total_accessions_count} were processed and ")
    print(f"{total_accessions_count-removed_accessions_count} were retained!\n")   
    # Close files
    input_file.close()
    output_fasta_file.close()
    output_metadata_file.close()


def main():
    """
    Main method
    """

    # Input files
    list_of_accessions = str(snakemake.input.accession_list) 
    # Params
    length_threshold = (
        int(str(snakemake.params.genome_size_threshold_lower)), 
        int(str(snakemake.params.genome_size_threshold_upper))
    )
    # Empty String to None Conversion
    None_conversion = lambda i : i or None
    desired_segment = None_conversion(str(snakemake.params.desired_segment))

    # Output files
    fasta_output = str(snakemake.output.fasta_sequences) 
    metadata_output = str(snakemake.output.metadata) 

    read_and_process_accession_list(
        list_of_accessions, 
        fasta_output, 
        metadata_output, 
        length_threshold, 
        desired_segment
        )


if __name__ == "__main__":
    main()