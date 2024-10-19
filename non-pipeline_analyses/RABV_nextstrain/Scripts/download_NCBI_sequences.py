# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr

# Imports
import re
import datetime
import pandas as pd
from Bio import Entrez
Entrez.email = "RABVexample@RABVexample.com"     # Always tell NCBI who you are

# Functions
def read_and_process_accession_list(
    input_file_name, 
    output_fasta_file_name, 
    output_metadata_file_name, 
    length_threshold, 
    remove_duplicates,
    desired_gene,
    output_log_file,
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    def check_if_genbank_feature(feature):
        """
        Function to check if a feature is
        a genbank feature
        """

        if feature in snakemake.params.genbank_features:
            return True
        else:
            return False

    def country_extraction(location):
        """
        Function to extract country
        """
        # lower case location and extract before colon
        location = location.lower().split(":")[0]

        # Intialize results
        extracted_country = None
        region = None

        # Reformat specific countries
        if location == "viet nam":
            extracted_country = "Vietnam"
        elif location == "usa":
            extracted_country = "USA"
        elif location == "zaire":
            extracted_country = "Democratic Republic of the Congo"
        elif location == "democratic republic of the congo":
            extracted_country = "Democratic Republic of the Congo"
        elif location == "cote d'ivoire":
            extracted_country = "Cote d'Ivoire"
        elif location == "bosnia and herzegovina":
            extracted_country = "Bosnia and Herzegovina"
        elif location == "ussr":
            extracted_country = "Russia"
        else:
            # Capitalize first char in string
            formatted_list = []
            for substring in location.split(" "):
                formatted_list.append(substring.capitalize())
            
            extracted_country = " ".join(formatted_list)

        # Group by region
        if extracted_country in snakemake.params.asia:
            region = "Asia"
        elif extracted_country in snakemake.params.oceania:
            region = "Oceania"
        elif extracted_country in snakemake.params.africa:
            region = "Africa"
        elif extracted_country in snakemake.params.europe:
            region = "Europe"
        elif extracted_country in snakemake.params.south_america:
            region = "South America"
        elif extracted_country in snakemake.params.north_america:
            region = "North America"
        else:
            region = "?"
        
        # Return results
        return (region, extracted_country)
    
    
    def host_grouper(host):
        """
        Function to group host species
        """

        # Initialize dict and results
        host_dict = snakemake.params.host_grouper
        result = "UNKNOWN"

        # Iterate through host dict
        for grouped_host, host_list in host_dict.items():
            if host.lower() in host_list:
                result = grouped_host
                break

        # Return results
        return result
    

    def parse_CDS(raw_CDS_string, sequence):
        """
        Function to extract CDS
        """

        # Extract CDS coordinates
        temp = re.findall(r"\d+", raw_CDS_string)
        CDS_coordinates = list(map(int, temp))

        # Check to make sure there is an even number of coordinates
        assert len(CDS_coordinates) % 2 == 0, "Odd number of coordinates"

        # Intialize CDS string
        CDS_string = ""
        n = len(sequence)

        # Make reverse complement of sequence
        sequence_RC = sequence[::-1] 
        sequence_RC = sequence_RC.upper()
        sequence_RC = (
            sequence_RC
            .replace("A", "t")
            .replace("C", "g")
            .replace("T", "a")
            .replace("G", "c")
            .replace("N", "n")
        )

        # Process each CDS coordinates
        if "join" in raw_CDS_string and "complement" in raw_CDS_string:
            for i in range(len(CDS_coordinates)-1, -1, -2):
                CDS_string += sequence_RC[n - CDS_coordinates[i]:n - CDS_coordinates[i-1] + 1]
        elif "join" in raw_CDS_string:
            for i in range(0, len(CDS_coordinates)-1, 2):
                CDS_string += sequence[CDS_coordinates[i]-1:CDS_coordinates[i+1]]
        elif "complement" in raw_CDS_string:
            CDS_string += sequence_RC[n - CDS_coordinates[1]:n - CDS_coordinates[0] + 1]
        else:
            CDS_string += sequence[CDS_coordinates[0]-1:CDS_coordinates[1]]

        return CDS_string
    
    def virus_grouper(virus):
        """
        Function to group virus
        """
        if "Taiwan bat" in virus:
            return ("I", "Taiwan bat lyssavirus")
        elif "Gannoruwa bat" in virus:
            return ("I", "Gannoruwa bat lyssavirus")
        elif "Lleida bat" in virus:
            return ("III", "Lleida bat lyssavirus")
        elif "khujand" in virus:
            return ("I", "Khujand lyssavirus")
        elif "Lyssavirus Ozernoe" in virus:
            return ("I", "Ozernoe lyssavirus")
        elif "Lyssavirus shimoni" in virus:
            return ("II", "Shimoni lyssavirus")
        elif "Lyssavirus caucasicus" in virus:
            return ("III", "Caucasicus lyssavirus")
        elif "Lyssavirus bokeloh" in virus:
            return ("I", "Bokeloh lyssavirus")
        elif "Lyssavirus lagos" in virus:
            return ("II", "Lagos lyssavirus")
        elif "Lyssavirus aravan" in virus:
            return ("I", "Aravan lyssavirus")
        elif "Lyssavirus irkut" in virus:
            return ("I", "Irkut lyssavirus")
        elif "Lyssavirus duvenhage" in virus:
            return ("I", "Duvenhage lyssavirus")
        elif "Ikoma lyssavirus" in virus:
            return ("III", "Ikoma lyssavirus")
        elif "European bat lyssavirus 1" in virus:
            return ("I", "European bat lyssavirus 1")
        elif "European bat lyssavirus 2" in virus:
            return ("I", "European bat lyssavirus 2")
        elif "Lyssavirus mokola" in virus:
            return ("II", "Mokola lyssavirus")
        elif "Lyssavirus australis" in virus:
            return ("I", "Australian lyssavirus")
        elif "Lyssavirus rabies" in virus:
            return ("I", "Rabies lyssavirus")
        elif "Phala bat lyssavirus" in virus:
            return ("I", "Phala bat lyssavirus")
        elif "Matlo bat lyssavirus" in virus:
            return ("III", "Matlo bat lyssavirus")
        elif "Kotalahti bat lyssavirus" in virus:
            return ("I", "Kotalahti bat lyssavirus")
        elif "European bat lyssavirus" == virus:
            # Based on phylogeny these sequences seem to be EBL 1 
            return ("I", "European bat lyssavirus 1") 
        else:
            return virus



     # Open output files
    output_fasta_file = open(output_fasta_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")
    log_file = open(output_log_file, "w")
    total_accessions_count = 0
    removed_accessions_count = 0
    seen_sequences = []

    # Write header for metadata file
    header = [
        "strain",
        "virus",
        "phylogroup",
        "gene",
        "host",
        "accession",
        "date",
        "location",
        "region",
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
            phylogroup = "?",
            gene = "?"
            host = "?"
            accession = "?"
            date = "?"
            location = "?"
            region = "?"
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
            curr_CDS = False
            CDS = ""
            reference_count = 0
            length = 0
            backup_date = ""

            # Extract current accession ID
            accession_from_list = line.split()[0]
            total_accessions_count += 1

            # Check if accession is to be excluded
            if accession_from_list[:-2] in snakemake.params.accesstions_to_exclude:
                log_file.write(f"{accession_from_list} excluded based on config file!\n")
                removed_accessions_count += 1
                continue

            # Retrieve genbank file for accession ID
            entrez_genbank = Entrez.efetch(
                db="nucleotide", 
                id=accession_from_list, 
                rettype="genbank", 
                retmode="text"
                )
            log_file.write("\n")
            log_file.write(f"Processing {accession_from_list}\n")
            # Parse genbank file line by line to retrieve all metadata
            for line in entrez_genbank:

                # Process line by removing spaces
                line = " ".join([ele for ele in line.split(" ") if ele != ""])
                split_line = line.split(" ")

                # Check if not in CDS feature 
                if curr_CDS and check_if_genbank_feature(split_line[0]):
                    curr_CDS = False

                # Extract current CDS if segment not found and 
                # CDS does not contain a complement sequence b/c rabies 
                # should be all from 5' direction
                if split_line[0] == "CDS" and gene != desired_gene:
                    CDS = line.split(" ")[1]
                    curr_CDS = True

                # Extract CDS products and determine which segment they come from
                if desired_gene != None and ("/product" in line or "/gene" in line) and curr_CDS and gene != desired_gene:
                    gene = line.replace("\n", "").replace("\"", "").split("=")[1].upper()
                    def check_list(string_list, word):
                        """
                        Function to check if a word has any
                        matching strings in a list
                        """
                        for substring in string_list:
                            if substring == word:
                                return True
                        return False

                    if check_list(snakemake.params.gene_names, gene):
                        log_file.write(f"{gene} is {desired_gene}!\n")
                        gene = desired_gene
                    else:
                        log_file.write(f"ERROR: {gene} product not known or not {desired_gene}\n")

                # Extract feature information from genbank file
                if features_flag == True:
                    if "/isolate" in line or "/strain" in line:
                        strain = line.replace("\n", "").replace("\"", "").split("=")[1]
                        # Remove special characters ()[]{}|#>< not desired in nextstrain
                        strain = strain.replace("@", "-")
                        for char in ["(",")","[","]","{","}","|","#",">","<",";"]:
                            strain = strain.replace(char, "")
                    if "/organism" in line:
                        phylo_and_virus = virus_grouper(line.replace("\n", "").replace("\"", "").split("=")[1])
                        virus = phylo_and_virus[1]
                        phylogroup = phylo_and_virus[0]
                    if "/host" in line or "/lab_host" in line:
                        host = line.replace("\n", "").replace("\"", "").split("=")[1].split(";")[0]
                        grouped_host = host_grouper(host)
                        if grouped_host == "UNKNOWN":
                            log_file.write(f"{host} not unknown!\n")
                            host = "?"
                        else:
                            host = grouped_host
                    if "/geo_loc_name" in line:
                        # For now, just using a single location field
                        location = line.replace("\n", "").replace("\"", "").split("=")[1]
                        region_and_country = country_extraction(location)
                        region = region_and_country[0]
                        country = region_and_country[1]
                    if "/collection_date" in line:
                        unformatted_date = line.replace("\n", "").replace("\"", "").split("=")[1]
                        if len(unformatted_date.split("-")) == 1:
                            if "/" in unformatted_date:
                                date = datetime.datetime.strptime(unformatted_date.split("/")[0], "%Y").strftime("%Y") + "-XX-XX"
                            else:
                                date = datetime.datetime.strptime(unformatted_date, "%Y").strftime("%Y") + "-XX-XX"
                        elif len(unformatted_date.split("-")) == 2:
                            for fmt in ["%b-%Y", "%Y-%m"]:
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                except:
                                    log_file.write(f"ERROR! {unformatted_date} not recongized as {fmt}\n")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                    break
                        elif len(unformatted_date.split("-")) == 3:
                            for fmt in ["%d-%b-%Y", "%Y-%m-%d"]:
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                except:
                                    log_file.write(f"ERROR! {unformatted_date} not recongized as {fmt}\n")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                    break
                        else:
                            log_file.write("Datetime not between 1 and 3\n")

                # Extract nucleotide sequence
                if sequence_flag == True and split_line[0] != "//": 
                    nucleotide_sequence += "".join(split_line[1:]).replace("\n", "")
                
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
                                log_file.write(f"ERROR! {unformatted_backup_date} not recongized as {fmt}\n")
                                continue
                            else:
                                backup_date = datetime.datetime.strptime(unformatted_backup_date, fmt).strftime("%Y-%m-%d")
                                break
                    else:
                        log_file.write("Datetime not between 1 and 3\n")
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
                    authors = authors.replace("\n", "")
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

            # Check if sequence is within the desired phylogroups
            if phylogroup not in snakemake.params.phylogroups and accession_from_list[:-2] not in snakemake.params.accessions_to_include:
                log_file.write(f"{accession_from_list} excluded because virus is in phylogroup {phylogroup}!\n")
                removed_accessions_count += 1
                continue
            
            # Check length of sequence and do not add if below threshold
            if length < length_threshold[0] or length > length_threshold[1]:
                log_file.write(f"{accession_from_list} excluded because {length} not within length limits!\n")
                removed_accessions_count += 1
                continue

            # Check if only specific genes should be kept
            if desired_gene != None and gene != desired_gene:
                log_file.write(f"{accession_from_list} excluded because {gene} is not the desired gene ({desired_gene})!\n")
                removed_accessions_count += 1
                continue

            # Check if CDS was found else extract CDS from full sequence
            if gene == desired_gene and CDS == "":
                log_file.write(f"{accession_from_list} excluded because no CDS was found!\n")
                removed_accessions_count += 1
                continue
            else:
                # Extract CDS sequence
                if desired_gene != None:
                    log_file.write(f"CDS found: {CDS}\n")
                    nucleotide_sequence = parse_CDS(CDS, nucleotide_sequence)

                if remove_duplicates == "Yes":
                    # Check if nulceotide sequence is a duplicate
                    if nucleotide_sequence in seen_sequences and accession_from_list[:-2] not in snakemake.params.accessions_to_include:
                        log_file.write(f"{accession_from_list} excluded because duplicate sequence!\n")
                        removed_accessions_count += 1
                        continue
                    else:
                        seen_sequences.append(nucleotide_sequence)

            # If not collection date is found, add date from top of genbank file
            if date == "?":
                date = backup_date

            # Check if sequence is below ambiguous base threshold
            if nucleotide_sequence.upper().count("N")/len(nucleotide_sequence) > snakemake.params.max_frac_N:
                log_file.write(f"{accession_from_list} excluded because of high N count!\n")
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

            # Create new metadata line
            new_metadata_line = "\t".join([
                strain,
                virus,
                phylogroup,
                gene,
                host,
                accession,
                date,
                location,
                region,
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

    log_file.write("\n")
    log_file.write(f"A total of {total_accessions_count} were processed and ")
    log_file.write(f"{total_accessions_count-removed_accessions_count} were retained!\n")   
    # Close files
    input_file.close()
    output_fasta_file.close()
    output_metadata_file.close()
    log_file.close()


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
    desired_gene = None_conversion(str(snakemake.params.desired_gene))
    remove_duplicates = str(snakemake.params.remove_duplicates)

    # Output files
    fasta_output = str(snakemake.output.fasta_sequences) 
    metadata_output = str(snakemake.output.metadata)
    output_log_file = str(snakemake.log)

    read_and_process_accession_list(
        list_of_accessions, 
        fasta_output, 
        metadata_output, 
        length_threshold, 
        remove_duplicates,
        desired_gene,
        output_log_file,
        )


if __name__ == "__main__":
    main()