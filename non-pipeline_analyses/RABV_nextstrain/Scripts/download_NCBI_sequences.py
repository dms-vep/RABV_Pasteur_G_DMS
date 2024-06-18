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
    desired_segment,
    remove_duplicates,
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    def country_extraction(location):
        """
        Function to extract country
        """
        # lower case location and extract before colon
        location = location.lower().split(":")[0]

        extracted_country = None

        # Reformat specific countries
        if location == "viet nam":
            extracted_country = "Vietnam"
        elif location == "usa":
            extracted_country = "USA"
        elif location == "zaire":
            extracted_country =  "Democratic Republic of the Congo"
        else:
            # Capitalize first char in string
            formatted_list = []
            for substring in location.split(" "):
                formatted_list.append(substring.capitalize())
            
            extracted_country = " ".join(formatted_list)

        asia = [
            "Israel",
            "Palestine",
            "Lebanon",
            "Jordan",
            "Syria",
            "Georgia",
            "Iraq",
            "Armenia",
            "Saudi Arabia",
            "Kuwait",
            "Yemen",
            "Bahrain",
            "Qatar",
            "Iran",
            "United Arab Emirates",
            "Oman",
            "Uzbekistan",
            "Afghanistan",
            "Kazakhstan",
            "Pakistan",
            "Tajikistan",
            "Maldives",
            "Kyrgyzstan",
            "India",
            "Sri Lanka",
            "Nepal",
            "Asia",
            "Bangladesh",
            "Myanmar",
            "Thailand",
            "Malaysia",
            "Mongolia",
            "Laos",
            "Singapore",
            "China",
            "Cambodia",
            "Vietnam",
            "Hong Kong",
            "Indonesia",
            "Brunei",
            "Taiwan",
            "Philippines",
            "Timor-Leste",
            "South Korea",
            "Japan",
        ]

        oceania = [
            "French Polynesia",
            "Palau",
            "Australia",
            "Papua New Guinea",
            "Oceania",
            "Micronesia",
            "Solomon Islands",
            "Vanuatu",
            "Marshall Islands",
            "New Zealand",
            "Kiribati",
            "Fiji",
        ]

        africa = [
            "Cabo Verde",
            "Gambia",
            "Guinea-Bissau",
            "Senegal",
            "Sierra Leone",
            "Guinea",
            "Liberia",
            "Mauritania",
            "CÃ´te d'Ivoire",
            "Morocco",
            "Mali",
            "Burkina Faso",
            "Ghana",
            "Togo",
            "Benin",
            "Algeria",
            "Sao Tome and Principe",
            "Nigeria",
            "Tunisia",
            "Niger",
            "Equatorial Guinea",
            "Gabon",
            "Cameroon",
            "Republic of the Congo",
            "Libya",
            "Namibia",
            "Angola",
            "Chad",
            "Central African Republic",
            "Democratic Republic of the Congo",
            "Africa",
            "Botswana",
            "South Africa",
            "Zambia",
            "Lesotho",
            "Sudan",
            "South Sudan",
            "Zimbabwe",
            "Burundi",
            "Rwanda",
            "Egypt",
            "Eswatini",
            "Uganda",
            "Malawi",
            "Mozambique",
            "Tanzania",
            "Kenya",
            "Ethiopia",
            "Djibouti",
            "Union of the Comoros",
            "Madagascar",
            "Somalia",
            "Seychelles",
            "Mauritius",
        ]

        europe = [
            "Iceland",
            "Portugal",
            "Ireland",
            "Spain",
            "United Kingdom",
            "Andorra",
            "France",
            "Belgium",
            "Netherlands",
            "Luxembourg",
            "Monaco",
            "Switzerland",
            "Norway",
            "Denmark",
            "Liechtenstein",
            "Germany",
            "Europe",
            "Italy",
            "Malta",
            "Slovenia",
            "Sweden",
            "Austria",
            "Croatia",
            "Czech Republic",
            "Bosnia and Herzegovina",
            "Poland",
            "Slovakia",
            "Montenegro",
            "Albania",
            "Hungary",
            "Serbia",
            "Kosovo",
            "North Macedonia",
            "Greece",
            "Lithuania",
            "Romania",
            "Bulgaria",
            "Latvia",
            "Estonia",
            "Finland",
            "Belarus",
            "Moldova",
            "Ukraine",
            "Cyprus",
            "Turkey",
            "Azerbaijan",
            "Russia",
        ]

        south_america = [
            "Argentina",
            "Uruguay",
            "Chile",
            "Paraguay",
            "Bolivia",
            "South America",
            "Brazil",
            "Peru",
            "Ecuador",
            "Colombia",
            "Suriname",
            "Guyana",
            "Venezuela",
            "Trinidad and Tobago",
            "Bonaire",
            "Curacao",
            "Aruba",
        ]

        north_america = [
            "Panama",
            "Costa Rica",
            "Grenada",
            "Nicaragua",
            "Saint Vincent and the Grenadines",
            "Barbados",
            "El Salvador",
            "Saint Lucia",
            "Honduras",
            "Dominica",
            "Guatemala",
            "Guadeloupe",
            "Belize",
            "Antigua and Barbuda",
            "Saint Kitts and Nevis",
            "Saint BarthÃ©lemy",
            "Sint Maarten",
            "Saint Martin",
            "Jamaica",
            "Dominican Republic",
            "Haiti",
            "Mexico",
            "Cuba",
            "Bahamas",
            "North America",
            "Bermuda",
            "USA",
            "Canada",
        ]

        region = None
        if extracted_country in asia:
            region = "Asia"
        elif extracted_country in oceania:
            region = "Oceania"
        elif extracted_country in africa:
            region = "Africa"
        elif extracted_country in europe:
            region = "Europe"
        elif extracted_country in south_america:
            region = "South America"
        elif extracted_country in north_america:
            region = "North America"
        else:
            region = "?"
        
        return (region, extracted_country)
    
    def check_if_genbank_feature(feature):
        """
        Function to check if a feature is
        a genbank feature
        """

        genbank_features = [
            "assembly_gap",
            "C_region",
            # "CDS",
            "centromere",
            "D-loop",
            "D_segment",
            "exon",
            "gap",
            "gene",
            "iDNA",
            "intron"
            "J_segment",
            "mat_peptide",
            "misc_binding",
            "misc_difference",
            "misc_feature",
            "misc_recomb",
            "misc_RNA",
            "misc_structure",
            "mobile_element",
            "modified_base",
            "mRNA",
            "ncRNA",
            "N_region",
            "old_sequence",
            "operon",
            "oriT",
            "polyA_site",
            "precursor_RNA",
            "prim_transcript",
            "primer_bind",
            "propeptide",
            "protein_bind",
            "regulatory",
            "repeat_region",
            "rep_origin",
            "rRNA",
            "S_region",
            "sig_peptide",
            "source",
            "stem_loop",
            "STS",
            "telomere",
            "tmRNA",
            "transit_peptide",
            "tRNA",
            "unsure",
            "V_region",
            "V_segment",
            "variation",
            "3'UTR",
            "5'UTR",
        ]

        if feature in genbank_features:
            return True
        else:
            return False

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


    def host_grouper(host):
        """
        Function to group host
        """

        human_list = ["Homo", "rabies patient"]

        canidae_list = [
            "dog", 
            "Dog",
            "Nyctereutes", 
            "fox", 
            "Otocyon", 
            "Pseudalopex", 
            "wolf", 
            "Dusicyon", 
            "Canis", 
            "canine", 
            "jackel", 
            "jaskal", 
            "jackal",
            "Jackal", 
            "coyote", 
            "Cerdocyon", 
            "vulpes", 
            "Urocyon",
        ]

        bat_list = [
            "Myotis", 
            "myotis",
            "Nycticeinops", 
            "Pipistrellus", 
            "Pteropus", 
            "Miniopterus", 
            "Hipposideros", 
            "Eidolon", 
            "Eptesicus", 
            "Artibeus",
            "Tadarida", 
            "Desmodus", 
            "Epomophorus",
            "Stenodermatinae", 
            "Nyctinomops", 
            "Lasiurus",
            "Nycticeius", 
            "Lasionycteris", 
            "Parastrellus",
            "bat",
            "Bat",
            "Corynorhinus", 
            "Rousettus", 
            "Histiotus",
            "Antrozous", 
            "Chiroptera", 
            "Molossus",
            "Eumops", 
            "Saccolaimus",
            "Nyctalus", 
            "Aeorestes",
            "Carollia",
            "Dasypterus",
        ]

        feline_list = [
            "feline",
            "Felis",
            # "cat" == host:
            "cat;",
            "Lynx",
            "domestic cat",
            "wild cat",
            "tiger",
            "Panthera",
            "bobcat",
            "bob cat"
            "puma",
            "Felidae"
        ]

        bovidae_list = [
            "bovine", 
            "goat", 
            "cattle",
            "calf",
            "Calf",
            "caprine",
            "Tragelaphus",
            "yak",
            "sheep",
            "taurus",
            "cow",
            "Bubalus",
            "Ovis",
            "buffalo",
            "Buffalo",
            "Bubaline",
            "Capra",
            "nilgai",
            "ovine",
            "bull"
        ]

        musteloidea_list = [
            "Raccoon",
            "raccon",
            "raccoon",
            "skunk",
            "marten",
            "Melogale",
            "Meles",
            "Procyon",
            "Mephitis",
            "kinkajou",
            "Bassariscus",
            "Nasua",
            "Lontra",
            "Mephitidae",
            "Martes",
            "badger",
            "Mellivora",
            "Potos",
            "Ferret",
        ]

        rodent_list = [
            "rodent",
            "beaver",
            "woodchuck",
            "mice",
            "Mus", 
            "mouse",
            "Marmota",
            "squirrel",
        ]

        feliformia_list = [
            "hyena",
            "Hyena",
            "mongoose",
            "mongose",
            "Galerella",
            "Civettictis",
            "Civet",
            "civet",
            "Paradoxurus",
            "Genetta",
            "Hyaena",
            "Proteles",
            "Suricata",
            "Cynictis",
            "Atilax",
        ]

        equidae_list = [
            "Equus", 
            "horse", 
            "mule", 
            "donkey", 
            "zebra",
            "equino", 
            "equine"
        ]

        cervidae_list = [
            "deer",
            "Deer",
            "Odocoileus",
            "Axis",
            "doe",
            "elk",
        ]

        camelidae_list = [
            "camel",
            "dromedary"
        ]

        ursoidea_list = [
            "bear",
            "Melursus",
            "Ursus"
        ]

        callitrichidae_list = [
            "Callithrix"
        ]

        suidae_list = [
            "Sus"
        ]

        tayassuidae_list = [
            "Pecari"
        ]

        soricidae_list = [
            "Crocidura"
        ]

        ameridelphia_list = [
            "Didelphis"
        ]

        lab_based_list = [
            "cell culture", 
            "vaccine strain",
        ]

        def check_list(substring_list, word):
            """
            Function to check if a word has any
            matching substrings in a list
            """
            for substring in substring_list:
                if substring in word:
                    return True
            return False

        if check_list(human_list, host):
            return "Human"
        elif check_list(canidae_list, host):
            return "Canine"
        elif check_list(bat_list, host):
            return "Bat"
        elif check_list(feline_list, host) or "cat" == host:
            return "Feline"
        elif check_list(bovidae_list, host):
            return "Bovine"
        elif check_list(feliformia_list, host):
            return "Feliformia"
        elif check_list(musteloidea_list, host):
            return "Musteloidea"
        elif check_list(rodent_list, host):
            return "Rodent"
        elif check_list(equidae_list, host):
            return "Equine"
        elif check_list(cervidae_list, host):
            return "Cervine"
        elif check_list(camelidae_list, host):
            return "Camelidae"
        elif check_list(ursoidea_list, host):
            return "Ursoidea"
        elif check_list(callitrichidae_list, host):
            return "Callitrichidae" 
        elif check_list(suidae_list, host):
            return "Suidae" 
        elif check_list(tayassuidae_list, host):
            return "Tayassuidae" 
        elif check_list(soricidae_list, host):
            return "Soricidae" 
        elif check_list(ameridelphia_list, host):
            return "Ameridelphia" 
        elif check_list(lab_based_list, host):
            return "Lab based" 
        else:
            print(f"Unknown host: {host}")
            return "?"
    


     # Open output files
    output_fasta_file = open(output_fasta_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")
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
        # For now, just using a single location field
        # "region",
        # "country",
        # "division",
        # "city",
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
            phylogroup = "?"
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
            CDS = (0,0)
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
            print()
            print(f"Processing {accession_from_list}")
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
                if split_line[0] == "CDS" and segment != desired_segment and "complement" not in line:
                    CDS = (
                        int(line.split(" ")[1].split("..")[0].replace(">", "").replace("<", "")), 
                        int(line.split(" ")[1].split("..")[1].replace(">", "").replace("<", ""))
                    )
                    curr_CDS = True

                # Extract CDS products and determine which segment they come from
                if "/product" in line and curr_CDS and segment != desired_segment:
                    curr_product = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "glycoprotein" in curr_product:
                        segment = "G"
                    elif "Glycoprotein" in curr_product:
                        segment = "G"
                    elif "G protein" in curr_product:
                        segment = "G"
                    elif "g protein" in curr_product:
                        segment = "G"
                    elif "Protein G" in curr_product:
                        segment = "G"
                    elif "G" == curr_product:
                        segment = "G"
                    else:
                        print(f"ERROR: {curr_product} product not known or not {desired_segment}")



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
                    if "/host" in line:
                        host = host_grouper(line.replace("\n", "").replace("\"", "").split("=")[1])
                    # if "/segment" in line:
                    #     segment = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/country" in line:
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

            # Check if CDS was found else extract CDS from full sequence
            if segment == desired_segment and CDS == (0,0):
                print(f"{accession_from_list} excluded because no CDS was found!\n")
                removed_accessions_count += 1
                continue
            else:
                print(f"CDS found between {CDS}!")
                nucleotide_sequence = nucleotide_sequence[CDS[0]-1 : CDS[1]]
                # Remove duplicate sequences
                if remove_duplicates == "Yes":
                    # Check if nulceotide sequence is a duplicate
                    if nucleotide_sequence in seen_sequences:
                        print(f"{accession_from_list} excluded because duplicate sequence!\n")
                        removed_accessions_count += 1
                        continue
                    else:
                        seen_sequences.append(nucleotide_sequence)

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
            assert segment != "?", "Segment is missing!"

            # Create new metadata line
            new_metadata_line = "\t".join([
                strain,
                virus,
                phylogroup,
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

    print()
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
    remove_duplicates = str(snakemake.params.remove_duplicates)

    # Output files
    fasta_output = str(snakemake.output.fasta_sequences) 
    metadata_output = str(snakemake.output.metadata) 

    read_and_process_accession_list(
        list_of_accessions, 
        fasta_output, 
        metadata_output, 
        length_threshold, 
        desired_segment,
        remove_duplicates,
        )


if __name__ == "__main__":
    main()