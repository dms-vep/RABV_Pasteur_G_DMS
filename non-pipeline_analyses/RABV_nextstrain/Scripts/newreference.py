# Description:
# Python script adapted from https://github.com/nextstrain/dengue/blob/main/phylogenetic/scripts/newreference.py

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Seq
import shutil
import argparse
import sys

def new_reference(referencefile, outgenbank, new_name, gene):
    ref = SeqIO.read(referencefile, "genbank")
    startofgene = None
    endofgene = None
    gene_found = False
    for feature in ref.features:
        if feature.type == "source":
            ref_source_feature = feature
        if feature.type == "CDS":
            a = list(feature.qualifiers.items())[0][-1][0]
            if a == gene:
                gene_found = True
                startofgene = int(list(feature.location)[0])
                endofgene =  int(list(feature.location)[-1])+1
                break


    # If user provides a --gene 'some name' that is not found, error out as this may indicate that
    # the gene name is misspelled or the user may be using the wrong GenBank file.
    if(gene is not None and not gene_found):
        print(f"ERROR: No '{gene}' was found under 'CDS' features in the GenBank file.", file=sys.stderr)
        sys.exit(1) 
    record = SeqRecord(
        feature.location.extract(ref).seq, 
        id=new_name,
        annotations={
            "molecule_type": ref.annotations["molecule_type"]
        }
    )
    source_feature = SeqFeature(FeatureLocation(start=0, end=len(record)), type="source",
                                 qualifiers=feature.qualifiers)
    record.features.append(source_feature)
    CDS = SeqFeature(FeatureLocation(start=0, end=len(record)), type="CDS",
                                 qualifiers=feature.qualifiers)
    record.features.append(CDS)

    # Rename gene features
    record.features[-1].qualifiers["gene"] = new_name
    record.features[-1].qualifiers["locus_tag"] = new_name

    # Write new genbank
    SeqIO.write(record, outgenbank, "genbank")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--reference", required=True, help="GenBank file with reference sequences")
    parser.add_argument("--output-genbank", required=True, help="GenBank new reference file")
    parser.add_argument("--new-name", required=True, help="new name for gene")
    parser.add_argument("--gene", required=True, help="gene name or genome for entire genome")
    args = parser.parse_args()

    if args.gene=='genome':
        shutil.copy(args.reference, args.output_genbank)
        SeqIO.write(SeqIO.read(args.reference, 'genbank'), args.output_fasta, 'fasta')
    else:
        new_reference(args.reference, args.output_genbank, args.new_name, args.gene)

