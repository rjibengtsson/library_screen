from typing import Dict
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation


email = ""

def fetch_gbk(accession_number) -> Path:
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def fetch_fasta(accession_number) -> Path:
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record


def get_feature_locations(genbank_file_path) -> Dict:

    """Return dictionary of feature and location."""

    feature_types = {}
    count = 0
    exclude_feature_types = ['gene', 'mat_peptide', 'source']

    for record in SeqIO.parse(genbank_file_path, "genbank"):
        for feature in record.features:
            if feature.type not in exclude_feature_types:
                gene = feature.qualifiers.get("gene", [None])[0]
                if gene is None:
                    gene = feature.type
                    start = feature.location.start.position
                    end = feature.location.end.position
                    feature_types[gene] = {"start": start, "end": end}
                if len(feature.location.parts) > 1:
                    for part in feature.location.parts:
                        start = part.start.position
                        end = part.end.position
                        count += 1
                        gene_name = f"{gene}_{count}"
                        feature_types[gene_name] = {"start": start, "end": end}
                else:
                    start = feature.location.start.position
                    end = feature.location.end.position
                    feature_types[gene] = {"start": start, "end": end}
    return feature_types



def fetch_gene_sequence(feature, start, end, fasta_file_path, output_dir) -> Path:
    
    """Fetch gene sequences from FASTA files."""
    output_file_path = output_dir / f"{feature}.fasta"
    with open(output_file_path, "w") as f:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            record_seq = record.seq[start:end]
            extracted_record = SeqRecord(seq=record_seq, 
                                         id=f"{feature}",
                                         description="")
            SeqIO.write(extracted_record, f, "fasta")
