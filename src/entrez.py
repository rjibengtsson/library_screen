from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

email = "rebecca.bengtsson@unimelb.edu.au"

def fetch_gbk(accession_number):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def fetch_fasta(accession_number):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record


def fetch_gene_coordinates(genbank_file_path, gene_name=None, feature_type="gene"):
    """Return feature coordinates from a local GenBank file.

    Parameters
    ----------
    genbank_file_path : str or Path
        Path to a local GenBank file.
    gene_name : str, optional
        Filter by the 'gene' qualifier value. Leave as None to return all
        features of the given type.
    feature_type : str
        GenBank feature type to search (default ``"gene"``). Pass e.g.
        ``"5'UTR"``, ``"3'UTR"``, or ``"CDS"`` for other feature types.

    Returns
    -------
    list of dict with keys: feature_type, gene, start, end, strand
    """
    with open(genbank_file_path) as f:
        record = SeqIO.read(f, "genbank")
    results = {}
    for feature in record.features:
        if feature.type != feature_type:
            continue
        name = feature.qualifiers.get("gene", [None])[0]
        if gene_name is not None and name != gene_name:
            continue
        results = {
            "accession": record.id,
            "feature_type": feature.type,
            "gene": name,
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "strand": feature.location.strand,
        }
    return results


def save_record(record, type, fasta_file_path, output_dir):
    """
    Extract fasta sequence from GenBank record and save to file.
    """
    gene_name = record.get(type)
    record_id = record.get("accession")
    start = record.get('start')
    end = record.get('end')

    output_file_path = output_dir / f"{record_id}_{gene_name}.fasta"
    
    # Write the extracted sequence to a new FASTA file
    with open(output_file_path, "w") as f:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            record_seq = record.seq[start:end]
            extracted_record = SeqRecord(seq=record_seq, 
                                         id=f"{record.id}_{gene_name}",
                                         description="")
            SeqIO.write(extracted_record, f, "fasta")
