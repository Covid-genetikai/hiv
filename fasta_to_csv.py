import sys
import pandas as pd
from pathlib import Path


def main(fasta_filename):

    path = Path(fasta_filename)

    sequences = []
    with open(path, "r") as fasta_file:
        sequence = {}
        accession = None
        for line in fasta_file:
            if '>' in line:
                if sequence and accession:
                    sequences.append({
                        "accession": accession,
                        "sgene_nucleotide": "".join(sequence)
                    })

                parts = line.split(" ")
                accession = parts[0][1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    
    df = pd.DataFrame.from_dict(sequences)
    df.to_csv(path.with_suffix(".csv"), index=False)

if  __name__ == "__main__":
    main(sys.argv[1])