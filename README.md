# Download and extract data
```
wget https://ftp.ncbi.nlm.nih.gov/genbank/gbvrl{1..297}.seq.gz


# unzip and delete archive
gzip -d *.gz

# unzip and keep archive
gzip -dk *.gz
```


# Read all input files, extract useful data and produce output files per gene
```
python3 01-build-dataset.py
```

# Align 
```
muscle -maxiters 2 -in ../data/1-pol.fasta -clwout ../data/2-pol-aligned.fasta 
```

# Convert fasta to csv
```
python3 fasta_to_csv.py ../data/2-pol-aligned.fasta
```