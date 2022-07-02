### Download and extract data
```
wget https://ftp.ncbi.nlm.nih.gov/genbank/gbvrl{1..297}.seq.gz


# unzip and delete archive
gzip -d *.gz

# unzip and keep archive
gzip -dk *.gz
```

### (optional) Install MAFFT
https://mafft.cbrc.jp/alignment/software/linux.html

```
wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb
sudo dpkg -i mafft_7.490-1_amd64.deb
```

### Install Muscle 3.8 (for biopython)
```
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar xvf muscle3.8.31_i86linux64.tar.gz 
sudo cp muscle3.8.31_i86linux64 /usr/bin/muscle
sudo chmod +x /usr/bin/muscle
```

### Install Muscle5

```
wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
sudo cp muscle5.1.linux_intel64 /usr/bin/muscle5
sudo chmod +x /usr/bin/muscle5
```

### Build dataset: read all input files, extract useful data and produce output files per gene
```
python3 01-build-dataset.py
```

### BLAST

Download latest from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-1.x86_64.rpm`
`sudo alien -i ncbi-blast-2.11.0+-1.x86_64.rpm`

### Blastp (search only)
```
makeblastdb -in 1-pol.fasta -dbtype prot
blastp -query query.seq -db 1-pol.fasta -out result.txt -outfmt 6
```

### Join blast output with sequences
3-align-detect notebook


### Slice part of sequences for alignment. 2 lines in fasta = 1 accession and 1 sequence
```
head -n 40000 2-pol.fasta > 2-pol-20000.fasta
```

### (optional) Align with MAFFT
```
mafft --amino --thread -1 1-pol.fasta > 2-pol-aligned.fasta
```

### Align with Muscle5 (not 3.8!)
```
cd ../data
muscle5 -super5 2-pol-20000.fasta -output 2-pol-20000-aligned.fasta 
```

### (optional) Convert fasta to csv
```
python3 fasta_to_csv.py ../data/2-pol-20000-aligned.fasta
```


### RAxML
```
cd /hiv/tools
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
sudo apt install make gcc -y
make -f Makefile.SSE3.PTHREADS.gcc

./../tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 90 -f a -x 860647 -p 860647 -N 2 -m PROTGAMMADAYHOFFX -O -n 3-pol-phyl-tree.tre -s /hiv/mantas/2-pol-20000-aligned.fasta -w /hiv/mantas/raxml20000

```

### Calc distances and make binary tree
```
python3 -m pip install anytree biopython

# change input file name (tree name)
python3 binaryTreeGen_multiproc.py /hiv/mantas/raxml20000/RAxML_bestTree.3-pol-phyl-tree.tre

```

### Install Graphviz
```
sudo apt-get install graphviz 
```

### Convert tree.dot to tree.json
```
dot -Txdot_json -o tree.json tree.dot
```

### Make pairs
```
python3 make_pairs.py /data/hiv/data/pol/tree.json
```