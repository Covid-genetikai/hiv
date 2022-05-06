### Download and extract data
```
wget https://ftp.ncbi.nlm.nih.gov/genbank/gbvrl{1..297}.seq.gz


# unzip and delete archive
gzip -d *.gz

# unzip and keep archive
gzip -dk *.gz
```

### Install MAFFT
https://mafft.cbrc.jp/alignment/software/linux.html

```
wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb
sudo dpkg -i mafft_7.490-1_amd64.deb
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

### (optional) Align with MAFFT
```
mafft --amino --thread -1 1-pol.fasta > 2-pol-aligned.fasta
```

### (optional) Align with Muscle5 (not 3.8!)
```
cd ../data
muscle5 -super5 2-pol.fasta -output 2-pol-aligned.fasta 
```

### Convert fasta to csv
```
python3 fasta_to_csv.py ../data/2-pol-aligned.fasta
```

---

### BLAST

Download latest from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-1.x86_64.rpm`
`sudo alien -i ncbi-blast-2.11.0+-1.x86_64.rpm`

### Blastp (search only)
```
makeblastdb -in 1-pol.fasta -dbtype prot
blastp -query query.seq -db 1-pol.fasta -out result.txt -outfmt 6
```



### RAxML
```
cd /hiv/tools
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
sudo apt install make gcc -y
make -f Makefile.SSE3.PTHREADS.gcc

./raxmlHPC-PTHREADS-SSE3 -T 8 -f a -x 860647 -p 860647 -N 2 -m PROTGAMMADAYHOFFX -O -n output_file.tre -s /full/path/input.fasta -w /full/path/output/directory

```