Download and extract data

```
wget https://ftp.ncbi.nlm.nih.gov/genbank/gbvrl{1..297}.seq.gz


# unzip and delete archive
gzip -d *.gz

# unzip and keep archive. Requires lots of storage!
gzip -dk *.gz
```

Install Muscle 3.8 (for biopython)

```
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar xvf muscle3.8.31_i86linux64.tar.gz 
sudo cp muscle3.8.31_i86linux64 /usr/bin/muscle
sudo chmod +x /usr/bin/muscle
```

Install Muscle5 (for command line tools)
```
wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
sudo cp muscle5.1.linux_intel64 /usr/bin/muscle5
sudo chmod +x /usr/bin/muscle5
```

Install BLAST
```
sudo apt install alien
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-1.x86_64.rpm 
sudo alien -i ncbi-blast-2.11.0+-1.x86_64.rpm
```

1. Build dataset and do BLAST

Run `1-build-dataset-and-blast.ipynb` notebook

Change `gene_name` variable to `integrase` or `transcriptase`

2. Align sequences
```
cd /data/hiv/data/integrase
muscle5 -super5 2-integrase.fasta -output 3-integrase-aligned.fasta
```

3. Build philogenetic tree. Requires AWS multicore VM
```
cd /hiv/tools
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
sudo apt install make gcc -y
make -f Makefile.SSE3.PTHREADS.gcc

mkdir /data/hiv/data/integrase/raxml

./hiv/tools/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 90 -f a -x 860647 -p 860647 -N 2 -m PROTGAMMADAYHOFFX -O -n 4-phyl-tree.tre -s /hiv/data/integrase/3-integrase-aligned.fasta -w /hiv/data/integrase/raxml
```


4. Calc distances and make binary tree
```
python3 -m pip install anytree biopython

# change input file name (tree name)
python3 binaryTreeGen_multiproc.py /hiv/data/integrase/raxml/RAxML_bestTree.4-integrase-phyl-tree.tre
```

5. Convert tree.dot to tree.json
```
sudo apt-get install graphviz 

dot -Txdot_json -o tree.json tree.dot
dot -Tpdf -o tree.pdf tree.dot
dot -Tpng -o tree.png tree.dot
```

6. Make pairs
```
python3 /hiv/tools/make_pairs.py /data/hiv/data/integrase/tree.json
```