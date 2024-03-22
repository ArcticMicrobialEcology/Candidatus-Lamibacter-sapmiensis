# *Candidatus* Lamibacter sapmiensis

Bioinformatics workflow used in the article:

> Pessi IS, Delmont T, Zehr JP, Hultman J. Discovery of Eremiobacterota with *nifH* homologs in tundra soil. bioRxiv. v2, March 19 2024. doi: [10.1101/2023.06.30.547195](https://doi.org/10.1101/2023.06.30.547195).


## Contacts

**Igor S Pessi**  
Postdoctoral Researcher  
[E-mail](mailto:igor.pessi@gmail.com)

**Jenni Hultman**  
Principal Investigator  
[E-mail](mailto:jenni.hultman@helsinki.fi)


## Before starting

You will need to have these softwares installed and added to your path:  

* [GNU parallel](https://gnu.org/software/parallel)
* [anvi'o v7.1](https://anvio.org)
* [Prodigal v2.6.3](https://github.com/hyattpd/Prodigal)
* [HMMER v3.3.2](http://hmmer.org)
* [blastp v2.14.0](https://ncbi.nlm.nih.gov/books/NBK52640) 
* [MAFFT v7.520](https://mafft.cbrc.jp)
* [trimAl v1.4](http://trimal.cgenomics.org)
* [IQ-TREE v2.2.2.7](http://www.iqtree.org)
* [DIAMOND v2.1.7.161](https://github.com/bbuchfink/diamond)
* [GTDB-Tk v1.5.0](https://github.com/Ecogenomics/GTDBTk) 
* [pyANI v0.2.12](https://github.com/widdowquinn/pyani)
* [Cutadapt v1.16](https://github.com/marcelm/cutadapt)
* [bowtie2 v2.3.5](https://github.com/BenLangmead/bowtie2)
* [SAMtools v1.10](https://github.com/samtools/samtools)
* [CoverM v0.6.1](https://github.com/wwood/CoverM)
* [R v4.2.2](https://r-project.org), with packages [vegan](https://cran.r-project.org/web/packages/vegan) and [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html)
* [Kaiju v1.9.2](https://bioinformatics-centre.github.io/kaiju)
* [ncbi-genome-download v0.3.3](https://github.com/kblin/ncbi-genome-download)

You will need to define some environmental variables:  

```bash
# Number of jobs to run in parallel; change to the number of cores available in your system
export NUM_THREADS=40 

# Working directory; change to where you want to run the analyses
export WD=~/Lamibacter-sapmiensis

# Create working directory if it doesn't exist
if [[ ! -d ${WD} ]]; then mkdir ${WD}; fi
```

Finally, you will need to get the 796 Kilpisjärvi MAGs from [Pessi *et al*. 2022](https://doi.org/10.1186/s40793-022-00424-2).  
Here we will download them from [FigShare](https://doi.org/10.6084/m9.figshare.19722505), but they are also available in [ENA](https://ebi.ac.uk/ena/browser/view/PRJEB41762).    
Then we will import them to `anvi'o`.  

```bash
mkdir ${WD}/KILPISJARVI_MAGs && cd $_

# Download tarball and unpack
wget https://figshare.com/ndownloader/files/35025601 -O kilpisjarvi_MAGs.tar.gz

# Unpack tarball
mkdir RAW && cd $_
tar zxfv ../kilpisjarvi_MAGs.tar.gz

# Clean up file and contig names
cd ..
mkdir CLEAN

for f in `ls RAW/*.fa.gz | grep -E -o 'K[UW]L-[0-9]+'`; do
  g=`echo ${f} | tr '-' '_'`
  gunzip -c RAW/${f}.fa.gz | sed 's/-/_/g' > CLEAN/${g}.fa
done
```

## Analysis of nifH homologs

Here we will search the 796 Kilpisjärvi MAGs for *nifH* homologs.  
We will use the hidden Markov Model (HMM) K02588 from [KOfam](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907).  
Then we will import the MAGs to `anvi'o`, where we will do the search:  

```bash
cd ${WD}/KILPISJARVI_MAGs
mkdir CONTIGSDB

# List MAGs
ls CLEAN/*.fa | grep -E -o 'K[UW]L_[0-9]+' > MAG-LIST.txt

# Create anvi'o contig databases
cat MAG-LIST.txt |
parallel -I % --bar -j ${NUM_THREADS}  'anvi-gen-contigs-database --contigs-fasta CLEAN/%.fa \
                                                                  --output-db-path CONTIGSDB/%.db 2> CONTIGSDB/%_contigs.log'
```

Now we can do the search:  

```bash
mkdir ${WD}/NIFH_ANALYSIS && cd $_
mkdir HMM

# Run HMM search
cat ${WD}/KILPISJARVI_MAGs/MAG-LIST.txt |
parallel -I % --bar -j ${NUM_THREADS} 'anvi-run-hmms --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                                                     --hmm-profile-dir ${WD}/K02588 2> HMM/%.log'

# Get hits
cat ${WD}/KILPISJARVI_MAGs//MAG-LIST.txt |
parallel -I % --bar -j ${NUM_THREADS}  'anvi-script-get-hmm-hits-per-gene-call --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                                                                               --output-file HMM/%.txt \
                                                                               --hmm-source K02588'
```

In the last step, some will throw an error:  
`Config Error: Your contigs database does not have any HMM hits for the HMM source K02588 :/`  
But that is OK, it just means that no hits were found for that particular MAG.    

Let's now retrieve the *nifH* sequences:  

```bash
# List MAGs with nifH
for f in `cat ${WD}/KILPISJARVI_MAGs/MAG-LIST.txt`; do
  if [ -f HMM/${f}.txt ]; then
    echo ${f}
  fi
done > MAG-LIST.txt

# Get nifH sequences
mkdir FASTA

cat MAG-LIST.txt |
parallel -I % --bar --max-args 1 'anvi-get-sequences-for-hmm-hits --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                                                                  --output-file FASTA/%.faa \
                                                                  --hmm-source K02588 \
                                                                  --get-aa-sequences'
```

The next step is to search publc databases to see to what our *nifH*  sequences are similar to.  

Now we can do the search wth `blastp`:  

```bash
# Concatenate sequences
for f in `cat MAG-LIST.txt`; do
  cat FASTA/${f}.faa
done > nifH_MAGs.faa

# BLAST against nr
blastp -query nifH_MAGs.faa \
       -out nifH_blast_nr.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${NR_PATH} \
       -num_threads ${NUM_THREADS}

# BLAST against RefSeq
blastp -query nifH_MAGs.faa \
       -out nifH_blast_refseq.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${RF_PATH} \
       -num_threads ${NUM_THREADS}

# BLAST against SWISS
blastp -query nifH_MAGs.faa \
       -out nifH_blast_swiss.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${SWISS_PATH} \
       -num_threads ${NUM_THREADS}
```

And, for the phylogenetic analysis, we need to get some references.  

Now let's get a phylogenetic tree:  

```bash
# Concatenate sequences
cat nifH_MAGs.faa ${NORTH_FAA} ${ZEHR_FAA} > nifH_MAGs_north_zehr.faa

# Align
mafft --auto \
      --reorder \
      --thread ${NUM_THREADS} \
      nifH_MAGs_north_zehr.faa > nifH_MAGs_north_zehr.aln.faa

# Trim alignment
trimal -in nifH_MAGs_north_zehr.aln.faa \
       -out nifH_MAGs_north_zehr.aln.trimal.faa \
       -gappyout

# Build tree
iqtree -s nifH_MAGs_north_zehr.aln.trimal.faa \
       -m LG+R10 \
       -alrt 1000 \
       -bb 1000 \
       -T ${NUM_THREADS} \
       --prefix nifH_MAGs_north_zehr_iqtree 
```

## Analysis of MAGs containing nifH homologs

Now we will continue working with a subset of the Kilpisjärvi MAGs that contain *nifH* homologs.  
We will import the *nifH* MAGs to `anvi'o` and annotate them.  

```bash
mkdir ${WD}/NIFH_MAGs && cd $_
mkdir FASTA

# List MAGs with nifH
cat ${WD}/NIFH_ANALYSIS/MAG-LIST.txt > MAG-LIST.txt

# Copy nifH MAGs and also create a concatenated fasta
for f in `cat MAG-LIST.txt`; do
  cat ${WD}/KILPISJARVI_MAGs/CLEAN/${f}.fa > FASTA/${f}.fa
  cat ${WD}/KILPISJARVI_MAGs/CLEAN/${f}.fa >> CONTIGS.fa
done

# Create a unified anvi'o contig database
anvi-gen-contigs-database --contigs-fasta CONTIGS.fa \
                          --output-db-path CONTIGS.db \
                          --project-name nif_MAGs \
                          --num-threads ${NUM_THREADS}

# Find single-copy genes
anvi-run-hmms --contigs-db CONTIGS.db \
              --num-threads ${NUM_THREADS}

# Setup GTDB, KEGG and COG databases in anvi'o
anvi-setup-scg-taxonomy
anvi-setup-kegg-kofams
anvi-setup-ncbi-cogs

# Assign taxonomy to single-copy genes
anvi-run-scg-taxonomy --contigs-db CONTIGS.db \
                      --num-threads ${NUM_THREADS}

# Annotate against KOfam
anvi-run-kegg-kofams --contigs-db CONTIGS.db \
                     --num-threads ${NUM_THREADS}

# Annotate against COG
anvi-run-ncbi-cogs --contigs-db CONTIGS.db \
                   --num-threads ${NUM_THREADS}
```

Now we classify the MAGs using `GTDB-Tk`:  

```bash
gtdbtk classify_wf --genome_dir FASTA \
                   --out_dir GTDB \
                   --extension fa \
                   --cpus ${NUM_THREADS}
```

## Read recruitment analysis

```bash
mkdir ${WD}/MAPPING && cd $_

# List MAGs with nifH
cat ${WD}/NIFH_MAGs/MAG-LIST.txt > MAG-LIST.txt

# Make 'external-genomes' file
printf '%s\t%s\n' name contigs_db_path > external_genomes.txt

for MAG in `cat MAG-LIST.txt`; do
  printf '%s\t%s\n' ${MAG} ${WD}/KILPISJARVI_MAGs/CONTIGSDB/${MAG}.db
done >> external_genomes.txt

# Compute ANI with pyANI
anvi-compute-genome-similarity --internal-genomes external_genomes.txt \
                               --output-dir PYANI \
                               --program pyANI \
                               --num-threads ${NUM_THREADS}
```

All MAGs are already < 95% ANI to each other, in fact maximum pairwise ANI is 84.1%.  
