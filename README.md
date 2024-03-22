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


## Table of contents

1. [Before starting](#before-starting)
2. [Analysis of nifH homologs](#analysis-of-nifh-homologs)
3. [Analysis of MAGs containing nifH homologs](#analysis-of-mags-containing-nifh-homologs)
4. [Read recruitment analysis](#read-recruitment-analysis)
5. [Analysis of the Eremiobacterota MAG KWL-0264](#analysis-of-the-eremiobacterota-mag-kwl-0264)


## Before starting

We will need to have these softwares installed and added to your path:  

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

We can define some environmental variables to simplify some things:  

```bash
# Number of jobs to run in parallel
# Change to the number of cores available in your system
export NUM_THREADS=40 

# Working directory
# Change to where you want to run the analyses
export WD='~/Lamibacter-sapmiensis'

# Create working directory if it doesn't exist
if [[ ! -d ${WD} ]]; then mkdir ${WD}; fi
```

Finally, we will need to get the 796 Kilpisjärvi MAGs from [Pessi *et al*. (2022)](https://doi.org/10.1186/s40793-022-00424-2).  
Here we will download them from [FigShare](https://doi.org/10.6084/m9.figshare.19722505), but they are also available in [ENA](https://ebi.ac.uk/ena/browser/view/PRJEB41762):  

```bash
mkdir ${WD}/KILPISJARVI_MAGs && cd $_

# Download tarball
wget https://figshare.com/ndownloader/files/35025601 -O kilpisjarvi_MAGs.tar.gz

# Unpack tarball
mkdir RAW && cd $_
tar zxfv ../kilpisjarvi_MAGs.tar.gz
cd ..

# Clean up file and contig names
mkdir CLEAN

for f in `ls RAW/*.fa.gz | grep -E -o 'K[UW]L-[0-9]+'`; do
  g=`echo ${f} | tr '-' '_'`
  gunzip -c RAW/${f}.fa.gz | sed 's/-/_/g' > CLEAN/${g}.fa
done
```

## Analysis of nifH homologs

Now that we have the FASTA files for the 796 Kilpisjärvi MAGs, we will search them for *nifH* homologs.   
The first thing that we will do is import them to `anvi'o`, which is where most of the analyses will be carried out:  

```bash
cd ${WD}/KILPISJARVI_MAGs
mkdir CONTIGSDB

# List MAGs
ls CLEAN/*.fa | grep -E -o 'K[UW]L_[0-9]+' > MAG-LIST.txt

# Create the anvi'o contig databases
cat MAG-LIST.txt |
parallel -I % --bar -j ${NUM_THREADS}  'anvi-gen-contigs-database --contigs-fasta CLEAN/%.fa \
                                                                  --output-db-path CONTIGSDB/%.db 2> CONTIGSDB/%_contigs.log'
```

To search for *nifH* homologs, we will use the hidden Markov Model (HMM) K02588 from [KOfam](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907).  
To do this search inside `anvi'o`, we have to create a [custom HMM source](https://anvio.org/help/main/artifacts/hmm-source/#user-defined-hmm-sources).  
You can either do this from scratch or you can get this [tarball](https://github.com/ArcticMicrobialEcology/Candidatus-Lamibacter-sapmiensis/blob/main/K02588_anvio.tar.gz) and unpack it in you working directory:   

```bash
cd ${WD}

wget https://github.com/ArcticMicrobialEcology/Candidatus-Lamibacter-sapmiensis/raw/main/K02588_anvio.tar.gz
tar zxf K02588_anvio.tar.gz
```

Now we can do the search in `anvi'o`:  

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

In the last step, some MAGs will throw an error:  
`Config Error: Your contigs database does not have any HMM hits for the HMM source K02588 :/`  
But that's OK, it just means that no hits were found for this particular MAG.    

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

# Concatenate sequences
for f in `cat MAG-LIST.txt`; do
  cat FASTA/${f}.faa
done > nifH_MAGs.faa
```

To check that these sequences are indeed *nifH*, we can use `blastp` to find similar sequences in public databases such as *nr*, *RefSeq* and *Swiss-Prot*.  
The simplest way to do this would be to take the file `nifH_MAGs.faa` and run it through the `blastp` [web interface](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CLIENT=web&DATABASE=nr&NCBI_GI=on&PAGE=Proteins&PROGRAM=blastp&QUERY=IDQILETNRIACRFNHSNQKYAFSITFQEECAHVTLVVYGRNLHKHFFYWKLHKQLIDLIANPNDMFFF&END_OF_HTTPGET=Y).  
We can also do this locally if we have a copy of the database:   

```bash

# Setup the databases
# Change to where they are located in your system
export NR_PATH='~/Lamibacter-sapmiensis/dbs/nr'
export RS_PATH='~/Lamibacter-sapmiensis/dbs/refseq'
export SP_PATH='~/Lamibacter-sapmiensis/dbs/swissprot'

# Search against nr
blastp -query nifH_MAGs.faa \
       -out nifH_blast_nr.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${NR_PATH} \
       -num_threads ${NUM_THREADS}

# Search against RefSeq
blastp -query nifH_MAGs.faa \
       -out nifH_blast_refseq.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${RS_PATH} \
       -num_threads ${NUM_THREADS}

# Search against SWISS
blastp -query nifH_MAGs.faa \
       -out nifH_blast_swissprot.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${SP_PATH} \
       -num_threads ${NUM_THREADS}
```

Another way to check our *nifH* sequences is with a phylogenetic analysis.
Here we will use reference sequences from the [Zehr lab](https://jzehrlab.com/nifh) and from [North *et al*. (2020)](https://doi.org/10.1126/science.abb6310).  
The *Zehr* sequences are available as a FASTA file [here](https://wwwzehr.pmc.ucsc.edu/Genome879/genome879.fasta), and the accession numbers for the *North* sequences can be found in this Excel file [here](https://www.science.org/doi/suppl/10.1126/science.abb6310/suppl_file/abb6310_tables5.xlsx).  
Once we have obtained the two sets of FASTA files, we can proceed to the phylogenetic analysis:  

```bash
# Set up the reference sequences
# Change to where they are located in your system
export NORTH_PATH='~/Lamibacter-sapmiensis/dbs/north.faa'
export JZEHR_PATH='~/Lamibacter-sapmiensis/dbs/jzehr.faa'

# Concatenate the Kilpisjärvi, North and Zehr sequences
cat nifH_MAGs.faa ${NORTH_PATH} ${JZEHR_PATH} > nifH_MAGs_north_zehr.faa

# Align with MAFFT
mafft --auto \
      --reorder \
      --thread ${NUM_THREADS} \
      nifH_MAGs_north_zehr.faa > nifH_MAGs_north_zehr.aln.faa

# Trim alignment with trimAl
trimal -in nifH_MAGs_north_zehr.aln.faa \
       -out nifH_MAGs_north_zehr.aln.trimal.faa \
       -gappyout

# Build tree with IQ-TREE
iqtree -s nifH_MAGs_north_zehr.aln.trimal.faa \
       -m LG+R10 \
       -alrt 1000 \
       -bb 1000 \
       -T ${NUM_THREADS} \
       --prefix nifH_MAGs_north_zehr_iqtree 
```

## Analysis of MAGs containing nifH homologs

Now we will continue working only with the Kilpisjärvi MAGs that contain *nifH* homologs.  
We will import them to `anvi'o`, this time as a unified `CONTIGS.db`.  
We will then search for a set of 71 bacterial and 76 archaeal single-copy genes, which will be searched against the *GTDB* database.  
And coding sequences will be annotated against the *KOfam* and *COG* databases:  

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

# Setup GTDB, KOfam and COG databases in anvi'o
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

To classify the MAGs, we will use `GTDB-Tk` and the GTDB r202.  
First, see [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data) how to setup the reference data.  
Then let's run `GTDB-Tk`:  

```bash
gtdbtk classify_wf --genome_dir FASTA \
                   --out_dir GTDB \
                   --extension fa \
                   --cpus ${NUM_THREADS}
```

## Read recruitment analysis

In this section, we will use read recrutiment analysis to estimate the abundance of the *nifH* MAGs in the metagenomes from which they originated.  
Before doing this, we need to see if the MAGs are sufficiently different from each other; if not they need to be dereplicated before mapping.  
We can do this with `pyANI` inside `anvi'o`:  

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
anvi-compute-genome-similarity --external-genomes external_genomes.txt \
                               --output-dir PYANI \
                               --program pyANI \
                               --num-threads ${NUM_THREADS}
```

All MAGs have < 95% ANI with each other, so we don't need to dereplicate them.  
