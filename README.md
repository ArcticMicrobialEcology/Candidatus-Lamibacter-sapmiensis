# *Candidatus* Lamibacter sapmiensis

Bioinformatics workflow used in the article:

> Pessi IS, Delmont TO, Zehr JP, Hultman J. Discovery of Eremiobacterota with *nifH* homologues in tundra soil. 2024. Environmental Microbiology Reports 16: e13277. doi: [10.1111/1758-2229.13277](https://doi.org/10.1111/1758-2229.13277).


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
* [SRA Toolkit v2.11.3](https://github.com/ncbi/sra-tools)
* [Cutadapt v1.16](https://github.com/marcelm/cutadapt)
* [bowtie2 v2.3.5](https://github.com/BenLangmead/bowtie2)
* [SAMtools v1.10](https://github.com/samtools/samtools)
* [CoverM v0.6.1](https://github.com/wwood/CoverM)
* [Kaiju v1.9.2](https://bioinformatics-centre.github.io/kaiju)
* [ncbi-genome-download v0.3.3](https://github.com/kblin/ncbi-genome-download)

We can define a few variables to simplify some things:  

```bash
# Number of jobs to run in parallel
# Change to the number of cores available in your system
export NTHREADS=40 

# Working directory
# Change to where you want to run the analyses
export WD='~/Lamibacter-sapmiensis'

# Create working directory if it doesn't exist
if [[ ! -d ${WD} ]]; then mkdir ${WD}; fi
```

Finally, we need to get the 796 Kilpisjärvi MAGs from [Pessi *et al*. (2022)](https://doi.org/10.1186/s40793-022-00424-2).  
Here we will download them from [FigShare](https://doi.org/10.6084/m9.figshare.19722505), but they are also available in [ENA](https://ebi.ac.uk/ena/browser/view/PRJEB41762).    

```bash
mkdir ${WD}/KILPISJARVI_MAGs && cd $_

# Download tarball
wget https://figshare.com/ndownloader/files/35025601 \
     -O kilpisjarvi_MAGs.tar.gz

# Unpack tarball
mkdir RAW && cd $_
tar zxfv ../kilpisjarvi_MAGs.tar.gz
cd ..

# Clean up file and contig names,
# so we don't run into issues donwstream
mkdir CLEAN

for i in `ls RAW/*.fa.gz | grep -E -o 'K[UW]L-[0-9]+'`; do
  o=`echo ${f} | tr '-' '_'`
  gunzip -c RAW/${i}.fa.gz | sed 's/-/_/g' > CLEAN/${o}.fa
done
```

## Analysis of nifH homologs

Now that we have the FASTA files for the 796 Kilpisjärvi MAGs, we will search them for *nifH* homologs.   
Here we will use [anvi'o](https://anvio.org) for most of the analyses, so the first thing that we need to do is to create a `CONTIGS.db` for each MAG.      

```bash
cd ${WD}/KILPISJARVI_MAGs
mkdir CONTIGSDB

# List MAGs
ls CLEAN/*.fa | grep -E -o 'K[UW]L_[0-9]+' > MAG-LIST.txt

# Create anvi'o contig databases
cat MAG-LIST.txt |
parallel -I % -j ${NTHREADS}  'anvi-gen-contigs-database \
                               --contigs-fasta CLEAN/%.fa \
                               --output-db-path CONTIGSDB/%.db 2> CONTIGSDB/%_contigs.log'
```

To search for *nifH* homologs, we will use the hidden Markov Model (HMM) K02588 from [KOfam](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907).  
To do this search inside `anvi'o`, we first have to create a [custom HMM source](https://anvio.org/help/main/artifacts/hmm-source/#user-defined-hmm-sources).  
You can either do this from scratch or you can get this [tarball](https://github.com/ArcticMicrobialEcology/Candidatus-Lamibacter-sapmiensis/blob/main/K02588_anvio.tar.gz) and unpack it in you working directory:   

```bash
cd ${WD}

url='https://github.com/ArcticMicrobialEcology/Candidatus-Lamibacter-sapmiensis/raw/main/K02588_anvio.tar.gz'
wget ${url}
tar zxf K02588_anvio.tar.gz
```

Now we can do the search in `anvi'o`:  

```bash
mkdir ${WD}/NIFH_ANALYSIS && cd $_
mkdir HMM

# Run HMM search
cat ${WD}/KILPISJARVI_MAGs/MAG-LIST.txt |
parallel -I % -j ${NTHREADS} 'anvi-run-hmms \
                              --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                              --hmm-profile-dir ${WD}/K02588 2> HMM/%.log'

# Get hits
cat ${WD}/KILPISJARVI_MAGs//MAG-LIST.txt |
parallel -I % -j ${NTHREADS}  'anvi-script-get-hmm-hits-per-gene-call \
                               --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                               --output-file HMM/%.txt \
                               --hmm-source K02588'
```

In the last step, some MAGs will throw an error:  
`Config Error: Your contigs database does not have any HMM hits for the HMM source K02588 :/`  
But that's OK, it just means that no hits were found for this particular MAG.    

Let's now retrieve the *nifH* homologs:  

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
parallel -I % --max-args 1 'anvi-get-sequences-for-hmm-hits \
                            --contigs-db ${WD}/KILPISJARVI_MAGs/CONTIGSDB/%.db \
                            --output-file FASTA/%.faa \
                            --hmm-source K02588 \
                            --get-aa-sequences'

# Concatenate sequences
for f in `cat MAG-LIST.txt`; do
  cat FASTA/${f}.faa
done > nifH_MAGs.faa
```

To check that these sequences are indeed *nifH* homologs, we will use `blastp` to match them against sequences in the *nr*, *RefSeq* and *Swiss-Prot* databases.   
The simplest way to do this would be to take the file `nifH_MAGs.faa` and run it through the `blastp` [web interface](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CLIENT=web&DATABASE=nr&NCBI_GI=on&PAGE=Proteins&PROGRAM=blastp&QUERY=IDQILETNRIACRFNHSNQKYAFSITFQEECAHVTLVVYGRNLHKHFFYWKLHKQLIDLIANPNDMFFF&END_OF_HTTPGET=Y).  
You can also do this locally if you have a copy of the database:   

```bash
# Set up the reference databases
# Change to where the formatted dabatases are located in your system, for example:
export NR_PATH='~/Lamibacter-sapmiensis/dbs/nr'        # NCBI nr
export RS_PATH='~/Lamibacter-sapmiensis/dbs/refseq'    # NCBI RefSeq
export SP_PATH='~/Lamibacter-sapmiensis/dbs/swissprot' # Swiss-Prot

# Search against nr
blastp -query nifH_MAGs.faa \
       -out nifH_blast_nr.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${NR_PATH} \
       -num_threads ${NTHREADS}

# Search against RefSeq
blastp -query nifH_MAGs.faa \
       -out nifH_blast_refseq.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${RS_PATH} \
       -num_threads ${NTHREADS}

# Search against Swiss-Prot
blastp -query nifH_MAGs.faa \
       -out nifH_blast_swissprot.txt \
       -outfmt "6 qseqid sseqid stitle sscinames pident qcovs evalue bitscore" \
       -db ${SP_PATH} \
       -num_threads ${NTHREADS}
```

Another way to check our *nifH* homologs is with a phylogenetic analysis.
Here we will use reference sequences from the [Zehr lab](https://jzehrlab.com/nifh) and from [North *et al*. (2020)](https://doi.org/10.1126/science.abb6310).  
The *Zehr* sequences are available as a FASTA file [here](https://wwwzehr.pmc.ucsc.edu/Genome879/genome879.fasta), and the accession numbers for the *North* sequences can be found in this Excel file [here](https://www.science.org/doi/suppl/10.1126/science.abb6310/suppl_file/abb6310_tables5.xlsx).  
Once you have obtained the two sets of FASTA files, you can proceed to the phylogenetic analysis:  

```bash
# Set up the reference sequences
# Change to where they are located in your system, for example:
export NORTH_PATH='~/Lamibacter-sapmiensis/dbs/north.faa'
export JZEHR_PATH='~/Lamibacter-sapmiensis/dbs/jzehr.faa'

# Concatenate the Kilpisjärvi, North and Zehr sequences
cat nifH_MAGs.faa ${NORTH_PATH} ${JZEHR_PATH} > nifH_MAGs_north_zehr.faa

# Align with MAFFT
mafft --auto \
      --reorder \
      --thread ${NTHREADS} \
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
       -T ${NTHREADS} \
       --prefix nifH_MAGs_north_zehr_iqtree 
```

## Analysis of MAGs containing nifH homologs

From now on we will work only with the Kilpisjärvi MAGs that contain *nifH* homologs.  
We will bring them to `anvi'o`, this time as a unified `CONTIGS.db`.  
We will then find a set of 71 bacterial and 76 archaeal single-copy genes, which will give us an idea of how complete and redundant they are.  
We will also search these single-copy genes against the *GTDB* database, which will give us an idea of their taxonomic placement.  
And we will also annotate the genes against the *KOfam* and *COG* databases.  

```bash
mkdir ${WD}/NIFH_MAGs && cd $_
mkdir FASTA

# List MAGs with nifH
cat ${WD}/NIFH_ANALYSIS/MAG-LIST.txt > MAG-LIST.txt

# Copy nifH MAGs and also create a concatenated fasta
for f in `cat MAG-LIST.txt`; do
  cat ${WD}/KILPISJARVI_MAGs/CLEAN/${f}.fa >  FASTA/${f}.fa
  cat ${WD}/KILPISJARVI_MAGs/CLEAN/${f}.fa >> CONTIGS.fa
done

# Create a unified anvi'o contig database
anvi-gen-contigs-database --contigs-fasta CONTIGS.fa \
                          --output-db-path CONTIGS.db \
                          --project-name nif_MAGs \
                          --num-threads ${NTHREADS}

# Find single-copy genes
anvi-run-hmms --contigs-db CONTIGS.db \
              --num-threads ${NTHREADS}

# Assign taxonomy to single-copy genes
anvi-run-scg-taxonomy --contigs-db CONTIGS.db \
                      --num-threads ${NTHREADS}

# Annotate against KOfam
anvi-run-kegg-kofams --contigs-db CONTIGS.db \
                     --num-threads ${NTHREADS}

# Annotate against COG
anvi-run-ncbi-cogs --contigs-db CONTIGS.db \
                   --num-threads ${NTHREADS}
```

Finally, to classify the MAGs properly we will use `GTDB-Tk`:   

```bash
gtdbtk classify_wf --genome_dir FASTA \
                   --out_dir GTDB \
                   --extension fa \
                   --cpus ${NTHREADS}
```

## Read recruitment analysis

In this section, we will use read recrutiment analysis to estimate the abundance of the *nifH* MAGs in the metagenomes from which they originated.  
Before doing this, we need to see if the MAGs are sufficiently different from each other; if not, they need to be dereplicated before mapping.  
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
                               --num-threads ${NTHREADS}
```

All MAGs have < 95% ANI with each other, so we don't need to dereplicate them.  

---

Now, the next step is to get the metagenomes from [Pessi *et al*. (2022)](https://doi.org/10.1186/s40793-022-00424-2).  
To do this, follow the instructions here:  

1. [Download raw data from ENA with fasterq-dump](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/01-pre-processing.md#download-raw-data-from-ena-with-fasterq-dump) 
2. [Trim adaptors and do quality filtering with Cutadapt](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/01-pre-processing.md#trim-adaptors-and-do-quality-filtering-with-cutadapt)
3. [Pool NextSeq and NovaSeq runs](https://github.com/ArcticMicrobialEcology/Kilpisjarvi-MAGs/blob/main/01-pre-processing.md#pool-nextseq-and-novaseq-runs)

In the end you should have the trimmed and pooled sequences in the folder `${WD}/MAPPING/POLLED_ILLUMINA`.  
Now we will map the reads to the MAGs:  

```bash
# List samples
ls POOLED_ILLUMINA | cut -f 1 -d '.' | sort | uniq > SAMPLE-LIST.txt

# Build bowtie index
bowtie2-build ${WD}/NIFH_MAGs/CONTIGS.fa \
              contigs

# Map reads with bowtie
cat SAMPLE-LIST.txt |
parallel -I % --bar --max-args 1 'bowtie2 \
                                  -1 POOLED_ILLUMINA/%.R1.fastq \
                                  -2 POOLED_ILLUMINA/%.R2.fastq \
                                  -S %.sam \
                                  -x contigs \
                                  --no-unal'

# Index reads with samtools
cat SAMPLE-LIST.txt |
parallel -I % --bar --max-args 1 'samtools view -F 4 -bS %.sam | samtools sort > %.bam && \
                                  samtools index %.bam && \
                                  rm %.sam'

# Compute MAG abundance
coverm genome --bam-files *.bam \
              --genome-fasta-files ${WD}/NIFH_MAGs/FASTA/*.fa \
              --output-file coverM_rpkm.txt \
              --output-format sparse \
              --methods rpkm \
              --min-read-percent-identity 95 \
              --min-read-aligned-percent 75
```
