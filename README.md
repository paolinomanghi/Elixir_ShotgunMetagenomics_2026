# Elixir_ShotgunMetagenomics_2026

- [Hands-on n.0 - Understading how to work on a public server using Conda](#Hands-on-n.0---Understading-how-to-work-on-a-public-server-using-Conda)
- [Hands-on n.1 - Preprocessing of standard metagenomic data](#Hands-on-n.1---Preprocessing-of-standard-metagenomic-data)
- [Hands-on n.2 - Taxonomic profiling: quantifying which species and taxa are there](#Hands-on-n.2---Taxonomic-profiling-quantifying-which-species-and-taxa-are-there)
- [Hands-on n.3 - Scripting in R to compare the results](#Hands-on-n.3---Scripting-in-R-to-compare-the-results)
- [Hands-on n.4: functional profiling at the community level using HUMAnN 4](#Hands-On-4-functional-profiling-at-the-community-level-using-HUMAnN-4)
- [Hands-on n.5 - Metagenome assembly and binning](#Hands-on-n.5---Metagenome-assembly-and-binning)
  * [APPROACH n.1: THE METAGENOMIC ASSEMBLY PROTOCOL STEP BY STEP](#APPROACH-1-THE-METAGENOMIC-ASSEMBLY-PROTOCOL-STEP-BY-STEP)
  * [APPROACH n.2: NEXTFLOW PIPELINES](#APPROACH-2-NEXTFLOW-PIPELINES)
- [Hands-on n.6 - Genome annotation](#Hands-on-n.6---Genome-annotation)

# Hands-on n.0 - Understading how to work on a public server using Conda
## Step n.0: log in into your machine and explore the configuration
```
ssh YOUR-NAME@212.189.202.106
```
In my case, this is:
```
ssh cdonati@212.189.202.106
```

Press Enter
Check our your current location
```
pwd
```

Did it return /data/YOURNAME/?

This is a public server, it means that several people can work in this environment without step on each other toes:
multiple programs must be installed at the same time, and it must be possible to install new ones without intefering with the existing program.
Personnel, for most, do not have root privilegies, so being able to install software packages without compiling is also a requirement.

There are several ways to achieve this situation, the most widely adopted are Docker, NextFlow, and Anaconda. These are different instruments, but they achieve 
similar goals: having an infrastructure that coordinates multiple softwares with multiple version able so that an entire lab can work simultanously.

We'll use Anaconda, as it is the one requiring less specific knowledge.

## WHAT IS CONDA:
* [A GUIDE ON ANACONDA](https://www.anaconda.com/docs/getting-started/main)

Briefly, Anaconda is a widely used open-source platform and distribution for the Python and R programming languages. Purpose: It is designed to simplify package management and deployment for data science, machine learning, and AI projects.

Why it’s used: Instead of downloading thousands of data analysis tools one by one, users install Anaconda, which comes pre-loaded with popular libraries (like NumPy and Pandas) and tools (like Jupyter Notebooks).

First, when you log into the server, you are already in a Conda Environment:
```
Last login: Tue Jun 16 10:17:30 2026 from 151.38.154.153
(base) cdonati@mbc1:~$
```

You see the *(base)* flag? It means that the system recognize your conda ecosystem.
You can deactivate this:

Type:
```
conda deactivate
```

Naturally you may have different anaconda platforms. For instance, we have set up an independent anaconda environment:
```
source /data/anaconda2026/.conda
```

Now, see the (base) string in the bottom-left corner of the screen ?
Type:

```
which python
which conda
```

It should return:
```
(base) cdonati@mbc1:~$ which python
/data/anaconda2026/bin/python

(base) cdonati@mbc1:~$ which conda
/data/anaconda2026/bin/conda
```

## Step n.1: check if your environments are all set up and ready
```
conda info --envs
```

You should see:
```
(base) cdonati@mbc1:/data/cdonati$ conda info --envs
# conda environments:
#
base                  *  /data/anaconda2026
assembly                 /data/anaconda2026/envs/assembly
genome_annotation        /data/anaconda2026/envs/genome_annotation
humann4                  /data/anaconda2026/envs/humann4
preprocessing            /data/anaconda2026/envs/preprocessing
r_notebooks              /data/anaconda2026/envs/r_notebooks
tax_profiling            /data/anaconda2026/envs/tax_profiling
```

Did it work ?

## HOW DO WE RECREATE ANYTHING OF THAT?
In case it doesn't work, the following part can be used to recreate the same exact set of environment from scratch. As you'll see, it is very easy to install all Conda environments.

First ensure that you have entered the Conda ecosystem:
```
conda deactivate
source /data/anaconda2026/.conda
```

Then you can execute the following steps:
0) Be sure the system knows where to search for correct program versions:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
 
1) To create the environment for Data Preprocessing, type:
```
conda create -n preprocessing -c bioconda -c conda-forge trimmomatic bowtie2 samtools
```

2) To create the environment for the Assembly (the genome reconstruction), type:
```
conda create -n assembly -c bioconda -c conda-forge megahit bowtie2 metabat2 checkm2 samtools biopython

conda activate assembly
checkm2 database --download --path /data/CheckM2_database/
conda deactivate
```

3) To create one environment for MetaPhlAn and for Kraken/Bracken (taxonomic profiling), type:
```
conda create -n tax_profiling -c bioconda -c conda-forge metaphlan kraken bracken
mkdir -p /data/metaphlan_databases

metaphlan --install --db_dir /data/metaphlan_databases/ -x mpa_vJan21_CHOCOPhlAnSGB_202103
metaphlan --install --db_dir /data/metaphlan_databases/ -x mpa_vOct22_CHOCOPhlAnSGB_202403

wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
mkdir -p /data/kraken_DB && tar -xvzf k2_standard_08gb_20250402.tar.gz -C /data/kraken_DB

cd /data/
git clone https://github.com/jenniferlu717/KrakenTools.git
chmod +x KrakenTools/*
```

4) To create one environment for running jupiter notebooks with analyses in R, type:
```
conda create -n r_notebooks -c bioconda -c conda-forge r-base=4.3 jupyter r-irkernel r-microeco r-tidyverse=2.0.0
```

Microeco must be installed directly in R:
```
(r_notebooks) cdonati@mbc1:~# R

R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-conda-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Previously saved workspace restored]

> install.packages("microeco", dependencies = TRUE)
```



5) To create one environment for HUMAnN4 (functional profiling of communities), type:
```
conda create -n humann4 -c bioconda python=3.12
conda activate humann4
conda config --add channels biobakery
conda install humann=4.0.0a1 -c biobakery -c bioconda -c conda-forge
conda install metaphlan=4.1.1 bowtie2 -c bioconda

humann_config --update database_folders nucleotide /data/humann_databases/chocophlan/
humann_config --update database_folders protein /data/humann_databases/uniref/
humann_config --update database_folders utility_mapping /data/humann_databases/utility_mapping/
```

Normally, you should execute the last three commands AFTER:
```
#### humann_databases --download chocophlan full humann_databases
#### humann_databases --download uniref uniref90_ec_filtered_diamond humann_databases
#### humann_databases --download utility_mapping full humann_databases
```
DON'T: we did it already.

6) To create one environment for Genomic Annotations (using the software Bakta), type:
```
conda create -n genome_annotation -c bioconda bakta phylophlan

bakta_db download --output /data/bakta_database/ --type light
```

## HOW DOES THIS WORK PRACTICALLY?


To understand how this works, Type:
```
conda activate tax_profiling
```

See the (tax_profiling) string in the bottom-left corner?
Now type:

```
which metaphlan; which kraken; which bracken
conda deactivate
which metaphlan; which kraken; which bracken
```

As you can see, the programs inside each environment are protected, meaning that they are visible only in that environment to not interfere with other installations.

# Hands-on n.1 - Preprocessing of standard metagenomic data
## Step n.1: raw data pre-processing on fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata
In this step, we just download a toy fastq sample. Normally, this step may take a few weeks!

```
mkdir ~/preprocessing
cd ~/preprocessing

wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
unzip input.zip

rm input.zip
cd input
```

## Step n.2: Define variable "s" with the sampleID and run TRIMMOMATIC
```
conda activate preprocessing
s="seq"
```
Conda has the ability to install part of the relevant files it needs. For instance, the software Trimmomatic, which is needed (among other things) to
eliminate Illumina adapters from raw fastq, is installed toghether with all revelant Illumina adapters. This can create some problems, sometimes, when 
one is not familiar with Anaconda package structure. Type:

```
ls /data/anaconda2026/envs/preprocessing/share/trimmomatic/adapters/
```
See multiple fasta files with different type of adapters? This reflect and adapts to different type of sequencing procedure applied, allowing you to use Trimmomatic in
different experimental settings.

Now set the correct adapter sequences for this project:
```
truseq_adap="/data/anaconda2026/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
```

Now launch Trimmomatic to 1) remove adapters 2) perform quality control
```
trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}1.fastq ${s}2.fastq \
${s}_filtered_1.fastq ${s}_unpaired_1.fastq ${s}_filtered_2.fastq ${s}_unpaired_2.fastq \
ILLUMINACLIP:${truseq_adap}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75
```
You should see something like:

```
TrimmomaticPE: Started with arguments:
 -threads 8 -phred33 -trimlog seq_trimmomatic.log seq1.fastq seq2.fastq seq_filtered_1.fastq seq_unpaired_1.fastq seq_filtered_2.fastq seq_unpaired_2.fastq ILLUMINACLIP:/mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75
ILLUMINACLIP: Using adapter file from user-specified absolute path: /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE-2.fa
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 42473 Both Surviving: 41006 (96.55%) Forward Only Surviving: 1259 (2.96%) Reverse Only Surviving: 178 (0.42%) Dropped: 30 (0.07%)
TrimmomaticPE: Completed successfully

```
While normally ignored, stdout prompts like this one can communitcate you a lot of information, establishing a bridge in between your approach and the developer mind (who has, 
typically, conceived the prompt messages to verify the meaningfullness of the software passages)

But what did we produced ?

```
for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done
```

Now that we have cleaned our reads, we can focus on some biological aspect of the project: metagenomic samples can be host-derived or not, for instance they can be associated with a 
human being, a dog, or a plant. In all these cases the characterization of associated microbial communities benefit from a step of removal of the host genome. The host genome can be present in small traces, like in the case of the humann intestine, or it can be constitute large part of the metagenomic sequences (like for instance in the case of dog saliva metagenomic samples).

Environmental samples (water, soil, air) are not associated to a host. In principle, this implies that is not necessary to exclude the host. However, practically ANY sample should be, in addition to any preprocessing phase performed, treated as a human-derived sample, for the simple reason that this will remove putative human DNA contamination in ìntroduced during the sample 
preparation.

We'll now see how to perform the removal of the human genome from the previously polished sequences.

## Step n. 3: Aligning reads against the humann genome, and retrieve only what does not align.

First of all it is common practive to chose a single instance of the human genome. For this purpose we'll chose a modern one, sequenced using long-reads-based technologies
(the 'thelomer-to-thelomer'). A copy of this genome can be find at: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0.fna.

This, like any other host-associated genome used in this phase, is normally retrieved from NCBI as it is important to be able to demonstrate the origin of the genomes used.
We'll use bowtie2 to align our DNA sequences the human genome. The first step is thus to create the bowtie2 indexes to perform such a mapping.
As it is a timy process, we have predisposed the ready-for-use bowtie2 indexes already. 

Let's thus assign the bowtie2 indexes to a variable:
```
human_gen_bowtie_indexes="/data/human_genome/GCF_009914755.1_T2T-CHM13v2.0"
```

Now we'll just align the fastq against the humann h
Run bowtie alignment against the human genome:
```
## VERSION 4 HOURS LONG:
## bowtie2-build ${human_gen_bowtie_indexes/_genomic/.fna} /data/human_genome/GCF_009914755.1_T2T-CHM13v2.0 ## DON'T RUN IT: IT TAKES A FEW HOURS TO BE EXECUTED

bowtie2 -x ${human_gen_bowtie_indexes} -1 ${s}_filtered_1.fastq -2 ${s}_filtered_2.fastq -S ${s}.sam --very-sensitive-local -p 8
```
See something like:

```
41006 reads; of these:
  41006 (100.00%) were paired; of these:
    41006 (100.00%) aligned concordantly 0 times
    0 (0.00%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    41006 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    41006 pairs aligned 0 times concordantly or discordantly; of these:
      82012 mates make up the pairs; of these:
        82007 (99.99%) aligned 0 times
        3 (0.00%) aligned exactly 1 time
        2 (0.00%) aligned >1 times
0.01% overall alignment rate
```
? While often disregarded, these stats are among the most important information you may want to store!
Next, we must operate on the .bam file produced by bowtie2 (Remember? We still want a clean and usable fastq!)

```
samtools view -bS ${s}.sam > ${s}.bam
samtools view -b -f 12 -F 256 ${s}.bam > ${s}.bothunmapped.bam
samtools sort -n -m 5G -@ 2 ${s}.bothunmapped.bam -o ${s}.bothunmapped.sorted.bam
samtools fastq ${s}.bothunmapped.sorted.bam -1 >(gzip > ${s}_filtered.final_1.fastq.gz) -2 >(gzip > ${s}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n
```
Look at the last prompt:
```
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 82002 reads
```

With the last command, we have split the bam file resulting from alignment of both Reverse and Forward reads into two files: the program claims
we have obtained 82002 reads. Let's do a brief check. We count the lines of the files produced:

```
zcat seq_filtered.final_1.fastq.gz | wc
zcat seq_filtered.final_2.fastq.gz | wc
```
The result is in both cases 164004: as one read occupies EXACTLY 4 ROWS, the numbers add up.

# Hands-on n.2 - Taxonomic profiling: quantifying which species and taxa are there
The following part covers the topic that is mostly what people often refer to when tehy name "metagenomics".
In reality, this step covers two highly important but distinct aspects of the characterization of a given microbial community: 
1) Detecting and naming species that are found in the community ("which species are there")
2) Quantifying them (i.e. in relative terms).

Why relative terms?
Taxonomic profiling is, by definition, bound to assign a quantity to a given taxon on the basis of what has been sequenced. 
However, what has been sequenced is, by definition, the result of a set of biases that will impact the final results, in short:

1) The total number of clean reads obtained
2) The general structure of the community.

Why the general structure of the community? Because metagenomic data suffer from a mathematical problem called "compositionality", meaning that
more abundant species tend, by the chemistry of the sequencing machinery, to produce more DNA, and consequently to attract more nucleotide and more
sequencing signallings in the sequencing phase. The result? The species that are more abundant tend to be sequenced more, and result to be more abundant that
they actually are. But they were already more abundant! Therefore sequencing tend to exacerbate the community structure, with more abundant species being
exaggerately more abundant, and lowly abundant species being found often at negligible abundances.

To perform taxonomic profiling, we will use two different methodologies.
The first, implemented in the software MetaPhAn, is called marker-based. It works by assessing whether metagenomic reads map, by the entirety of their length, on 
marker-genes: these marker genes are recorded in pre-constituted databases, and are species-specific. A species presence and abundance is therefore quantified by
the number of marker genes that are found in the community of interest (after a step of normalization accounting the lenght of the genome of the species in question and 
the total number of species-specific marker genes that are available). 

The second approach, implemented in the software Kraken, is called k-mer-based. It works by assessing the presence of k-mers both in the query reads and in a database of genomes
which taxonomy is known. K-mers are combinations of nucleotides strings of lenght 7, 13, 21 or a variable lenght. The software database is constituted in such a way that each k-mer 
points directly to a genome in the taxonomic tree: therefore a read is decomposed into k-mers, and the read k-mers are searched directly in this database.

Which approach is the best?
The approaches display trade-offs and peculiarities.

1) During database creation, MetaPhlAn must align each gene against each other, to ensure the marker-specifity. This massive computation cannot be performed by an individual user,
therefore MetaPhlAn is restricted to the environments that are fully covered by its proprietary database. For human and mice it is considered slightly more accurate.
2) Kraken databases can be easily created for almost any environment, which make it the standard choice for environmental and non-common metagenomes.
3) MetaPhlAn integrates its own taxonomy, while Kraken must be coupled with some pre-existing taxonomy (either pure NCBI of GTDB). While there is no strict "better" in these two approaches,
they remain different: MetaPhlAn has produced important advacements in amending NCBI taxonomy for this reason. Kraken is obliged to adhere to NCBI taxonomy, but remains the sole
transportable to fully uncharacterized communities.
4) To date, MetaPhlAn and Kraken database are similar in disk occupancy and consider a similar number of species (~120000). However MetaPhlAn stores information for 4M genomes, while Kraken stores in the typical use-case 1 genome per species. Therefore in the future MetaPhlAn risk to be inapplicable in environmental settings, while Kraken risks not be scalable to a relevant
number of species.

Conclusion? By most viewpoints, they leverage the same genomic databases and deply a similar species diversity, so to some extent they result comparable. It is recommandable to know both 
and apply the most suitable on a per-case basis.

## Step n.1: Perform MetaPhlAn profiling

We will start by activating the right environment. We have previously created it, is called *tax_profiling*:
```
mkdir metaphlan_taxonomic_profiling
cd metaphlan_taxonomic_profiling

conda deactivate
conda activate tax_profiling
```

## **Step n.2: download metagenomic samples**

We will perform the taxonomic profiling of five samples from human-associated microbiome. Spefically, we will analyzed
microbiomes from:

1) stool
2) buccal mucosa
3) tongue dorsum
4) gingival plaque
5) posterior fornix

Let's launch the command to download these samples:

```
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz
```

## **Step n.3: Run MetaPhlAn 4**

The first thing is to take a brief look at MetaPhlAn parameters. As for any other versatile software, there are numerous paramaters, but we will need only a few.
```
metaphlan -h
```
See ?
```
usage: metaphlan --input_type {fastq,fasta,mapout,sam} [--force] [--db_dir METAPHLAN_DB] [-x INDEX] [--mapout FILE_NAME]
                 [--min_mapq_val MIN_MAPQ_VAL] [--no_map] [--tmp_dir ] [--bt2_ps BowTie2 presets] [--bowtie2_exe BOWTIE2_EXE]
                 [--bowtie2_build BOWTIE2_BUILD] [--tax_lev TAXONOMIC_LEVEL] [--min_alignment_len ] [--ignore_eukaryotes] [--ignore_bacteria]
                 [--ignore_archaea] [--ignore_ksgbs] [--ignore_usgbs] [--stat_q ] [--perc_nonzero ] [--ignore_markers IGNORE_MARKERS]
                 [--avoid_disqm] [--stat ] [-t ANALYSIS TYPE] [--nreads NUMBER_OF_READS] [--pres_th PRESENCE_THRESHOLD] [-o output file]
                 [--sample_id_key name] [--use_group_representative] [--sample_id value] [-s sam_output_file] [--CAMI_format_output]
                 [--skip_unclassified_estimation] [--biom_format_output] [--biom_mdelim mdelim] [--profile_vsc] [--vsc_out VSC_OUT]
                 [--vsc_breadth VSC_BREADTH] [--long_reads] [--max_gcsd MAX_GCSD] [--minimap2_exe MINIMAP2_EXE] [--minimap2_ps minimap2 presets]
                 [--nbases NUMBER_OF_BASES] [--split_reads] [--split_readlen SPLIT_READLEN] [--nproc N] [--subsampling SUBSAMPLING]
                 [--mapping_subsampling] [--subsampling_seed SUBSAMPLING_SEED] [--subsampling_output SUBSAMPLING_OUTPUT]
                 [--subsampling_paired SUBSAMPLING_PAIRED] [-1 FORWARD_READS] [-2 REVERSE_READS] [--install] [--offline] [--force_download]
                 [--read_min_len READ_MIN_LEN] [--verbose] [-v] [-h]
                 [INPUT_FILE]

 MetaPhlAn version 4.2.4 (21 Oct 2025): 
 METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.

 AUTHORS: Aitor Blanco-Miguez (aitor.blancomiguez@unitn.it), Francesco Beghini (francesco.beghini@unitn.it), Nicola Segata (nicola.segata@unitn.it), Duy Tin Truong, Francesco Asnicar (f.asnicar@unitn.it), Claudia Mengoni (claudia.mengoni@unitn.it), Linda Cova (linda.cova@unitn.it)

positional arguments:
  INPUT_FILE            the input file can be:
                        * a fastq file containing metagenomic reads
                        OR
                        * a BowTie2 produced SAM file. 
                        OR
                        * a minimap2 produced SAM file (for long reads). 
                        OR
                        * an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run (mapout)
                        If the input file is missing, the script assumes that the input is provided using the standard 
                        input, or named pipes.
                        IMPORTANT: the type of input needs to be specified with --input_type

[...]

```
MetaPhlAn can take multiple file types as input, and outputs different files. In short, when run on metagenomic reads,
the procedure entails:

1) mapping reads (using Bowtie2) against the marker gene database
2) calculus of average coverage per species
estimated by the coverage on the marker genes
3) pooling of the coverages and estimation of species relative abundance
4) grouping of the taxonomic level, starting from species-level, up to the kingdom

We will start setting MetaPhlAn database paths and running it on the gingival plaque metagenomic sample:
```
mpa_db="/data/metaphlan_databases/"
db_version="mpa_vJan21_CHOCOPhlAnSGB_202103"

s="SRS014476-Supragingival_plaque"

metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --nproc 8 --db_dir ${mpa_db} --index ${db_version}

for s in SRS014459-Stool SRS014472-Buccal_mucosa SRS014470-Tongue_dorsum SRS014494-Posterior_fornix SRS014476-Supragingival_plaque; do
    cp /data/pmanghi/metaphlan_taxonomic_profiling/${s}_profile.txt ${s}_profile.txt; done
```

***** DON'T RUN THE FOLLOWING COMMANDS *****
```
#### s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
   --nproc 8 --db_dir ${mpa_db} --index ${db_version}
#### s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
   --nproc 8 --db_dir ${mpa_db} --index ${db_version}
#### s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
  --nproc 8 --db_dir ${mpa_db} --index ${db_version}
#### s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
  --nproc 8 --db_dir ${mpa_db} --index ${db_version}
```

...AND merge the profile in a single table
```
merge_metaphlan_tables.py *_profile.txt | grep -P "clade_name|UNCLASSIFIED|t__" > metaphlan_table.tsv
```

## **Step n.4: Perform Kraken + Bracken taxonomic profiling (of the same samples!)**

Kraken is normally coupled with Bracken. The two softwares are considered companions and have been designed by the same developers.

Although considered among the top-performers, metagenomic people are aware of a potential kraken limitation: mapping onto the whole genomic sequence, k-mers may be ambigously located on two highly-similar species. To solve this problem, Kraken assigns reads that may ambigously map against species A and species B their Lowest Common Ancestor (LCA). As species A and B are similar, their LCA will likely be their genus: Kraken thus display a bias, in that it assigns reads to superior taxonomic level more easily.

Bracken is, in words, the piece of software designed to address this issue: Bracken uses Baysian statistics to redistribute ambigous reads, and as such the trick is hydden in the correct estimation of the priors. Each species has two key-numbers:

1) the number of unique k-mers (which is a charcteristic of the database)
2) the number of uniquely assigned reads.

Clearly, the ** random ** distribution of reads between the two species is not derived from having the same number of reads among the two, but having the same proportion of
Number of reads / Number of unique K-mers

Or, in other words, species that have very little unique k-mers in their genomes, are expected to show a greater abundance if we assumed that the number of uniquely assigned reads
is the same. In short, Bracken redistribute ambiguous reads based the probability of observing the observed number of non-ambiguous calls, considering the number of unique k-mers for
the species.

Albeit pratically this methodology works, and is the de-facto standard for non-humann communities, precautions must be taken when extensively running Bracken:

1) The resulting number of species with standard parameters is notoriously too-high (based on consideration on human being, in which we know precisely what to expect)
2) In environmental samples, the Bracken rate of False Positives may be perceived as sensitivity (since finding a lot of species in an environment which should have a lot of species seems like the right thing: in reality the great majority of species found in such environment remains uncharacterized also by Kraken).

When running Bracken, it is therefore practically advised to apply an additional filter at > 1000 reads, meaning lowering to zero everything that, after the Baysian statistical algorithm, maintain a count of reads below 1000. Despite this potential pitfall, Bracken remains the standard for profiling environmental or non-common communities (for example marine metagenomes or yet-to-be-surveyed animals).

Let's now create the folder and run the two softwares:
```
mkdir -p ~/Kraken_Bracken_taxonomic_profiling
cd ~/Kraken_Bracken_taxonomic_profiling
```

Let us have a look a Kraken parameters:
```
kraken2 -h
```
See ?
```
Usage: kraken2 [options] <filename(s)>

Options:
  --db NAME               Name for Kraken 2 DB
                          (default: none)
  --threads NUM           Number of threads (default: 1)
  --quick                 Quick operation (use first hit or hits)
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "-" will
                          suppress normal output
  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                          in [0, 1].
  --minimum-base-quality NUM
                          Minimum base quality used in classification (def: 0,
                          only effective with FASTQ input).
  --report FILENAME       Print a report with aggregrate counts/clade to file
  --use-mpa-style         With --report, format report output like Kraken 1's
                          kraken-mpa-report
  --report-zero-counts    With --report, report counts for ALL taxa, even if
                          counts are zero
  --report-minimizer-data With --report, report minimizer and distinct minimizer
                          count information in addition to normal Kraken report
  --memory-mapping        Avoids loading database into RAM
  --paired                The filenames provided have paired-end reads
  --use-names             Print scientific names instead of just taxids
  --gzip-compressed       Input files are compressed with gzip
  --bzip2-compressed      Input files are compressed with bzip2
  --minimum-hit-groups NUM
                          Minimum number of hit groups (overlapping k-mers
                          sharing the same minimizer) needed to make a call
                          (default: 2)

```

We are now going to run kraken on the same set of samples as before:
```
for s in SRS014459-Stool.fasta.gz SRS014470-Tongue_dorsum.fasta.gz SRS014472-Buccal_mucosa.fasta.gz SRS014476-Supragingival_plaque.fasta.gz SRS014494-Posterior_fornix.fasta.gz;
do
kraken2 --db /data/kraken_DB/ --threads 8 --report `basename ${s%.fasta.gz}`.kraken2_report.txt --output `basename ${s%.fasta.gz}`.kraken2_output.txt ../metaphlan_taxonomic_profiling//${s};
done
```

See?
```
Loading database information... done.
20000 sequences (1.97 Mbp) processed in 0.274s (4383.7 Kseq/m, 430.85 Mbp/m).
  6529 sequences classified (32.65%)
  13471 sequences unclassified (67.36%)
Loading database information... done.
20000 sequences (1.90 Mbp) processed in 0.206s (5820.9 Kseq/m, 551.89 Mbp/m).
  4430 sequences classified (22.15%)
  15570 sequences unclassified (77.85%)
Loading database information... done.
20000 sequences (1.95 Mbp) processed in 0.206s (5827.7 Kseq/m, 566.89 Mbp/m).
  6786 sequences classified (33.93%)
  13214 sequences unclassified (66.07%)
Loading database information... done.
20000 sequences (1.85 Mbp) processed in 0.134s (8956.7 Kseq/m, 827.78 Mbp/m).
  4803 sequences classified (24.02%)
  15197 sequences unclassified (75.98%)
Loading database information... done.
20000 sequences (1.90 Mbp) processed in 0.127s (9420.9 Kseq/m, 895.69 Mbp/m).
  8515 sequences classified (42.58%)
  11485 sequences unclassified (57.42%)
```

As said, we are now going to correct the ambiguity bias in Kraken using Bracken:
```
for s in SRS014459-Stool.fasta.gz SRS014470-Tongue_dorsum.fasta.gz SRS014472-Buccal_mucosa.fasta.gz SRS014476-Supragingival_plaque.fasta.gz SRS014494-Posterior_fornix.fasta.gz;
do

bracken -d /data/kraken_DB/ -i `basename ${s%.fasta.gz}`.kraken2_report.txt -o `basename ${s%.fasta.gz}`.bracken_abundance.txt -w `basename ${s%.fasta.gz}`.bracken_report.txt -l S -t 150;

done
```

Similarly to what done for MetaPhlAn, we will now merge all the Bracken outputs into a MetaPhlAn-like taxonomy, to make them comparable. We will make use of the suite KrakenTools, which we have previously downloaded from github.

```
for s in *.bracken_report.txt; do  /data/KrakenTools/kreport2mpa.py --display-header -r ${s} -o ${s%.txt}.mpa.tsv; done
```

This last command has created files in a format that we can handle very easily. Type:
```
cat SRS014470-Tongue_dorsum.bracken_report.mpa.tsv
```
See?

```
#Classification SRS014470-Tongue_dorsum.bracken_report.txt
d__Bacteria     1152
d__Bacteria|k__Pseudomonadati   1152
d__Bacteria|k__Pseudomonadati|p__Bacteroidota   880
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia    880
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales   880
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae 880
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella   880
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_melaninogenica      475
d__Bacteria|k__Pseudomonadati|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_histicola   404
d__Bacteria|k__Pseudomonadati|p__Campylobacterota       271
d__Bacteria|k__Pseudomonadati|p__Campylobacterota|c__Epsilonproteobacteria      271
d__Bacteria|k__Pseudomonadati|p__Campylobacterota|c__Epsilonproteobacteria|o__Campylobacterales 271
d__Bacteria|k__Pseudomonadati|p__Campylobacterota|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae   271
d__Bacteria|k__Pseudomonadati|p__Campylobacterota|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter  271
d__Bacteria|k__Pseudomonadati|p__Campylobacterota|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_concisus        271
```

We now want to combine all Bracken table together, like we did for MetaPhlAn:

```
/data/KrakenTools/combine_mpa.py -i *.bracken_report.mpa.tsv -o merged_bracken_table.tsv
sed 's/.bracken_report.txt//g' merged_bracken_table.tsv | grep -P 'Classification|s__' | sed 's/Bacillati/Bacteria/g' | sed 's/Pseudomonadati/Bacteria/g' > bracken_table.tsv
```

# Hands-on n.3 - Scripting in R to compare the results

## Exercises:
First, set the correct conda environment:
```
conda deactivate 
conda activate r_notebooks

cd ~
```

Next, open jupiter notebook:
```
jupyter notebook /data/Jupyter/Alpha_diversity.ipynb --allow-root --no-browser --port=8888 --ip=127.0.0.1
```

Next, open up a new TERMINAL and Type:
```
ssh -N -L 8888:127.0.0.1:8888 cdonati@212.189.202.106
```

Note it holds but doesn't do anything. It means that a tunnel is open.
Then, copy the URL that the server prompt is suggesting you into YOUR BROWSER:

Es.:
```
http://212.189.202.106:8888/tree?token=398cde02036d5c0c4e8162b5e21758c5d7f9fa90dc4eabad
```
press ctrl + c, type y and exit jupyter


# Hands-On 4: functional profiling at the community level using HUMAnN 4

The next tool we are going apply represents one unique alternative to perform genetic/functional analysis in microbiome project. As a general concept, HUMAnN 4 solves a easy-to-understand problem, i.e.: it assigns metagenomic reads to functions in UniRef90.

Practically, there are two things to consider. 
First, it is one of a kind, in the sense that just a few other, less popular tools perform a similar task, while on the contrary a good number of tools exists that quantify or predict genes starting from reconstructed genomes. HUMANnN quantifies gene abundances starting from reads, which makes it sort of an unicum.

Second, HUMAnN employs an algoritmical trick to make this computation faster and more efficient. In short, mapping raw reads against the whole UniRef90 is an expensive task. HUMAnN thus relies on the following pipeline:

1) It runs MetaPhlAn
2) It selects the species found at an abundance ≥0.01%.
3) It filters out the pangenomes of such species (the pangenome of those species is available simply because is part of the MetaPhlAn database. Therefore if it runs MetaPhlAn it can also access its pangenomes). It also filters out the reads that were confidently mapped by MetaPhlAn.
4) **It maps the reads that are confidently mapped by MetaPhlAn to a certain species only against the pangenome of THAT species**
5) It maps the remaining raw reads against UniRef90 in translated research.

What does HUMAnN do with the results from this mapping?
HUMAnN has two key-outputs:

1) the UniRef90 quantification, i.e.: for each protein identified (~20000-40000 per sample), it computes the RPK (read per kilobase), meaning that it natively normalizes the average number of reads hitting a single amino acid position on a protein by the lenght of the protein (NOTE: but it doesn't normalize over the sample sequencing depth)
2) the pathway table: starting from the previous one, HUMAnN runs an algorithm known as "MinPath" on the MetaCyc database, applying a maximum parsimony algorithm: basically, it reconstructs, based on the single protein result, a pathway-based quantification of the genetic repertoire of a given community. 

Last thing to remember: the algorithmical trick of dividing the total reads into i) mappable onto a specific pangenome and ii) not mappable against a specific pangenome, allows one more result: BOTH the tables aforementioned contain, in reality, two different results each:

a) The community-level quantification of a single protein, or a pathway
b) The PER-SPECIES ("stratified") output, were the protein or the pathway total RPKs are stratified by species of origin.

A few last notes of this tool: if you followed up to here, you should have learned that MetaPhlAn, at least natively, does not provide a wonderful performance when speaking about non-human, and in general non-mammalian (i.e. soil, water...) microbiomes. ** Does this mean that HUMAnN cannot run on these metagenomes? No, technically not: it can; the sole difference is that the number of reads that will be conveniently mapped against MetaPhlAn-identified pangenomes will much lower, and likely the total computational time will increase.

## Step n.1: Get into the right directory & install download the necessary files
```
mkdir -p ~/functional_profiling
cd ~/functional_profiling
```

## Step n.2: test that HUMAnN runs properly and have a look at the HUMAnN parameters
```
conda activate humann4

humann_test
humann_config
humann -h
```

## Step n.3: get a sample from EBI
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/096/SRR15408396/SRR15408396.fastq.gz
```

## Step n.4: RUN humann
```
s="SRR15408396"
metaphlan_options="--bowtie2db /data/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202403 --index mpa_vOct22_CHOCOPhlAnSGB_202403 -t rel_ab_w_read_stats"
```

***** DON'T RUN THE FOLLOWING COMMANDS *****
```
#### \humann --input ${s}.fastq.gz --output ${s} --threads 8  --count-normalization RPKs --metaphlan-options "${metaphlan_options}"
#### rm -r ${s}/${s}_humann_temp/
```

It would take 4-5 HOURS:
RUN INSTEAD:
```
mkdir -p ${s}

cp /data/pmanghi/functional_profiling/${s}/${s}_2_genefamilies.tsv ${s}/${s}_genefamilies.tsv
cp /data/pmanghi/functional_profiling/${s}/${s}_4_pathabundance.tsv ${s}/${s}_pathabundance.tsv
```

## Step n.5: Regrouping genes to other functional categories
```
humann_regroup_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_pfam.tsv --groups uniref90_pfam
```

While executing, it will communicate some of information of the conversion from UniRef90 to PFAM:
```
(humann4) cdonati@mbc1:/data/cdonati/functional_profiling# humann_regroup_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_pfam.tsv --groups uniref90_pfam
Loading table from: SRR15408396/SRR15408396_genefamilies.tsv
  Treating SRR15408396/SRR15408396_genefamilies.tsv as stratified output, e.g. ['UniRef90_Q45125', 's__Parabacteroides_distasonis.t__SGB1934']
Loading mapping file from: /data/humann_databases/utility_mapping/map_pfam_uniref90.txt.gz
  This is a large file, one moment please...
Original Feature Count: 165701; Grouped 1+ times: 80349 (48.5%); Grouped 2+ times: 24801 (15.0%)
(humann4) cdonati@mbc1:/data/cdonati/functional_profiling#
```

## Step n.6: Run HUMAnN on a second sample
```
s="SRR15408398"
```

***** AGAIN, DON'T RUN THE FOLLOWING COMMANDS !! *****

```
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/098/SRR15408398/SRR15408398.fastq.gz
#### \humann --input ${s}.fastq.gz --output ${s} --threads 8  --count-normalization RPKs --metaphlan-options "${metaphlan_options}"
#### rm -r ${s}/${s}_humann_temp/
```

RUN INSTEAD:
```
mkdir -p ${s}
cp /data/pmanghi/functional_profiling/${s}/${s}_4_pathabundance.tsv ${s}/${s}_pathabundance.tsv
```

## Step n.7: Merge together community profiles under different ontologies
```
mkdir -p merged

cp SRR15408396/SRR15408396_pathabundance.tsv merged/
cp SRR15408398/SRR15408398_pathabundance.tsv merged/

humann_join_tables -i merged -o merged_pathabundance.tsv --file_name pathabundance
```

You should see something like:
```
(humann4) cdonati@mbc1:/data/cdonati/functional_profiling# humann_join_tables -i merged -o merged_pathabundance.tsv --file_name pathabundance
Gene table created: /data/cdonati/functional_profiling/merged_pathabundance.tsv
(humann4) root@mbc1:/data/cdonati/functional_profiling#
```

You can have a look at this table:
```
(humann4) cdonati@mbc1:/data/cdonati/functional_profiling# less -S merged_pathabundance.tsv
```
Press q


# Hands-on n.5 - Metagenome assembly and binning
This part of the tutorial introduces one of the foundations of shotgun metagenomics, i.e. reference-free metagenomic assembly.

Practically, this topic is centred around reconstructing good-quality draft genomes from metagenomic samples.

Clearly, the procedure as described hydes some hurdles that is worth name:

1) The reconstruction of genomes or chunks of genomes involves, to simplify, 2 main steps: assembly properly called (1), and binning (2)
2) The procedure allows to define genomes of species, but *which species each genome is from* is a meaningfull question that is answered *after* the metagenomic assembly procedure.
3) The procedure allows to detect single genes. This step, a.k.a. genome annotation, is also perormed, usually, later on, but in both cases involves more operations.

All these aspects are normally covered during metagenomic assembly. We will follow, in this tutorial, a dual approach. If you wonder: "Has ever someone put all these pieces into a larger program/wrapper to performed them automatically?" The answer is yes. As the goal of this tutorial is to illustrate some of the key procedures metagenomics involves, we will follow either routes: we will now first show the listed assembly step following the standard protocol step by step; secondly, we will show how these aspects can be all covered via an automatic pipelines, *using NextFlow*. 

## APPROACH 1: THE METAGENOMIC ASSEMBLY PROTOCOL STEP BY STEP
Following the standard ("by hands") procedure, assembly comprises at least four phases:

1) reconstruction of genomic chunks from metagenomic reads
2) remapping of metagenomic reads against the obtained chunks
3) binning
4) quality control.

In brief, the overall procedure works like this: 1) genomic chunks are **assembled**. This literally means that reads that **end** with a specific sequence of, say, 50 nucleotides, are merged to reads that **start** with **those very 50 nucleotides**. Once this step is iterated for hundreds of thousands of times, it results in the reconstruction of a genome chunk. Practically, however, this phase requires vast algorithmical optimization. In fact it must align billions of reads, deciphering different degree of similarity and in general constituting an enormous graph in which nodes (reads) are link via more or less robust branches (overlaps). This is graph is termed De Brujin graph. 

Most notably, the result of a De Brujin graph is a sequence of nucleotides that represents a fraction of a genome. There are however at least 3 relevant details to keep in mind:

1) The De Brujin graph by definition works determining what is the **consensus sequence** of the reconstructed chunk. Think about it: is a logical consequence of determining overlapping regions despite single nucleotide differences. Is this a problem ? It really depends on your goal. Although assembly is generally a desirable starting point, this problem may arise in fields in which assembly seems the most suitable technique: soil metagenomes, for example, are characterized by a high number of uncharacterized taxa, making you think that "assembly is the best approach to start exploring the soil microbiome". However, soil metagenomes are often floaded with micro-niches and strain-heterogenenity, meaning that assembly may overlook substantial intra-species variation. This results in:
2) Assembly produces 1 nucleotide sequence, but you don't really known how important this sequence is in the community: it may be the most important species of the entire ecosystem or a low-abundance contaminant.
3) The produced sequence, termed contigs (for "contigous fragments"), are not genomes. In fact they are rather shorter. Contigs can be as long as half of a genome, or be very short, down to a few hundreds nucleotides.

Conceptually, reconstruction of contigs is followed by binning those contigs into genomes. In fact, binning is simply the procedure of "binning together" contigs from the same microbial species. While this looks intuituive, following from points 1-3 you may understand while this is not trivial. In fact, the reason why we perform also step n. 2 ("remapping of metagenomic reads against the obtained chunks") is that to bin contigs into genomes it is necessary to have the results from this remapping.

Binning means collecting toghether pieces of genome reconstructed from assembly, assuming that they belong to the same species. There are two popular algoritmical strategies to do so, which are normally combined: 

a) frequency of tetra-nucleotides (TNF): each combination of ACTG (256 in total) is counted. The frequency of such combination is informative, meaning that it tend to gather information such as codon usage, restriction sites
b) co-abundant chunks: this part uses the results from the bowtie2 remapping, telling us how deeply are covered genomic chunks. In binning this is highly informative: two moieties from the genome of a single species are clearly present with the same number of copies: n copies each, for n copies of the species. This implies that, mapping metagenomic reads onto all contigs, the two moieties will absorb approximately the same number of reads.

Binning normally combines these two strategies in a single similarity metric, that is then use to perform a hierarchical clustering to tell chunks from the same species apart. These metrics must lean of a certain contig lenght to work properly, reason practically contigs shorter than 1500 nt are normally excluded from any computation focused on genomes.
    
The last step is quality-control: clearly, the above procedure is not granted to produce perfect genomes. The binned contigs may be a) binned wrongly or simply b) too few to form a decent genome. Practically, what is normally done is that bins are generated, and suitable genomes are determined among these by screening them with dedicated softwares.

## Step n.1: Start Metagenomic Assembly

We will start by activating the right environment. We have previously created it, is called *assembly*:
```
mkdir -p ~/assembly
cd ~/assembly

conda deactivate
conda activate assembly
```

## Step n.2: check if software Megahit is present and run it:
```
megahit -h
```

See?
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ megahit -h
MEGAHIT v1.2.9

contact: Dinghua Li <voutcn@gmail.com>

Usage:
  megahit [options] {-1 <pe1> -2 <pe2> | --12 <pe12> | -r <se>} [-o <out_dir>]

  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
    -1                       <pe1>          comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
    -2                       <pe2>          comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
    --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
    -r/--read                <se>           comma-separated list of fasta/q single-end files

Optional Arguments:
  Basic assembly options:
    --min-count              <int>          minimum multiplicity for filtering (k_min+1)-mers [2]

```

We won't run it. It takes a couple of hours. In fact is by far the time-consuming step in assembly.
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz

s="SRR341725"
## MEGAHIT WILL TAKE A FEW HOURS:
## megahit -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -o ${s}.megahit_asm -t 8

## FOR NOW WE CAN COPY THE RESULTS FROM MEGAHIT
mkdir -p ${s}.megahit_asm/

cp /data/pmanghi/assembly/${s}.megahit_asm/final.contigs.fa ${s}.megahit_asm/
cp /data/pmanghi/assembly/${s}.megahit_asm/contigs.fasta ${s}.megahit_asm/
```

Type:
```
less SRR341725.megahit_asm/final.contigs.fa
```

You see:
```
>k99_61317 flag=0 multi=13.0000 len=254
ATAGGAGAAAAGAAAATGGAATTTATAATTGATACAGTAAACTTAGAAGATATAAAAGAAGCAGTAGAATATATGCCAATTGTAGGAGTAACAAGTAATCCTTCAATCGTAAAGAAAACAAGTCCAAAAGACTTCTTTCAACATATGAAAGAAGTAAGAAAAATCATTGGAAAAGAAAGAAGTCTACACATTCAAGTTATTTCTAAAGAATGTGATGAAATCGTTAAAGAAGCCCATCGTATTATTGAAGAAAT
>k99_18561 flag=1 multi=2.0000 len=334
TTCTAAAAACTCAAAAACGTATTTTGAATATTCATTACCATTTCCTCTTACTATATCTATTGTATCTTCATGCTCAGCACCTTTTATTATTTTTAAATAAGAATTATTGTTTACATTCTTTAAATTATTTATAAAGGATTTTGAGTTTTCTATATCTATTAGTGTATCACTATCACCATGTATACAAAGAGTACTTATAGTAGAGTTTTTATTTATTAGATTTATAGGATTACAAGTATCAGTATTACTTTTATATATATACCTTTTTATTAACCTTTTAGATCTATTTGAATTACATTTATTAAAATTTAAAACTCCAGATATGGATATAAAA
>k99_9684 flag=0 multi=1.0000 len=293
TATCTCCCAGCTGTTTTTACAATGCAGTTACCCGTCTCCGCAAAAAAGCCTGTCAGATACCGGATCCAGCCCCCAAAGCCAGCACTCTTGACCTTACATCCCATAAACAGGATGTAGTCCAGATTGCTATAGAACCGGAGTCTTCTCCGGCAGGACTGATCCCAGATAATGGGAACAGACCTATGCACCTTGACAATTCACATACGATTGAAATCGAAGCAGACGGATTGCTTATACGAATGAGTAATGAGATCAAACCGTTACTTCTGAAGATGCTTATGGATACGCTTAAG
>k99_16140 flag=1 multi=2.0000 len=336
TATTGGCGAGGGAGCTGGAGTTGTTGTATTAGAGGAATTAGAACATGCGAAAAAGCGTGGGGCAAAGATTCTTGCAGAGTTAGCCGGATATGGTTCTACCTGTGATGCCTTCCACATCACTTCTCCTGCAGAGGACGGAAGCGGTGCAGCAAAAGCAATGGAACTTGCCATGGAAGAAGCAGGTGTAAAACCGGAAGAAGTAGATTACATTAACGCACATGGAACAGGAACCCACCACAACGATTTATTTGAGACAAGAGCTGTAAGACGTGCCTTTGGTGAGAGCGCAGACCATTTAAAAATGAGTTCTACAAAGTCAATGATCGGACATTTGCT
>k99_13719 flag=1 multi=4.0000 len=378
TCGTATTATAGTCCAACGACTCTTTTGAGTAGATGAAGCGGTGAGCAAACTGACCCTTAGCTAACGAACTGTTGGGGTTATATGCCAGCGTACCCAAAACGTAATCATCCGGATTACCGTCTGTTGTATACTCGCTGACTTCCGGAATACCGGTGCGGTTGATATCGAACCATGCATCCCAAGCCTGGCAACGGGTGCTGGCAATCCATTTCTGGGTAAGGATACATTTCAGCATACTTTCTTCGTCGGTGTCTTGGAATTCGTAAGGTTTTCCTGCTGCAATAAATTCATCGATTTTATTGTCATGGTTCCATCTGGCGAAAGAGAGTTCAACGGCATCGTCATAGTGAGTTTTTGCATTACCTTTGTCTCCCAGTC
>k99_30666 flag=1 multi=3.0000 len=364
ATGTACATGTACAGCGCCCTTCCCGAATATGTTCCTCAAAATCTTCCCGATAACCGACAAGCGCCTTATATACCATGTGTGCCGCCTCGTATCCGATGGCACAGTCGGAGGACTGCATAATTGACAATGCGCTGTCCTCCATTAAATCAAGTGTATCTAAAGTTGCATTTCCGTTCAATACGTCCGTGAGTAAATTCTTCAATTGCCAGAGACCAACACGGCACGGTACACATTTTCCACAAGTCTGTGCGTGACAAAGTTCAAGAAAAGCCCTTGCCATATCTACAGGACATAATCCCGGAGGGCTTGCCTCGATTCTTCGTTCCAAATCTTTGTAAAGCCCTTCTACAACAACCTGTGCTTT
>k99_4035 flag=1 multi=3.0000 len=420
```

Press q.
Now, let's see if we can get some insight into the lenghts of the contigs that have been produced.
In the following command, be sure to put '>' into double quotes (">"), or you will just remove your assemblies:

```
grep ">" SRR341725.megahit_asm/final.contigs.fa | awk '{print $4}' | cut -f2 -d"=" | sort -n | tail -n 10
```

You get:
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ grep ">" SRR341725.megahit_asm/final.contigs.fa | awk '{print $4}' | cut -f2 -d"=" | sort -n | tail -n 10
126168
126295
126807
135316
140121
147158
147953
154844
173376
217021
```

So the longest contig obtained is 217021 nt. Conversely, the 10 shortest are long:
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ grep ">" SRR341725.megahit_asm/final.contigs.fa | awk '{print $4}' | cut -f2 -d"=" | sort -n | head -n 10
200
200
200
200
200
201
201
201
201
201
```

We will now use a short custom script. It only serves to cut out all cut out all contigs shorter than 1000. 
Binning softwares will do this operation as well, but we wan't to avoid mapping reads against too-short contigs in the previous phase:

We first copy this custom script:
```
cp /data/pmanghi/assembly/filter_contigs.py .
python filter_contigs.py ${s}.megahit_asm/final.contigs.fa ${s}.megahit_asm/contigs_filtered.fasta
```

## Step n.3: Mapping metagenomic reads against contigs: almost no assembly pipeline can work without this.
Remember, this step is to obtain a function called depth of coverage: pratically, a genome chunk from a very abundant species will absorb many more reads than a lowly-abundant one. Translating, this implies that the average number of reads hitting a nucleotide in the first chunk, will be much higher than the second.

The next commands are used to perform this additional mapping step. They are commented because they take time. You can copy their final results instead:
```
## bowtie2-build ${s}.megahit_asm/contigs_filtered.fasta ${s}.megahit_asm/contigs_filtered
## bowtie2 -x ${s}.megahit_asm/contigs_filtered -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -S ${s}.sam -p 8 2> ${s}.bowtie2.log
## samtools view -bS ${s}.sam > ${s}.bam
## samtools sort ${s}.bam -o sorted_${s}.bam

cp /data/pmanghi/assembly/sorted_${s}.bam .
```

We will then use an utility named jgi_summarize_bam_contig_depths, to summarize these alignments into a file storing the average depth per contigs: 

```
jgi_summarize_bam_contig_depths --outputDepth ${s}_depth.txt sorted_${s}.bam 2> ${s}_depth.log
```

You can have a look at this file:
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ head -n 30 ${s}_depth.txt
contigName      contigLen       totalAvgDepth   sorted_SRR341725        sorted_SRR341725-var
k99_58093       1076    4.8563  4.8563  4.8463
k99_47613       1094    5.06673 5.06673 4.00426
k99_19368       1040    4.30816 4.30816 1.45755
k99_25017       1179    3.90885 3.90885 2.19563
k99_17754       1179    4.15192 4.15192 3.37584
k99_25824       1519    16.0843 16.0843 40.379
k99_32280       1667    6.45302 6.45302 7.21057
k99_12106       1029    3.89886 3.89886 1.78317
k99_39544       1172    4.70594 4.70594 6.45623
k99_21790       1018    4.56054 4.56054 5.09504
k99_20177       1029    5.0289  5.0289  4.67684
k99_8878        1329    17.1174 17.1174 20.5106
k99_18563       1060    3.856   3.856   1.55883
k99_16141       1529    5.90197 5.90197 7.38685
k99_4037        1025    3.51503 3.51503 2.01145
k99_29052       1968    4.87683 4.87683 7.49146
k99_41964       1635    4.5746  4.5746  5.33099
k99_59706       1283    58.7923 58.7923 229.212
k99_26632       1272    4.62459 4.62459 3.69628
k99_22598       1123    4.83725 4.83725 4.2908
k99_43578       2206    4.86813 4.86813 3.21384
k99_10493       1521    4.33402 4.33402 3.413
k99_17755       1045    3.28629 3.28629 1.14966
k99_13722       1269    5.43921 5.43921 7.36239
k99_37929       1340    4.80391 4.80391 5.06005
k99_14526       1441    4.68863 4.68863 4.06821
k99_2422        1019    4.62982 4.62982 4.70518
k99_63735       2652    5.2338  5.2338  3.91677
k99_46807       1068    4.37302 4.37302 1.73361
```

## Step n.4: Binning, i.e. grouping assemblies into genomes using MetaBat2
```
metabat2 -i SRR341725.megahit_asm/contigs_filtered.fasta -a ${s}_depth.txt -o ${s}_bins/bin -m 1500 --unbinned -t 8 > ${s}_metabat2.log
```

This has practically distributed the contigs into "genomes", called "bins":
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ ls SRR341725_bins
bin.1.fa   bin.12.fa  bin.15.fa  bin.18.fa  bin.20.fa  bin.23.fa  bin.26.fa  bin.29.fa  bin.31.fa  bin.34.fa  bin.4.fa  bin.7.fa  bin.BinInfo.txt     bin.tooShort.fa
bin.10.fa  bin.13.fa  bin.16.fa  bin.19.fa  bin.21.fa  bin.24.fa  bin.27.fa  bin.3.fa   bin.32.fa  bin.35.fa  bin.5.fa  bin.8.fa  bin.BinMembers.txt  bin.unbinned.fa
bin.11.fa  bin.14.fa  bin.17.fa  bin.2.fa   bin.22.fa  bin.25.fa  bin.28.fa  bin.30.fa  bin.33.fa  bin.36.fa  bin.6.fa  bin.9.fa  bin.lowDepth.fa
```

## Step n.5: Estimate MAG quality using checkM2
```
checkm2_db="/data/CheckM2_database/uniref100.KO.1.dmnd"
checkm2 testrun --database_path ${checkm2_db} --threads 8
```
We now run it: checkM2 will estimate the two key-metrics to evaluate an assembled genomes:

1) Completeness
2) Contamination
   
```
checkm2 predict -i ${s}_bins -o ${s}_checkm2 -x .fa --database_path ${checkm2_db} --threads 8
```

You should see:
```
(assembly) cdonati@mbc1:/data/cdonati/assembly$ checkm2 predict -i ${s}_bins -o ${s}_checkm2 -x .fa --database_path ${checkm2_db} --threads 8
[06/17/2026 12:37:38 PM] INFO: Running CheckM2 version 1.1.0
[06/17/2026 12:37:38 PM] INFO: Custom database path provided for predict run. Checking database at /data/CheckM2_database/uniref100.KO.1.dmnd...
[06/17/2026 12:37:41 PM] INFO: Running quality prediction workflow with 8 threads.
[06/17/2026 12:37:41 PM] WARNING: Skipping bin bin.lowDepth.fa as it has a size of 0 bytes.
[06/17/2026 12:37:41 PM] WARNING: Skipping bin bin.tooShort.fa as it has a size of 0 bytes.
[06/17/2026 12:37:41 PM] INFO: Calling genes in 37 bins with 8 threads:
    Finished processing 32 of 37 (86.49%) bins.
```

Now this command has computed the metrics: normally, for bacterial genomes you want to genomes that have:

1) Completeness>50%
2) Contamination<5%

Where those found at Completeness<90 are considered (high-quality). Depending on the kind of application you may want to select only these: for instance, NCBI accepts high-quality (completeness <90%, contamination <5%) as correctly determined genomes (sort of, as if they had had been isolated).

Let's filter checkM2 results:
```
awk -F'\t' '$2 > 50 && $3 < 5' ${s}_checkm2/quality_report.tsv > ${s}_checkm2/quality_report_filtered.tsv

mkdir -p ${s}_bins_filtered

cut -f1 ${s}_checkm2/quality_report_filtered.tsv | while read -r value; do cp ${s}_bins/${value}.fa ${s}_bins_filtered/; done
```

## APPROACH 2: NEXTFLOW PIPELINES
Visit the link:
[Nextflow pipeline metagenomic tutorials](https://github.com/claudiodonati/Elixir_Bari_2026/edit/main/README.md)

# Hands-on n.6 - Genome annotation
In this part we cover two main questions in genome annotations. For each genome we can ask two questions:

1) Which genes are in this genome?
2) Which species is this genome?

## Step n.1: Annotate the bin genetic repertoire using Bakta

Bakta has replaced Prokka as the most popular genome annotator. Bakta can be for few minutes quite memory intensive, so we'll turn off short ORFs detection to limit RAM usage.

```
mkdir -p ~/annotations
cd ~/annotations

for genome in ../assembly/SRR341725_bins_filtered/*.fa; do bakta --db /data/bakta_database/db-light --threads 4 --skip-sorf --skip-crispr --skip-trna --skip-tmrna --skip-rrna ${genome}; done
```
Check Bakta results:
```
(genome_annotation) cdonati@mbc1:/data/cdonati/annotations# cat -S bin.18.tsv
# Annotated with Bakta
# Software: v1.12.0
# Database: v6.0, light
# DOI: 10.1099/mgen.0.000685
# URL: github.com/oschwengers/bakta
#Sequence Id    Type    Start   Stop    Strand  Locus Tag       Gene    Product DbXrefs
contig_1        cds     769     1284    -       GAANHE_00001            hypothetical protein    
contig_1        cds     1271    2353    -       GAANHE_00002            HNH endonuclease        SO:0001217, UniRef:UniRef50_A0A416ESJ6
contig_2        cds     137     658     -       GAANHE_00003            FMN reductase [NAD(P)H] SO:0001217, UniRef:UniRef50_A0A380LQ39
contig_2        cds     696     1229    -       GAANHE_00004            DUF based on E rectale Gene description SO:0001217, UniRef:UniRef50_C4ZAI0
contig_2        cds     1251    1709    -       GAANHE_00005            DUF4367 domain-containing protein       SO:0001217, UniRef:UniRef50_D4MXY1
contig_2        cds     1812    2237    -       GAANHE_00006            ABC transporter permease        SO:0001217, UniRef:UniRef50_A0A1Q9JR67
contig_2        cds     2465    3847    +       GAANHE_00007            Argininosuccinate lyase SO:0001217, UniRef:UniRef50_A0A3R6QUT5
contig_2        cds     3861    5582    +       GAANHE_00008            Lipoprotein     SO:0001217, UniRef:UniRef50_A0A1Q9JR25
contig_2        cds     5987    6142    -       GAANHE_00009            hypothetical protein    
contig_2        cds     6227    6571    -       GAANHE_00010            hypothetical protein    SO:0001217, UniRef:UniRef50_A0A173RSC7
contig_2        cds     6920    7234    +       GAANHE_00011    yqjI    Transcriptional regulator YqjI  SO:0001217, UniRef:UniRef50_A0A173RSA6
contig_3        cds     112     1215    -       GAANHE_00012            Transposase     SO:0001217, UniRef:UniRef50_A0A0M6WZX4
contig_3        cds     1714    1914    -       GAANHE_00013            hypothetical protein    
contig_3        cds     1926    2345    -       GAANHE_00014            hypothetical protein    PFAM:PF01381.28, PFAM:PF12844.13
contig_4        cds     283     1512    -       GAANHE_00015            Argininosuccinate synthase      SO:0001217, UniRef:UniRef50_C4Z4C1
contig_4        cds     1813    2853    +       GAANHE_00016            N-acetyl-gamma-glutamyl-phosphate reductase     SO:0001217, UniRef:UniRef50_C4Z4C2
contig_4        cds     2880    4106    +       GAANHE_00017    argJ    Arginine biosynthesis bifunctional protein ArgJ SO:0001217, UniRef:UniRef50_F7V5T9
contig_4        cds     4175    5107    +       GAANHE_00018            Acetylglutamate kinase  SO:0001217, UniRef:UniRef50_A0A7C7V1P3
contig_4        cds     5109    6302    +       GAANHE_00019            Acetylornithine aminotransferase        SO:0001217, UniRef:UniRef50_A0A173TPG2
```

While, if you need the sequence:
```
(genome_annotation) cdonati@mbc1:/data/cdonati/annotations# less -S bin.18.faa 
```




