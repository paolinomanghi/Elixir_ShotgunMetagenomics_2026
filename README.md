# Elixir_ShotgunMetagenomics_2026

- [Hands-on n.1 - Preprocessing of standard metagenomic data](#user-content-hands-on-n1---preprocessing-of-standard-metagenomic-data)
- [Hands-on n.2 - Taxonomic profiling using marker genes with MetaPhlAn 4](#user-content-hands-on-n2---taxonomic-profiling-using-marker-genes-with-metaphlan-4)
- [Hands-on n.3 - Taxonomic profiling using k-mers: Kraken + Bracken taxonomic profiling](#user-content-hands-on-n3---taxonomic-profiling-using-k-mers-kraken--bracken-taxonomic-profiling)
- [Hands-on n.4 - Functional profiling at the community level using HUMAnN 4](#user-content-hands-on-4-functional-profiling-at-the-community-level-using-humann-4)
- [Hands-on n.5 - Metagenome assembly and binning](#user-content-hands-on-n5---metagenome-assembly-and-binning)
  * [Approach 1: Following a protocol step-by-step](#user-content-approach-1-following-a-protocol-step-by-step)
  * [Approach 2: Trying Nextflow pipelines](#user-content-approach-2-trying-nextflow-pipelines)
- [BONUS: Hands-on n.6 - Taxonomic profiling beyond the level of species using StrainPhlAn](#user-content-bonus-hands-on-n6---taxonomic-profiling-beyond-the-level-of-species-using-strainphlan)

# Hands-on n.1 - Preprocessing of standard metagenomic data
## Step n.0: log in into your machine and explore the configuration
```
ssh YOUR-NAME@212.189.202.106
```

Press Enter
Check our yotu current location
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
In order to activate the right Anaconda settig, we must just set a series of environmental variables. Type:
```
source /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/.conda
```

Now, see the (base) string in the bottom-left corner of the screen ?
Type:

```
which conda
which python
```

Did it return:
```
/mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/bin/python 
/mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/bin/conda
```
?

## Step n.1: check if your environments are all set up and ready
```
which python
conda info --envs
```
You should see:

```
# conda environments:
#
base                 * /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda
assembly               /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/assembly
genome_annotation      /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/genome_annotation
humann4                /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/humann4
preprocessing          /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/preprocessing
tax_profiling          /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/tax_profiling
```

Did it work ?
In case it doesn't work is very easy to install all Conda environments.
0) Explain how to search for correct program versions:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
 
1) create the environment for Data Preprocessing:
```
conda create -n preprocessing -c bioconda -c conda-forge trimmomatic bowtie2 samtools
```
2) create the environment for the Assembly (the genome reconstruction):
```
conda create -n assembly -c bioconda -c conda-forge megahit bowtie2 metabat2 checkm2 samtools
```
3) Create one environment for MetaPhlAn and for Kraken/Bracken (taxonomic profiling)
```
conda create -n tax_profiling -c bioconda -c conda-forge metaphlan kraken bracken

mkdir -p ~/database/
cd ~/database/
metaphlan --install --db_dir metaphlan_databases --idx mpa_vJan21_CHOCOPhlAnSGB_202103

wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
mkdir -p kraken_databases && tar -xvzf k2_standard_08gb_20250402.tar.gz -C kraken_databases

cd ../src
git clone https://github.com/jenniferlu717/KrakenTools.git
chmod +x KrakenTools/*

```
4) Create one environment for HUMAnN4 (functional profiling of communities)
```
conda create -n humann4 -c bioconda python=3.12
conda activate humann4
conda config --add channels biobakery
conda install humann=4.0.0a1 -c biobakery -c bioconda -c conda-forge
conda install metaphlan=4.1 -c bioconda
#### metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202403
```
5) Create one environment for Genomic Annotations
```
conda create -n genome_annotation -c bioconda bakta
```

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

As you can see, the programs inside each environment are protected, meaning that they are visible only in that environment to not interefe with other installations.

## Step n.1: raw data pre-processing on fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata
In this step, we just download a toy fastq sample. Normally, this step may take a few weeks!

```
mkdir preprocessing
cd preprocessing

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
ls /mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/preprocessing/share/trimmomatic/adapters/
```
See multiple fasta files with different type of adapters? This reflect and adapts to different type of sequencing procedure applied, allowing you to use Trimmomatic in
different experimental settings.

Now set the correct adapter sequences for this project:
```
truseq_adap="/mnt/nfs2/manghip/projects/elixir_2026/src/elixir_conda/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
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
human_gen_bowtie_indexes='/mnt/nfs2/databases/bowtie2_indexes/GCF_009914755.1_T2T-CHM13v2.0_genomic'
```
Now we'll just align the fastq against the humann h
Run bowtie alignment against the human genome:
```
##VERSION 4 HOURS LONG:
## bowtie2-build ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0.fna human_genome/GCF_009914755.1_T2T-CHM13v2.0 ### DON'T RUN IT! IT TAKES A FEW HOURS TO BE EXECUTED

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
mkdir ~/metaphlan_taxonomic_profiling
cd ~/metaphlan_taxonomic_profiling

conda activate tax_profiling

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
mpa_db="../database/metaphlan_databases/"
db_version="mpa_vJan21_CHOCOPhlAnSGB_202103"

s="SRS014476-Supragingival_plaque"

metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --nproc 8 --db_dir ${mpa_db} --index ${db_version}

#### for s in SRS014459-Stool SRS014472-Buccal_mucosa SRS014470-Tongue_dorsum SRS014494-Posterior_fornix SRS014476-Supragingival_plaque; do
####     cp /data/course_backup/2_metaphlan/${s}_profile.txt ${s}_profile.txt; done

s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
   --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
   --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
  --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
  --nproc 8 --db_dir ${mpa_db} --index ${db_version}

merge_metaphlan_tables.py *_profile.txt | grep -P "clade_name|UNCLASSIFIED|t__" > metaphlan_table.tsv
```

## **Step n.4: Perform Kraken + Bracken taxonomic profiling (of the same samples!)**

Kraken is normally coupled with Bracken. The two softwares are considered companions and have been designed by the same developers.
Although considered among the top-performers, metagenomic people are aware of a potential kraken limitation: mapping onto the whole
genomic sequence, k-mers may be ambigously located on two highly-similar species. To solve this problem, Kraken assigns reads that
may ambigously map against species A and species B their Lowest Common Ancestor (LCA). As species A and B are similar, their LCA will
likely be their genus: Kraken thus display a bias, in that it assigns reads to superior taxonomic level more easily.

Bracken is, in words, the piece of software designed to address this issue: Bracken uses Baysian statistics to redistribute ambigous reads, and as such the trick is
hydden in the correct estimation of the priors. Each species has two key-numbers:

1) the number of unique k-mers (which is a  charcteristic of the database)
2) the number of uniquely assigned reads.

Clearly, the ** random ** distribution of reads between the two species is not derived from having the same number of reads among the two, but having the same proportion of
Number of reads / Number of unique K-mers

Or, in other words, species that have very little unique k-mers in their genomes, are expected to show a greater abundance if we assumed that the number of uniquely assigned reads
is the same. In short, Bracken redistribute ambiguous reads based the probability of observing the observed number of non-ambiguous calls, considering the number of unique k-mers for
the species.

Albeit pratically this methodology works, and is the de-facto standard for non-humann communities, precautions must be taken when extensively running Bracken:

1) The resulting number of species with standard parameters is notoriously too-high (based on consideration on human being, in which we know precisely what to expect)
2) In environmental samples, the Bracken rate of False Positives may be perceived as sensitivity

When running Bracken, it is therefore practically advised to apply an additional filter at > 1000 reads, meaning lowering to zero everything that, after the Baysian statistical
algorithm, maintain a count of reads below 1000. Despite this potential pitfall, Bracken remains the standard for profiling environmental or non-common communities (for example marine metagenomes or yet-to-be-surveyed animals).

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
kraken2 --db ../database/kraken_databases/ --threads 8 --report `basename ${s%.fasta.gz}`.kraken2_report.txt --output `basename ${s%.fasta.gz}`.kraken2_output.txt ../metaphlan_taxonomic_profiling//${s};
done
```

As said, we are now going to correct the ambiguity bias in Kraken using Bracken:
```
for s in SRS014459-Stool.fasta.gz SRS014470-Tongue_dorsum.fasta.gz SRS014472-Buccal_mucosa.fasta.gz SRS014476-Supragingival_plaque.fasta.gz SRS014494-Posterior_fornix.fasta.gz; do
bracken -d ../database/kraken_databases/ -i `basename ${s%.fasta.gz}`.kraken2_report.txt -o `basename ${s%.fasta.gz}`.bracken_abundance.txt -w `basename ${s%.fasta.gz}`.bracken_report.txt -l S -t 150;
done
```

Similarly to what done for MetaPhlAn, we will now merge all the Bracken outputs into a MetaPhlAn-like taxonomy, to make them comparable. We will make use
of the suite KrakenTools, which we have previously downloaded from github.

```
for s in *.bracken_report.txt; do ../src/KrakenTools/kreport2mpa.py --display-header -r ${s} -o ${s%.txt}.mpa.tsv; done
```

This last command has created files in a format that we can handle very easily. Type:
```
cat SRS014470-Tongue_dorsum.bracken_report.mpa.tsv
```
See ?

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
../src/KrakenTools/combine_mpa.py -i *.bracken_report.mpa.tsv -o merged_bracken_table.tsv
sed 's/.bracken_report.txt//g' merged_bracken_table.tsv | grep -P 'Classification|s__' | sed 's/Bacillati/Bacteria/g' | sed 's/Pseudomonadati/Bacteria/g' > bracken_table.tsv
```

# Hands-on n.3 - Scripting in R to compare the results

## Exercises:
First, set the correct conda environment:
```
conda deactivate 
conda activate r_notebook
```

Next, open jupiter notebook:
```
jupyter notebook /data/Jupyter/Alpha_diversity.ipynb --allow-root --no-browser --port=8888 --ip=127.0.0.1
```

Next, open up a new terminal and type:
```
ssh -N -L 8888:localhost:8888 <YOURNAME>@212.189.202.106
```

Note it holds but doesn't do anything. It means that a tunnel is open
Then, copy the URL that the server prompt is suggesting you into YOUR BROWSER:

Es.:
```
http://212.189.202.106:8888/tree?token=398cde02036d5c0c4e8162b5e21758c5d7f9fa90dc4eabad
```
press ctrl + c, type y and exit jupyter

# Hands-On 4: functional profiling at the community level using HUMAnN 4
#### Step n.1: Get into the right directory & install download the necessary files
```
cd ~
## conda create -n <humann4> -c bioconda python=3.12 ## DON'T DO IT. WE DID ALREADY

conda deactivate 
conda activate humann4

## conda config --add channels biobakery
## conda install humann=4.0 -c biobakery  ## DON'T DO IT. WE DID ALREADY
## conda install metaphlan=4.1 -c bioconda ## DON'T
## metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202403

mkdir -p 4_humann
cd 4_humann

## humann_databases --download chocophlan full humann_databases ## DON'T 
## humann_databases --download uniref uniref90_ec_filtered_diamond humann_databases ## DON'T
## humann_databases --download utility_mapping full humann_databases ## DON'T
```

#### Step n.2: test that HUMAnN runs properly and have a look at the HUMAnN parameters
```
humann_test
## humann_config
humann -h
```

#### Step n.3: get a sample from EBI
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/096/SRR15408396/SRR15408396.fastq.gz
```

#### Step n.4: RUN humann
```
s="SRR15408396"

metaphlan_params="--index mpa_vOct22_CHOCOPhlAnSGB_202403 -t rel_ab_w_read_stats"

## NOW YOU CAN RUN: ## --nucleotide-database /data/humann_databases/chocophlan/
## \humann --input ${s}.fastq.gz --output ${s} --threads 8  --count-normalization RPKs --metaphlan-options "${metaphlan_params}"
## rm -r ${s}/${s}_humann_temp/

## BUT IT TAKES THREE HOURS... OR YOU CAN RUN:
mkdir -p ${s}

cp /data/course_backup/4_humann/${s}/${s}_2_genefamilies.tsv ${s}/${s}_genefamilies.tsv
cp /data/course_backup/4_humann/${s}/${s}_4_pathabundance.tsv ${s}/${s}_pathabundance.tsv
```

#### Step n.5: Regrouping genes to other functional categories
```
humann_regroup_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_pfam.tsv --groups uniref90_pfam
```

#### Step n.6: Run HUMAnN on a second sample
```
s="SRR15408398"

## SAME:
## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/098/SRR15408398/SRR15408398.fastq.gz

## \humann --input ${s}.fastq.gz --output ${s} --threads 8 --count-normalization RPKs --metaphlan-options "${metaphlan_params}"
## rm -r ${s}/${s}_humann_temp/

## FOR NOW, RUN:
mkdir ${s}
cp /data/course_backup/4_humann/${s}/${s}_4_pathabundance.tsv ${s}/${s}_pathabundance.tsv
```

#### Step n.7: Merge together community profiles under different ontologies
```
mkdir -p merged

cp SRR15408396/SRR15408396_pathabundance.tsv merged/
cp SRR15408398/SRR15408398_pathabundance.tsv merged/

humann_join_tables -i merged -o merged_pathabundance.tsv --file_name pathabundance
```

Visualize the merged table:
```
less -S merged_pathabundance.tsv 
```
Press q

# Hands-on n.5 - Metagenome assembly and binning
## Approach 1: following a protocol step-by-step

#### Step n.1: check everything is set up, download a sample, and run Megahit
```
cd ~

## conda create -n <megahit> -c bioconda megahit ## DON'T DO IT. WE DID ALREADY
conda deactivate
conda activate megahit

mkdir 5_assembly
cd 5_assembly

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz
```

Let's have a look at the megahit parameters !
```
megahit -h
```

Let's start with the assembly
```
s="SRR341725"
## MEGAHIT WILL TAKE A FEW HOURS:
## megahit -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -o ${s}.megahit_asm -t 8

## FOR NOW WE CAN COPY THE RESULTS FROM MEGAHIT
mkdir -p ${s}.megahit_asm/

cp /data/course_backup/5_assembly/${s}.megahit_asm/final.contigs.fa  ${s}.megahit_asm/
cp /data/course_backup/5_assembly/${s}.megahit_asm/contigs.fasta  ${s}.megahit_asm/

## WE ALSO NEED TWO CUSTOM SCRIPT:
cp /data/course_backup/5_assembly/filter_contigs.py .
cp /data/course_backup/5_assembly/megahit2spades.py .

python megahit2spades.py ${s}.megahit_asm/final.contigs.fa ${s}.megahit_asm/contigs.fasta
python filter_contigs.py ${s}.megahit_asm/contigs.fasta ${s}.megahit_asm/contigs_filtered.fasta
```

#### Step n.2: Binning, i.e. grouping assemblies into genomes using MetaBat2
```
cd ~

mkdir 6_MAG-reconstruction
cd 6_MAG-reconstruction

conda deactivate 
## conda create -n <metabat2> -c bioconda metabat2 ## DON'T DO IT. WE DID ALREADY
conda activate metabat2

## conda install -c bioconda <bowtie2> ## DON'T DO IT. WE DID ALREADY
## conda install -c bioconda <samtools> ## DON'T DO IT. WE DID ALREADY

s="SRR341725"

cp ../5_assembly/${s}.megahit_asm/contigs_filtered.fasta ./
cp ../5_assembly/${s}_1.fastq.gz ./
cp ../5_assembly/${s}_2.fastq.gz ./
 
## bowtie2-build contigs_filtered.fasta contigs_filtered
## bowtie2 -x contigs_filtered -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -S ${s}.sam -p 8 2> ${s}.bowtie2.log
## samtools view -bS ${s}.sam > ${s}.bam
## samtools sort ${s}.bam -o sorted_${s}.bam

## COPY THE RESULT FOR NOW:
cp /data/course_backup/6_MAG-reconstruction/sorted_${s}.bam .

jgi_summarize_bam_contig_depths --outputDepth ${s}_depth.txt sorted_${s}.bam 2> ${s}_depth.log
```

#### Step n.3: Run MetaBat 2 for binning
```
metabat2 -i contigs_filtered.fasta -a ${s}_depth.txt -o ${s}_bins/bin -m 1500 --unbinned -t 8 > ${s}_metabat2.log
```

#### Step n.4: Estimate MAG quality using checkM2
```
conda deactivate
## conda create -n <checkm2> -c bioconda checkm2 ## DON'T DO IT. WE DID ALREADY

conda activate checkm2
## pip install absl-py==1.1.0 ## DON'T DO IT. WE DID ALREADY

## LET'S NOT DOWNLOAD THE DATABASE
## checkm2 database --download --path ./

## WE CAN USE A COPY
checkm2_db="/data/CheckM2_database/uniref100.KO.1.dmnd"
checkm2 testrun --database_path ${checkm2_db} --threads 8

checkm2 predict -i ${s}_bins -o ${s}_checkm2 -x .fa --database_path ${checkm2_db} --threads 8

awk -F'\t' '$2 > 50 && $3 < 5' ${s}_checkm2/quality_report.tsv > ${s}_checkm2/quality_report_filtered.tsv

mkdir -p ${s}_bins_filtered
cut -f1 ${s}_checkm2/quality_report_filtered.tsv | while read -r value; do cp ${s}_bins/${value}.fa ${s}_bins_filtered/; done
```
