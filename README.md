# Elixir_ShotgunMetagenomics_2026

- [Hands-on n.1 - Preprocessing of standard metagenomic data](#user-content-hands-on-n1---preprocessing-of-standard-metagenomic-data)
- [Hands-on n.2 - Taxonomic profiling using marker genes with MetaPhlAn 4](#user-content-hands-on-n2---taxonomic-profiling-using-marker-genes-with-metaphlan-4)
- [Hands-on n.3 - Taxonomic profiling using k-mers: Kraken + Bracken taxonomic profiling](#user-content-hands-on-n3---taxonomic-profiling-using-k-mers-kraken--bracken-taxonomic-profiling)
- [Hands-On 4: functional profiling at the community level using HUMAnN 4](#user-content-hands-on-4-functional-profiling-at-the-community-level-using-humann-4)
- [Hands-on n.5 - Metagenome assembly and binning](#user-content-hands-on-n5---metagenome-assembly-and-binning)
  * [Approach 1: following a protocol step-by-step](#user-content-approach-1-following-a-protocol-step-by-step)
  * [Approach 2: trying Nextflow pipelines](#user-content-approach-2-trying-nextflow-pipelines)
- [BONUS: Hands-on n.6 - Taxonomic profiling beyond the level of species using StrainPhlAn](#user-content-bonus-hands-on-n6---taxonomic-profiling-beyond-the-level-of-species-using-strainphlan)

# Hands-on n.1 - Preprocessing of standard metagenomic data

#### Step n.0: log in into your machine and explore the configuration
```
ssh YOUR-NAME@212.189.202.106
```

Check our yotu current location
```
pwd
```
did it return /data/YOURNAME/?

#### Step n.1: check if your environment is present
```
which python
conda info --envs
```
you should see:

```
# conda environments:
#
base                  *  /data/anaconda3
ai-microbiome            /data/anaconda3/envs/ai-microbiome
bowtie2                  /data/anaconda3/envs/bowtie2
checkm2                  /data/anaconda3/envs/checkm2
humann4                  /data/anaconda3/envs/humann4
kraken_+_bracken         /data/anaconda3/envs/kraken_+_bracken
megahit                  /data/anaconda3/envs/megahit
metabat2                 /data/anaconda3/envs/metabat2
mpa                      /data/anaconda3/envs/mpa
samtools                 /data/anaconda3/envs/samtools
school_notebooks         /data/anaconda3/envs/school_notebooks
trimmomatic              /data/anaconda3/envs/trimmomatic
workflows                /data/anaconda3/envs/workflows
```
did it work ?

#### Step n.2: raw data pre-processing on fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata

```
## conda create -n <trimmomatic> -c bioconda trimmomatic ## DON'T DO IT. WE DID ALREADY
## conda create -n <bowtie2> -c bioconda bowtie2 ## DON'T DO IT. WE DID ALREADY
## conda create -n <samtools> -c bioconda samtools ## DON'T DO IT. WE DID ALREADY

mkdir 1_pre-processing
cd 1_pre-processing

wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
unzip input.zip
cd input
```

#### Step n.3: Define variable "s" with the sampleID and run TRIMMOMATIC
```
s="seq"

conda activate trimmomatic

trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}1.fastq ${s}2.fastq \
${s}_filtered_1.fastq ${s}_unpaired_1.fastq ${s}_filtered_2.fastq ${s}_unpaired_2.fastq \
ILLUMINACLIP:/data/anaconda3/envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75

for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done
```

#### Step n. 4: Generate bowtie2 index of the human genome GCF_009914755.1_T2T-CHM13v2.0.fna (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0.fna)

Activate conda, after setting the path of your version of the human genome
```
human_gen_path="/data/human_genome/"
conda deactivate
conda activate bowtie2
```

Run bowtie alignment against the human genome:
```
##VERSION 4 HOURS LONG:
## mkdir -p human_genome/
## bowtie2-build ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0.fna human_genome/GCF_009914755.1_T2T-CHM13v2.0 ### DON'T RUN IT! IT TAKES A FEW HOURS TO BE EXECUTED

bowtie2 -x ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0 -1 ${s}_filtered_1.fastq -2 ${s}_filtered_2.fastq -S ${s}.sam --very-sensitive-local -p 8

conda deactivate
conda activate samtools

samtools view -bS ${s}.sam > ${s}.bam
samtools view -b -f 12 -F 256 ${s}.bam > ${s}.bothunmapped.bam
samtools sort -n -m 5G -@ 2 ${s}.bothunmapped.bam -o ${s}.bothunmapped.sorted.bam
samtools fastq ${s}.bothunmapped.sorted.bam -1 >(gzip > ${s}_filtered.final_1.fastq.gz) -2 >(gzip > ${s}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n
```

then remove the intermediate files, and check files size 
```
rm ${s}.sam; rm ${s}.bam; rm ${s}.bothunmapped.bam; rm ${s}.bothunmapped.sorted.bam ### REMOVE THE INTERMEDIATE FILES

for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done; echo; for i in *.gz; do echo -ne "${i}\t"; zcat "$i" | wc -l; done
```
Did the preprocessing produce the same exact number of reads in R1 and R2 ?

# Hands-on n.2 - Taxonomic profiling using marker genes with MetaPhlAn 4
#### Step n.1: Setup correct variables, activate environment and navigate to the right folders

We create the conda environment **we did it already**
```
conda deactivate
## conda create -n <mpa> -c conda-forge -c bioconda python=3.11 metaphlan=4.2.0
conda activate mpa
```

We move to use it
```
cd ~

mkdir 2_metaphlan
cd 2_metaphlan

## metaphlan --install --db_dir metaphlan_databases --idx mpa_vJan21_CHOCOPhlAnSGB_202103 ## DON'T RUN THIS
```

#### Step n.2: download metagenomic samples
```
mpa_db="/data/metaphlan_databases/"

## db_version="mpa_vJan25_CHOCOPhlAnSGB_202503"
db_version="mpa_vJan21_CHOCOPhlAnSGB_202103"

wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz

```

#### Step n.3: Run MetaPhlAn 4

Take look at the MetaPhlAn parameters
```
metaphlan -h
```

Then run it
```
s="SRS014476-Supragingival_plaque"
metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --nproc 8 --db_dir ${mpa_db} --index ${db_version}

for s in SRS014459-Stool SRS014472-Buccal_mucosa SRS014470-Tongue_dorsum SRS014494-Posterior_fornix SRS014476-Supragingival_plaque; do
    cp /data/course_backup/2_metaphlan/${s}_profile.txt ${s}_profile.txt; done

## s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
##    --nproc 8 --db_dir ${mpa_db} --index ${db_version}
## s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
##    --nproc 8 --db_dir ${mpa_db} --index ${db_version}
## s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
##    --nproc 8 --db_dir ${mpa_db} --index ${db_version}
## s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
##    --nproc 8 --db_dir ${mpa_db} --index ${db_version}

merge_metaphlan_tables.py *_profile.txt | grep -P "clade_name|UNCLASSIFIED|t__" > metaphlan_table.tsv
```

# Hands-on n.3 - Taxonomic profiling using k-mers: Kraken + Bracken taxonomic profiling
#### Step n.1: Check everything is set up and create kraken + bracken DB

```
conda deactivate
## conda create -n <kraken_+_bracken> -c bioconda kraken2
conda activate kraken_+_bracken

cd ~
mkdir 3_kraken
cd 3_kraken

## wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
## mkdir -p kraken_DB && tar -xvzf k2_standard_08gb_20250402.tar.gz -C kraken_DB
## git clone https://github.com/jenniferlu717/KrakenTools.git
## chmod +x KrakenTools/*
```

#### Step n.3: Let's have a look at Kraken parameters
```
kraken2 -h
```

Run Kraken
```
for s in SRS014459-Stool.fasta.gz SRS014470-Tongue_dorsum.fasta.gz SRS014472-Buccal_mucosa.fasta.gz SRS014476-Supragingival_plaque.fasta.gz SRS014494-Posterior_fornix.fasta.gz;

do kraken2 --db /data/kraken_DB/ --threads 8 --report `basename ${s%.fasta.gz}`.kraken2_report.txt --output `basename ${s%.fasta.gz}`.kraken2_output.txt ../2_metaphlan/${s}; done
```

Run Bracken
```
for s in SRS014459-Stool.fasta.gz SRS014470-Tongue_dorsum.fasta.gz SRS014472-Buccal_mucosa.fasta.gz SRS014476-Supragingival_plaque.fasta.gz SRS014494-Posterior_fornix.fasta.gz;

do bracken -d /data/kraken_DB/ -i `basename ${s%.fasta.gz}`.kraken2_report.txt -o `basename ${s%.fasta.gz}`.bracken_abundance.txt -w `basename ${s%.fasta.gz}`.bracken_report.txt -l S -t 150; done

for s in *.bracken_report.txt; do /data/KrakenTools/kreport2mpa.py --display-header -r ${s} -o ${s%.txt}.mpa.tsv; done

/data/KrakenTools/combine_mpa.py -i *.bracken_report.mpa.tsv -o merged_bracken_table.tsv

sed 's/.bracken_report.txt//g' merged_bracken_table.tsv | grep -P 'Classification|s__' | sed 's/Bacillati/Bacteria/g' | sed 's/Pseudomonadati/Bacteria/g' > bracken_table.tsv
```

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
