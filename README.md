# NGS-Pipeline
Next generation sequencing workflow pipline, this documentation outlines the steps involved in NGS analysis workflow from data acquisition to final visualization and interpretation step by step.

## Next generation sequencing with definitions and file formats
NGS it is a technology for determining the sequence of DNA or RNA to study the genetic variation associated with diseases, it is a sequencing method where millions of sequencing reactions are carried out in parallel, produces massive number of sequences from a single cell in a short period of time.

## Before we begin, it is better to recoginze some definitions and termenologies will be uesd in the following lines.
(High-throughput sequencing) : A key feature of NGS, refering to the ability for processing a massive number of DNA sequencing simultaneously with lower cost unlike the earlier methods.

(Read) : short DNA or RNA sequencies are generated during sequencing.

(Library) : DNA fragments collection have been prepared for analysis.

(Run) : the entire sequencing process untill the end.

(flowcell) : physical chip where the DNA is loaded for sequencer.

(Adapter) : short fragments of DNA that have a "T" base 'thymine' overhang.

(Alignment) : the proccess of comparing and matching the sequenced reads(short DNA fragments) to the referance genome.

(Tile) : specific and defined area on the flowcell or the sequenced chip used during the sequencing process.

(Single cell sequencing) : powerful  technology that analysis the genetic material from individual cells for high-resolution insights.

(Basecalling) : is the fundamental process of translating the raw data from sequencing machine and converting it into actual DNA sequence.

(Single-end sequencing) : the DNA or RNA fragments are sequenced from only one end, meaning the sequencer reads the fragment from one end.

(paired-end sequencing) : the DNA or RNA fragments are sequenced from two ends, meaning the sequencer reads the the fragment from both ends.

## File formats
Real-time analysis in sequencing technologies stores the individual base call data in intermediate files called 'BCL' these files are filtered, demultiplexed and converted into a sequence file format called 'FASTQ' when the sequencer completes the run.

(FASTQ files) : Human-readable uses simple text lines, consists of a number of records, each record having four lines of data
    -first line is the sequencing header strats with an '@'
    -second line is the sequence
    -third line is placeholder line begins with '+' separating the sequence from it’s scores.
    -fourth line is quality score line.

(FASTA files) : is text-based format for repressenting the DNA, RNA, and protein sequences, consists of two lines
    -first line starts with a '>' character followed by a unique identifier and description of the sequence
    -second line contians the sequence.

(SAM files) : sequence alignment map is the most basic human readable text-based format,it is commonly used for storing the biological sequence aligned to the reference genome.

(BAM files) : binary alignment map binaru version of bam file is desgined for more efficient storage, retrieval and processing of aligned sequence data.  



## How NGS works?
NGS involves in four main steps:
1. Sample collection and preparation.
2. Library Amplification.
3. Sequencing.
4. Basecalling.

**Sample collection and preparation**:
the process starts with the extraction of DNA or RNA form organism of your interet, to determine the collection method you should determine what the sample is.
fragmening the DNA or RNA into smaller pieces to load these fragments in the flowcell, adding adapter sequences to each fragment for binging DNA/RNA to the sequencung flowcell.
During the DNA fragmentation the ends of DNA pieces can be ragged 'unblunt' this means that DNA fragment is not ready to bind with other DNA fragments or adapters, to fix this issue 'en-repair process' is used, in which uneven ends are smooths out. After ends are blunt-ended 'A-tailing' is performed to add an 'A' base 'adenine' is added to the 3' ends of DNA fragments to ensure that the DNA fragments do not bind to each other and attach easily to the adapters as the 'A' base in DNA fragment will match up with the 'T' base on the adapter.

**Library Amplification**:
increase the quantity of DNA using polymerase chain reaction(PCR).

**Sequencing**:
loading the DNA fragments onto a flowcell, then the fragments are amplified to create clusters of identical DNA molecules.

**Basecalling**:
base calling step is performing directly after the sequencer machine output it’s result (raw data) and converting raw data signals 'durnig sequencing the machine detects signals such as flashes of lights or changes in electrical current' into a readable DNA sequence "A,T,C,G". during base calling step phred score (Quality score) is performed. it is pre-base estimates error emitted by the sequencer to express the level of confidence and measure the quality of each nucleotide base call in DNA sequence it is logaritmic scale of the base call error probability. the output of the base calling process is resulted in sequence of nucleotides along with their corresponding quality score stored in FASTQ files contains base sequence and it’s score for each base.
    

## NGS workflow





![NGS workflow](https://github.com/user-attachments/assets/3c69d0ed-204f-4575-a881-17032546096b)





## Source of the data
Data used in this workflow divided into two types, evolved and ancestral Ecoli raw data. we will use two ancestral sample and four evolved samples form two E.coli strains
you can download raw dataset from here (https://genomics.sschmeier.com/downloads)


## Scope and Objective
You will learn how to analyse Next generation sequencing data in this pipeline, which will focus on detecting the genomic variations of E.coli dataset to understand the functional impact on biological processes and adaptive traits compared to evolved strains.


## Installation
if you want to perform this pipeline on external data, you can install SRAtool for downloading dataset from NCBI.

    #SRA-toolkit downloading
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.1.1-ubuntu64.tar.gz

    #Extract toolkit
    tar -xvzf sratoolkit.3.1.1-ubuntu64.tar.gz

    #add to PATH
    echo 'export PATH=$PATH:~/Downloads/sratoolkit.3.0.0-ubuntu64/bin' >> ~/.bashrc

    #Reload the shell configuration
    source ~/.bashrc

Downloading FASTQ file to be well-organized, create the directory "fastqs" and download the file, run this in your terminal;

    #make new directory 
    mkdir fastqs

    #navigate to fastqs directory
    cd fastqs

    #make new directory for single-end files
    mkdir singel-end

    #download FASTQ file
    fasterq-dump --verbose SRR030834

    #in the following steps we will perform Quality control, so it is recommended to decompress FASTQ file if it is downloaded in compressed form, you can do it by running this command
    gzip -d SRR030834.fastq.gz


Post-Sequencing quality

Post-sequence quality check must be performed to assess the read quality to ensure that there are no back-ground noises, the data looks good and there is no biases leads to inacurrate results. FASTQC is the     most common programm use to check the quality of the sequence it is generate summary and simple graphical reports that show an overall idea about the quality of the raw data

   
   ### downloading fastqc tool
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.0.zip

    #decompress the zipped file with unzip
    fastqc_v0.12.0.zip

    #navigate to FaSTQC directory   
    cd FastQC

    #to make "fastqc" executable
    chmod 755 fastqc

    #topen ".bashrc" file (or any text editor) to add fastqc to the path
    nano ~/.bashrc

    #add the path at the bottom of the file and replace this "/home/mypath/FastQC" with the actual path where the fastqc file is located then save the changes and exit from text editor
    export PATH=$PATH:/home/mypath/FastQC

    #apply the changes
    source ~/.bashrc

    #for testing the fastqc, if the path set correctly, it will show help instructions
    fasqc -h

    #quality control
    fastqc SRR030834.fastq

## FastQC result 
will generate a '.html' file which contains a detailed report about the sequence data and .zip file, so how we determine if our data need to be filtered or trimmed? if we open the .html file, we will focus on graphs and plots that have warning sign, which is;

* Per base sequence quality: this ia show box plots of the quality distributions on each position across all bases of the reads, the background graph is divided into three regions representing the quality score, the green region (Q>20) are verry good, the orange region (20<Q<28) are acceptable, the red region (Q<20) are poor, in our case, we have sequence read’s score less than 20 

* Per sequence quality score: heatmap shows the quality score across different tiles of the sequencing chip with coloring code indicating the QS the colors ranges from green(high quality) to red(lowquality), if there are low quality score thats because of physical defetcs, issues with the sequencing chemistry or problem with sample preparation, so it should be filtered out.

 * per base sequence content: provide detailed view of the nucleotides compositions at each base position and this identify if there any biases in base position,   there are four colored lines representing the persentage of each base " A in green, T in red, G in black, and C in blue " it will be a good read if the lines are  relatively flat to each other, in our case the per base sequence content have warning sing that should be filtered.

* Per sequence GC content: we have two curves, blue curve that represent the normal distribution of the GC content in the sequence reads, and the purple curve represents the up-normal distribution of the GC content.

* Sequence duplication levels: high sequence duplication level indicating issues with PCR amplification biase or adapter containation.

* Adapter content: identify the extend of adapter contamination of the sequence data. After making FastQC report, we notice that our data need to be filtered or trimming step.

**Trimmig with fastp**

  

    #run fastqc on the trimmed result file, which will result HTML report
    fastqc SRR030834_trimmed.fastq

    #run multiqc
    multiqc .

