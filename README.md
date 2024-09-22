## NGS-Pipeline
Next generation sequencing workflow pipline, this documentation outlines the steps involved in NGS analysis workflow from data acquisition to final visualization and interpretation step by step.

## Next generation sequencing with definitions and file formats
NGS it is a technology for determining the sequence of DNA or RNA to study the genetic variation associated with diseases, it is a sequencing method where millions of sequencing reactions are carried out in parallel, produces massive number of sequences from a single cell in a short period of time.

## Before we begin, it is better to recoginze some definitions and termenologies will be uesd in the following lines.
(High-throughput sequencing) : A key feature of NGS, refering to the ability for processing a massive number of DNA sequencing simultaneously with lower cost unlike the earlier methods.

(Read) : short DNA or RNA sequencies are generated during sequencing.

(Library) : DNA fragments collection have been prepared for analysis.

(Run) : the entire sequencing process untill the end.

(flowcell) : physical chip where the DNA is loaded for sequencer.

(Adapter) : synthetic short fragments of DNA that have a "T" base 'thymine' overhang. are added to the DNA fragments during library preparation for sequencing, it is unique as it attaches the fragments. sometimes these adapters still be present in the raw sequencing data.

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

(BAM files) : binary alignment map binary version of bam file is desgined for more efficient storage, retrieval and processing of aligned sequence data.
for more information about sam and bam files format (https://samtools.github.io/hts-specs/SAMv1.pdf)



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
Data used in this workflow from E.coli, evolved and ancestral raw data. two ancestral samples and four evolved samples(evolved from two different strains)
you can download raw dataset from here (https://genomics.sschmeier.com/downloads)


## Scope and Objective
You will learn how to analyse Next generation sequencing data in this pipeline, which will focus on detecting the genomic variations of E.coli dataset to understand the functional impact on biological processes and adaptive traits compared to evolved strains.


## 'SRAtool' Installation
if you want to perform this pipeline on external data, you can install "SRAtool" for downloading dataset from NCBI.

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
    
    #download FASTQ file
    fasterq-dump --verbose "accession number for dataset"

    #decompress fastq data for further steps
    gzip -d data.fastq.gz

## Fastp Installation by Conda
Conda is an open-source package management system tool for bioinformatics, to be more informed look here(https://docs.conda.io/en/latest/), but substantially, conda has channels for storing and downloading tools and packages, so after installing conda, we will install conda channels if needed. for more reusability, organization and version control make an isolated environment for each package/tool to prevent any errors or conflicts with the base of conda environment (it is optional but highly recommended to create a separate environment). have a look at "fastp" usage (https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp)
Alternatively, you can just download "fastp" manualy from here(https://anaconda.org/bioconda/fastp), but it is recommended to install conda as we will set up more tools and packages in further analysis.

    #download miniconda latset version
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    #Run installer script
    bash Miniconda3-latest-Linux-x86_64.sh

    #open bashrc file via nano(or any text editor) and put package’s path at the end of the file, then save the changes
    nano ~/.bashrc

    #Add Path
    export PATH="~/miniconda3/bin:$PATH"

    #Reload shell configuration
    source ~/.bashrc

    #Ensure that conda is installed 
    conda --version

    #After restart your terminal, run this into your terminal to show the channels, it will show channels that are set up via the latest version of conda, if you want to add more channles, run the second command.
    conda config --show channels
    conda config --<CHANNEL NAME>

    #create new environmet
    conda create -n fastp-env

    #initialize conda in your shell and restart your terminal
    conda init

    #env activation, to navigate to fastp-env 
    conda activate fastp-env

    #install fastp
    conda install -c bioconda fastp

    #check if "fastp" is installed by list env packages, you will see the name, version and channel(where tools stored) of every installed tool.
    conda list


## Quality control and Data trimming  
In this step we will filter out data from any low-quality bases, adapters, overrepresented sequences, biases, for improving the accuracy of read mapping and overall data quality.
The most two common tool for trimmig are "fastp and trimmomatic" in our case we will use "fastp" because of:
* we do not have information about the adapters used in dataset used in this pipeline and "fastp" automaticaly detects and removes adapters based on the data without knowing the used adapters.
* "fastp" faster and more effecient, user-friendly with automatic optimization and fewer parameters and doesn’t require detailed configurations.

      #make new directory 
      mkdir data
    
      #make new directory for trimmed data
      mkdir trimmedData  
    
      #navigate to the directory where the data are located
      cd data
    
      #trimming data with fastp(ancestral)
      fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 2 --html trimmedData/anc.fastp.html --json trimmedData/anc.fastp.json -i anc_R1.fastq.gz -I anc_R2.fastq.gz -o        trimmedData/anc_R1.fastq.gz -O trimmedData/anc_R2.fastq.gz
    
      #trimming data with fastp(evolved sample1)
      fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 2 --html trimmedData/evol1.fastp.html --json trimmedData/evol1.fastp.json -i evol1_R1.fastq.gz -I                     evol1_R2.fastq.gz -o trimmedData/evol1_R1.fastq.gz -O trimmedData/evol1_R2.fastq.gz
    
      #trimming data with fastp(evolved sample2)
      fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread 2 --html trimmedData/evol2.fastp.html --json trimmedData/evol2.fastp.json -i evol2_R1.fastq.gz -I                     evol2_R2.fastq.gz -o trimmedData/evol2_R1.fastq.gz -O trimmedData/evol2_R2.fastq.gz
  
* --detect_adapter_for_pe: by default the auto-detection for adapter is for single-end data only. so in this case we need to turn on paired-end option but put that in mind specifies PE optionmay slow the performance and increase time of fastp since it will do additional analysis to identify and remove adapter. good to know it results a cleaner output as it is improve the accuracy of the overlap analysis and overall quality for data.
* overrepresentation_analysis: for identifying potentail PCR duplication.
* --correction: enabling this option as by default it is disabled in "fastp" so it perform overlap analysis only for paired-end data if a proper overlap is found, it can correct mismatched basepaires in overlapped regions.
* --cut_right: trim the nucleotide from the right-end of reads(3’end) based on the quality scores cutting until the quality value of the last base is above threshold.
* --thread 2: use 2 threads for parallel processing.
* --html: generating an html report showing results of trimming
* --json: json files for record-keeping and furthur analysis
* trimmedData/: the directory where the results will be saved
  
FastQC
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

    #create new directory for fastqc results
    mkdir fastqc-trimmed

    #navigate to the directory where the data are trimmed
    cd trimmedData

    #run the following command in all fastq files, "-o" specifies the output directory where FastQC will store its results, "*.fastq.gz" this is a wildcard that matches all fastq files in the navigated directory.
    fastqc -o fastqc-trimmed/ *.fastq.gz

## FastQC result 
will generate a '.html' file which contains a detailed report about the sequence data and .zip file, so how we determine if our data need to be filtered or trimmed? if we open the .html file, we will focus on graphs and plots that have warning sign. For more details about FastQC report look here (https://mugenomicscore.missouri.edu/PDF/FastQC_Manual.pdf)

## MultiQC
As we are working on many samples at once, we will run "multiqc", which is a tool used to aggregate and summarize results from versious quality control samples instead of reviewing the results from each samples one by one, multiQC compiles all these results into one easy-to-read report.

    #download multiQC
    conda install -c bioconda multiqc

    #navigate to the directory where fastqc results are stored
    cd fastqc-trimmed

    #run multiqc to combine all fastqc reports
    multiqc .
    
    


![fastqc_sequence_counts_plot](https://github.com/user-attachments/assets/925a04fc-d42f-4995-875f-60fc4662588b)

![fastqc_per_sequence_quality_scores_plot](https://github.com/user-attachments/assets/c7728f87-6e18-4a8b-91ba-dfb1c5494187)


![fastqc_per_sequence_gc_content_plot](https://github.com/user-attachments/assets/1f5919c7-f6ff-4e1d-8794-53efc9625220)
![fastqc_per_base_n_content_plot](https://github.com/user-attachments/assets/1034ac8c-f912-454c-92d1-c18077cb56f6)
![fastqc_sequence_duplication_levels_plot](https://github.com/user-attachments/assets/4ffbc9d1-f5d3-41b4-9e55-5df7e8eba79a)

## Genome Assembly
Sequencing technologies generates short fragments of DNA called 'reads', genome assembly process is for reconstructing the complete sequence of genome smaller fragments of DNA, the aim is to arrange these fragments in the correct order and orientation to produce a representation of the organism’s genome.(https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/genome-assembly)
we will create reference genome from ancestral samples, becuase of:
* less contaminated: ancestral samples are less likely to be contaminated, which can lead to errors in genome assembly
* preserved genetic diversity
* identification of ancient variants.
* resolve ambiguities in assembly.
  
look here for more info(https://www.nature.com/articles/s41559-022-01956-z#:~:text=Reconstructed%20ancestral%20genomes%20are%20similar,in%20silico%20reconstructions%20when%20available.) (https://research.pasteur.fr/fr/publication/reconstruction-of-hundreds-of-reference-ancestral-genomes-across-the-eukaryotic-kingdom/)

## SPAdes and Quast
'SPAdes' it is powerful genome assembler and It's frequently used for bacterial genomes due to its ability to handle small genomes efficiently.(https://ablab.github.io/spades/)
'Quast' assess the quality of genome assemblies. we will run quast on the two scaffolds.fasta files to compare the results.
Note: Spades runs multiple k-mer sizes by default since it generates second/another fasta file depends on the length of the reads.

    #downloading SPAdes, Quast and creating new environment
    conda create -n assembly spades quast

    #activate conda environment
    conda activate assembly

    #create new directory for assemblies
    mkdir assembly

    #navigate to the trimmed data directory
    cd trimmedData

    #run spades
    spades.py -1 anc_R1.fastq.gz -2 anc_R2.fastq.gz -o assembly

    #run Quast '-o' specifies te directory where quast’s results will be saved, and specifies the directory for two scaffolds.fasta files
    quast -o data/assembly assembly/scaffolds.fasta assembly/K77/scaffolds.fasta 

    
an overview of key assembly statistics for two scaffolds fasta files.
![Screenshot (179)](https://github.com/user-attachments/assets/1b6a15eb-e646-4b65-b08d-0c12c7b98ffd)
![Screenshot (182)](https://github.com/user-attachments/assets/2b79b2ff-3a22-46f8-8504-f443a23d23bf)
![Screenshot (181)](https://github.com/user-attachments/assets/ea806291-ab41-45de-b0f8-a46cdccc7fa3)
![Screenshot (180)](https://github.com/user-attachments/assets/6a547989-076c-41f7-8ef3-1d93d68877f4)
![Screenshot (184)](https://github.com/user-attachments/assets/b9c5020e-773c-41ea-ac40-eae31fa11790)
![Screenshot (183)](https://github.com/user-attachments/assets/a5c1b1e6-875a-4e38-9573-1f70ee6f319d)

By conclusion, assembly-scaffolds and K77-scaffolds represents to be very similar across all metrics: contig lengths, GC content, and coverage depth as well the key statistics as the dataset is high quality, The assembly methods are robust and reliable but this similarity doesn’t gaurantee our data is perfect and the report results is significant there could be subtle differences specially we was creating the referance genome from ancestral paired-end raw data, so we will do further analysis such as(alignment, variant calling, annotation) compared to evolved samples to see the variations.  
For more details about Quast report (https://quast.sourceforge.net/docs/manual.html)


## Alignment 
we will align the reads(evolved data) against reference genome(created from ancestral data), before mapping step we need to index our reference genome as its improve the efficiency, speed, and accuracy of the read mapping process without indexing, the read mapper would need to scan the entire reference genome for every read, which would be computationally expensive and time-consuming.
will install samtools for interacting and analysing bam files, bwa the aligner used in mapping the reads against the reference genome, and qualimap used for quality control of mapped sequencing data.


    #installations and create new environment
    conda create --yes -n mapping samtools bwa qualimap r-base

    #activate the environment
    conda activate mapping
    
    #navigate to trimmedData directory
    cd trimmedData
    
    #make new directory for mappings results
    mkdir mappings
    
    #indexing
    bwa index scaffolds.fasta

    #mapping the evolved samples to reference genome(we just created) with bwa
    bwa mem assemblyRef/scaffolds.fasta evol1_R1.fastq.gz evol1_R2.fastq.gz > mapping/evol1.sam 
    bwa mem assemblyRef/scaffolds.fasta evol2_R1.fastq.gz evol2_R2.fastq.gz > mapping/evol2.sam 

## BAM formatting and fixing SAM flags
After alignment/mapping step with paired-end sometimes SAM flags (bitwise values used to describe the characteristics of read and it’s alignment) can be wrong, in this case we use 'SAMtools' to fix this problem, to make sure paired-end reads relationships are correctly represented. and then compress SAM files to BAM to reduce their size for faster and smoother analysis. so we will performe three post-alignment steps:

* use 'fixmate' to fix any errors and compress to BAM
* sort by name

        #sort by name. 'sort -n': sort SAM by read name, it is required for 'fixmate' for ensuring paired-end reads are next to each other
        #'-O' specifies output format as sam, '|' this pipe for passing the sorted sam file from the first command part to the second part without intermediate file
        # 'samtools fixmate' for cleaning up any errors, '-m' ensure additional fields are included to the output such as('MC' mate coordinates, 'TLEN' template length, '-O bam' specifies the output file as bam
        #'-' tells 'fixmate' to take the input file from the previous command "evol1.sam', 'rm evol1.sam' remove uneeded file
        #Evol1
        samtools sort -n -O sam evol1.sam | samtools fixmate -m -O bam - evol1.fixmate.bam 


        #remove evol1.sam file as it is not needed
        rm evol1.sam
  
        #sort bam file in order
        samtools sort -O bam -o evol1.sorted.bam evol1.fixmate.bam

        #remove fixmate file to save space
        rm evol1.fixmate.bam


       #Evol2
          
     


## Remove Duplicates
During library preparation step in NGS, polymerase chain reaction(PCR) is used to amplify DNA, this preduce multiple copies of the same DNA fragment(PCR duplicates) this is not a real biological data that may lead to biases in furthur analysis, so we are going to remove it

    #remove duplicates, '-r' to remove dupicates, '-S' handle supplementary alignments 
    samtools markdup -r -S evol1.sorted.bam evol1.sorted.dedup.bam 


    #remove the original data to safe space
    rm evol1.sorted.bam  

    #statistics overview
    samtools flagstat evol1.sorted.dedup.bam

    
![Screenshot (185)](https://github.com/user-attachments/assets/dd7ab278-27e7-4ad4-a815-6d97d2a4c94c)

These statistics shows that:
* This BAM file contains 1,619,698 reads, of which 97.73% are mapped to the reference genome.
* The data is paired-end, with most pairs (96.31%) being properly paired.
* A small percentage (0.88%) of reads are singletons, and a small number of reads have mates mapped to a different chromosome.
* there are no duplicates as expected after removing them above by 'samtools markdup'                                                    

According to this statistics, there are some mates mapped to different chromosomes or contigs also known as(discordnat read pairs) can occur for several reasons:
* Structural variations ( translocations, inversions,duplications)
* chimeric DNA fragments: During library  preparation, DNA fragments from different chromosomes may accidentally ligated togethere, leading to chimeric reads.
* Sequencing or Assembly Errors
* Gene Duplications and paralogous regions: in regions where two or more genes are similar, reads can map to different contigs because of the aligner may split them between the most similar locations
* discordnat read pairs can be usefull in SVs detection, cancer genomics, genome assembly improvement, Evolutionary studies.

    

      
  
    
          
    




    
    




