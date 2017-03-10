IDP-APA (version 0.1)
Time-stamp: <2017-02-20 Yunhao Wang, Email: yunhaowang@126.com>

* Introduction
==============

Alternative cleavage and polyadenylation (APA), a common phenomenon in eukaryotes, is emerging as an important layer of gene regulation in eukaryotes and plays important regulatory roles in human development and diseases. 

Previous methods based on Second-Generation Sequencing (SGS) technology are difficult to match each identified polyA site within a gene to its derived isoform because of the limitation of sequencing length of SGS.

To achieve the isoform-resolved APA analysis, we present a tool called IDP-APA that constructs truly-expressed isoforms and identifies polyA sites by adequately integrating the respective strengths of Third-Generation Sequencing (e.g. PacBio and ONT) long reads (LRs) and SGS (e.g. Illumina) short reads (SRs). 

IDP-APA takes the aligned SRs and LRs and a known gene annotation libray as inputs, and outputs the constructed isoforms with assigned polyA sites.


* Prerequisite
==============

- Linux system

- python 2.7


* Install and Run IDP-APA
=========================

(1) Download the package (e.g. `wget https://github.com/yunhaowang/IDP-APA/archive/v0.1.1.tar.gz`) to a directory of your choice. (e.g. "/home/")

(2) Unpack it using the command `tar zxvf /home/v0.1.1.tar.gz`

(3) Now, you can run IDP-APA by the executable file `/home/IDP-APA-0.1.1/bin/idpapa`. Optionally, you can add IDP-APA into your [PATH environment variable] so that you can run `idpapa` from the command line without having to specify the entire path. For example, you can add one line `export PATH=/home/IDP-APA-0.1.1/bin:$PATH` to your "~/.bashrc".


* Inputs of IDP-APA 
===================

(1) aligned SRs (sam format).

- Currently, IDP-APA (version 0.1) only support the paired-end strand-specific RNA-seq data (e.g. prepared using Illumina TruSeq Stranded mRNA/total RNA Library Prep Kit).

- The suggested aligner is Hisat2 (verison 2.0.0-beta was used in our study, see https://ccb.jhu.edu/software/hisat2/index.shtml) with the parameter `-k 1 --rna-strandness RF --no-mixed --no-discordant`. For aligned SR pairs, the alignment flag (second column in sam file) should be 83 (read paired, read mapped in proper pair, read reverse strand, first in pair) vs 163 (read paired, read mapped in proper pair, mate reverse strand, second in pair); or 99 (read paired, read mapped in proper pair, mate reverse strand, first in pair) vs 147 (read paired, read mapped in proper pair, read reverse strand, second in pair).

(2) aligned LRs (sam format).

- The suggested aligner is GMAP (version 2016-06-09 was used in our study, see http://research-pub.gene.com/gmap/) with the parameter `-n 0 -f samse --split-output`. We suggest to use both uniquely-aligned LRs (sam format file suffixed by '.uniq' in GMAP output file) and multiply-aligned LRs (sam format file suffixed by  '.mult' in GMAP output file). 

- If using PacBio sequencing, we suggest to run PacBIo Iso-Seq pipeline to get ROI (reads of insert) and full-length non-chimera and non-full-length non-chimera LRs (the command file we used is `ConsensusTools.sh CircularConsensus --minFullPasses 0 --minPredictedAccuracy 70` and `pbtranscript.py classify --detect_chimera_nfl --min_seq_len 100`). 

(3) a known gene annotation library (gtf format).

- Considering the RefSeq annotation libray contains some genes/isoforms which have more than one copy in different genomic loci but use same gene/isoform ID, we do not suggest to use it. We suggest to use Ensembl or GENCODE. GENCODE version 25 was used in our study.

(4) an optional file (csv file) if using PacBio sequencing data.

- The output primer information file (suffixed by '.primer_info.csv', see "./example/SIRV_ROI.primer_info.csv"). This file is output when running PacBio Iso-Seq pipeline (getting full length reads, https://github.com/PacificBiosciences/cDNA_primer/wiki/RS_IsoSeq-%28v2.3%29-Tutorial-%231.-Getting-full-length-reads#commandline). In our study, it was produced by the command line version `pbtranscript.py classify`. This file is used to determinate if the polyA tail is detected for each PacBio long read (ROI, reads of insert).


* Outputs of IDP-APA
====================

(1) a modified gpd format file. 

- column 1: gene ID. The ID prefixed by 'novel_loci_sgt_' means that it is a novel singleton gene that has no overlap with any known gene loci; and the ID prefixed by 'novel_loci_mlt_' means that it is a novel multi-exon gene that has no the overlap with any known gene loci.

- column 2: isoform ID. The prefix 'novel_sgt_iso_' means that it is a novel singleton isoform; and the prefix 'novel_mlt_iso_' means that it is a novel multi-exon isoform.

- column 3: chromosome

- column 4: strand

- column 5: TSS (transcription start site, for '+' strand)

- column 6: TES (transcription end site, for '+' strand)

- column 7: CDS (coding sequence) start site (for '+' strand)

- column 8: CDS end site (for '+' strand)

- column 9: number of exons

- column 10: start positions of exons (for '+' strand)

- column 11: end positions of exons (for '+' strand)

- column 12: polyA sites identified by long reads, the delimiter is comma ','. Format is 'polyA site position' + '_' + 'number of long reads of different types'. If primer information csv file is used when running IDP-APA, 'F' means full-length non-chimera LRs with polyA tail, 'N' means non-full-length non-chimera LRs with polyA tail, and 'P' means non-full-length non-chimera LRs without polyA tail. If primer csv file is not used when running IDP-APA, 'L' means all types of LRs.

- column 13: polyA sites identified by short reads, the delimiter is comma ','. Format is 'polyA site position' + '_' + 'number of short reads of different types'. 'S' means that this polyA site is uniquely assigned to one isoform, 'M' means that this polyA site is multiply assigned to multiple isoforms.

- column 14: polyA sites identified by short reads. Two parts is splitted by semicolon ';'. The left part means polyA sites that pass the cutoff; and the right part means polyA sites that do not pass the cutoff. The cutoff is determinated by the parameter `--sc`, `--lc` and `--ds`.

- column 15: polyA type. 'NA' means no available polyA site for this isoform, 'PA' means only one available polyA site for this isoform, and 'APA' means multiple available polyA sites for this isoform.

(2) a statistical file (suffixed by ".stat") for all constructed isoforms with identified polyA sites


* Usage and Example
===================

You can use the test data in "./example/" to test IDP-APA.

`idpapa` `-s` ./example/input/SIRV_Illumina_short_reads.sam `-l` ./example/input/SIRV_PacBio_long_reads.sam `-a` ./example/input/SIRV_gene_annotation.gtf `-p` ./example/input/SIRV_PacBio_ROI.primer_info.csv `-o` ./example/SIRV_constructed_isoform_with_polyA_site.gpd `-t` ./example/temp
