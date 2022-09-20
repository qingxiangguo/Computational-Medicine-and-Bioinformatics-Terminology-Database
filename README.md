# Computational-Medicine-and-Bioinformatics-Terminology-Database
A comprehensive mapping database, including detailed information, of the terminology and concepts in computational medicine and bioinformatics field.

# Contributors
Qingxiang (Allen) Guo  
Postdoctoral Fellow  
Northwestern University, Feinberg School of Medicine
Qingxiang.guo@northwestern.edu

# Introduction
In the database, I provide a list of commonly used computational medicine and bioinformatic terminology for quick reference. This database will be kept updating. Feedback or experience is warmly welcomed.

# Terminology Content
## Cancer immune Evasion through loss of MHC Class I antigen presentation
Major histocompatibility class I (MHC I) molecules bind peptides derived from a cell's expressed genes and then transport and display this antigenic information on the cell surface. This allows CD8 T cells to identify pathological cells that are synthesizing abnormal proteins, such as cancers that are expressing mutated proteins. In order for many cancers to arise and progress, they need to evolve mechanisms to avoid elimination by CD8 T cells. MHC I molecules are not essential for cell survival and therefore one mechanism by which cancers can evade immune control is by losing MHC I antigen presentation machinery (APM). Not only will this impair the ability of natural immune responses to control cancers, but also frustrate immunotherapies that work by re-invigorating anti-tumor CD8 T cells, such as checkpoint blockade. 

##Castration-resistant prostate cancer (CRPC) 
A form of advanced prostate cancer. With CRPC, the cancer no longer completely responds to treatments that lower testosterone. It shows signs of growth, like a rising PSA (prostate-specific antigen), even with low levels of testosterone.

## Checkpoint blockade immunotherapy
Patients are treated with antibodies that block negative regulatory molecules, such as PD-1/PD-L1 or CTLA4, which normally restrain T cell responses. This kind of therapy can reinvigorate a patient's anti-tumor T cell responses, which then can cause tumors to shrink and even lead to cures in some patients

## Chimeric reads 
Chimeric reads occur when one sequencing read aligns to two distinct portions of the genome with little or no overlap. This could be like sequence A mapped to 85156-85257 bp of genome, while part of sequence A mapped to 85273-85320 bp of genome. Then, sequence A is a chimeric read. Chimeric reads are indicative of structural variation. Chimeric reads are also called split reads.

## Exitron
Exitrons are exon-like introns located within protein-coding exons. Removal or retention of exitrons through alternative splicing increases proteome complexity and thus adds to phenotypic diversity.Exitrons are defined as introns within protein-coding exons that, when retained, maintain the protein-coding potential of the transcript. Marquez and colleagues argue that four features distinguish exitrons from other introns: high GC content, absence of stop codons, overrepresentation of a size corresponding to multiples of three nucleotides, and prevalence of synonymous substitutions (as usually observed for exonic sequences).

Transcripts with exitrons in their sequences can be distinguished from those with retained introns in three ways. First, transcripts containing exitrons are transported out of the nucleus to be translated, whereas those containing introns are identified as incompletely processed and are kept in the nucleus where they cannot be translated. Second, only transcripts with exitrons of lengths not divisible by three have the potential to incorporate premature termination sequences, while sequences with introns normally result in premature termination. Third, exitron transcripts are usually the major isoform, but those with introns are only present in small amounts.

## F1-score
The F1-score combines the precision and recall of a classifier into a single metric by taking their harmonic mean. It is primarily used to compare the performance of two classifiers. Suppose that classifier A has a higher recall, and classifier B has higher precision. In this case, the F1-scores for both the classifiers can be used to determine which one produces better results.

## Gene fusions
A fusion gene is a hybrid gene formed from two previously independent genes. It can occur as a result of translocation, interstitial deletion, or chromosomal inversion. Fusion genes have been found to be prevalent in all main types of human neoplasia. The identification of these fusion genes play a prominent role in being a diagnostic and prognostic marker.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/e770f449513ba588b77ff7ed02f0cff742b64fcd/imgs/Gene_Fusion_Types.png">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/f419b91f0972151a4b70fdea8155eae54ab06905/imgs/Figures-to-explain-terminology-A-Intact-exon-IE-type-and-broken-exon-BE-type.png">
</div>

## Gene fusion - Intact exon (IE) type fusion  
The original exon is retained intact after the fusion, and the original exon structure is not affected. For example, in the above figure A, Exon2 of Gene A and Exon1 of Gene B are fused and the sequences of the two exons are retained intact.

## Gene fusion - Broken exon (BE) type fusion  
Which refers to the fusion without retaining the original intact exon sequence. For example, in the above figure A, part of the sequence of Exon3 of Gene A is fused with Exon2 of Gene B. In the new gene after fusion, part of the sequence of Exon3 from Gene A is lost.

## Gene fusion - Breakpoint  
The location on the genome where the fusion of two fused genes occurs, such as the site where Gene A (blue) and Gene B (green) are fused in Figure B above.

## Gene fusion - Spanning read
a paired-end read that matches across the fusion site to two fusion genes, respectively, such as the pair of reads matching to Gene A (blue) and Gene B (green) in Figure B above.

## Gene fusion - Split read
A read that matches exactly to the fusion site, as shown in the right panel in Figure B above.

## Gene fusion - Anchor length
The length of the left and right ends of the read spanning the fusion site, as shown in the right panel in Figure B above.

## Gene fusion - Short insert size
The shorter distance between two reads in double-end sequencing paired-end sequencing, generally a few hundred bp.

## Gene fusion - Long insert size
The longer distance between two reads in double-end sequencing mate-pair sequencing, generally several kilobases or even longer.

## Germline INDEL
INDEL that is identified in germline (i.e., blood‐extracted) DNA samples. For germline indels from healthy genomes, they are mainly genetic variants with the type and position of the indels presumably conserved in sub-populations or super populations. In other words, they are less random compared with somatic variants and usually do not lead to diseases.

## Germline mutation 
A genetic change in a germ cell (egg or sperm) that becomes the DNA of each cell in the offspring's body. a variant (or mutation) contained in the germline can be passed from parent to offspring and is therefore inherited. They are also called germline mutations. INDELs are identified by removing the germline mutations from INDELs.

## Homopolymer
Homopolymer, e.g. A-A-A-A-A-A-A-A-A-A-A, simple repetition

## Human genome assembly - b37
The Broad Institute created a human genome reference file based on GRCh37. When people at The Broad Institute's Genomics Platform refer to the hg19 reference, they are actually referring to b37.

## Human genome assembly - GRCh37
GRCh37 is the Genome Reference Consortium Human genome build 37. As of May 7, 2014 it has been replaced with GRCh38 as the standard reference assembly sequence used by NCBI. Unlike other sequences, GRCh37 is not from one individual's genome sequence, but is built from reference sequences of different individuals. In essence: GRCh37 is identical to hg19 on the main contigs (chr1-24), but differ on chrM.

## Human genome assembly - GRCh38
GRCh38 is an improved representation of the human genome compared to GRCh37, where many gaps were closed, sequencing errors corrected and centromere sequences modelled. For the state-of-the-art of the human genome and its annotation, go to GRCh38.

## Human genome assembly - hg19
There are a few minor differences between GRCh37 and hg19. The contig sequences are the same but the names are different, i.e. "1" may need to be converted to "chr1". In addition UCSC hg19 is currenly using the old mitochondrial sequence but NCBI and Ensembl have transitioned to NC_012920.

## Human genome assembly - hg38
The same with GRCh38.

## Human genome assembly - hs37d5
hs37d5 (known also as b37 + decoy) was released by The 1000 Genomes Project (Phase II), which introduced additional sequence (BAC/fosmid clones, HuRef contigs, Epstein-Barr Virus genome) to the b37 reference to help reduce false positives for mapping. Note that this one uses the primary assembly of GRCh37.

## Intergenic region
An intergenic region is a stretch of DNA sequences located between genes. Intergenic regions may contain functional elements and junk DNA. Intergenic regions should not be confused with intragenic regions (or introns), which are short, non-coding regions that are found within genes, especially within the genes of eukaryotic organisms.Intergenic regions contain a number of functional DNA sequences such as promoters and regulatory elements, enhancers, spacers, and (in eukaryotes) centromeres. They may also contain origins of replication, scaffold attachment regions, and transposons and viruses.[2] Non-functional DNA elements such as pseudogenes and repetitive DNA, both of which are types of junk DNA, can also be found in intergenic regions—although they may also be located within genes in introns.

## Intron retention 
An overview of the intron retention (IR) mechanism: different isoforms can be produced from a single gene through AS. (A), Isoforms with introns fully spliced are sent out of the nucleus for translation. Intron-retaining isoforms (IRIs) can be generated through IR (no intron retention): (B), In most cases, the IRIs are degraded by the nonsense-mediated decay (NMD) pathway, the reason being that retained introns often contain premature termination codons (PTCs) that can trigger NMD (with intron retention): (C), In some cases, the IRIs are detained in the nucleus, and in response to stimuli these IRIs can undergo further splicing to remove the retained intron, before being exported out of nucleus for translation (with intron retention): (D), In the case of cytoplasmic splicing, IRIs are shuttled to the cytoplasm for preservation and may be subject to further splicing (with intron retention): (E), In yet another case, IRIs escape from the NMD pathway and are translated into protein isoforms, which, compared with normal protein isoforms, are often truncated and may lose domains; however, it could also be that the alternative protein isoforms include extra domains formed by the amino acid sequences translated from retained introns (with intron retention).  

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/89c8759153191b960c4f4c064d2967a6965c0a87/imgs/fgene-11-00586-g001.jpg">
</div>

## LNCaP cells
LNCaP cells are a cell line of human cells commonly used in the field of oncology. LNCaP cells are androgen-sensitive human prostate adenocarcinoma cells derived from the left supraclavicular lymph node metastasis from a 50-year-old caucasian male in 1977. They are adherent epithelial cells growing in aggregates and as single cells.

## Melanoma
Melanoma, the most serious type of skin cancer, develops in the cells (melanocytes) that produce melanin — the pigment that gives your skin its color.

## MHC class I
MHC class I molecules are expressed in all nucleated cells and also in platelets—in essence all cells but red blood cells. It presents epitopes to killer T cells, also called cytotoxic T lymphocytes (CTLs). A CTL expresses CD8 receptors, in addition to T-cell receptors (TCR)s. When a CTL's CD8 receptor docks to a MHC class I molecule, if the CTL's TCR fits the epitope within the MHC class I molecule, the CTL triggers the cell to undergo programmed cell death by apoptosis. Thus, MHC class I helps mediate cellular immunity, a primary means to address intracellular pathogens, such as viruses and some bacteria, including bacterial L forms, bacterial genus Mycoplasma, and bacterial genus Rickettsia. In humans, MHC class I comprises HLA-A, HLA-B, and HLA-C molecules.

## Precision and Recall
Precision is calculated by dividing the true positives by anything that was predicted as a positive. Recall (or True Positive Rate) is calculated by dividing the true positives by anything that should have been predicted as positive.

## SAM (file format)
<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a4ddccceb15fd1d1ee05ae6bb0183febb48feae4/imgs/1.png">
</div>

Here is another example of SAM header part:

@HD	VN:1.0	SO:coordinate  
@SQ	SN:1	LN:249250621	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:1b22b98cdeb4a9304cb5d48026a85128  
@SQ	SN:2	LN:243199373	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:a0d9851da00400dec1098a9255ac712e  
@SQ	SN:3	LN:198022430	AS:NCBI37	UR:file:/data/local/ref/GATK/human_g1k_v37.fasta	M5:fdfd811849cc2fadebc929bb925902e5  
@RG	ID:UM0098:1	PL:ILLUMINA	PU:HWUSI-EAS1707-615LHAAXX-L001	LB:80	DT:2010-05-05T20:00:00-0400	SM:SD37743	CN:UMCORE  
@RG	ID:UM0098:2	PL:ILLUMINA	PU:HWUSI-EAS1707-615LHAAXX-L002	LB:80	DT:2010-05-05T20:00:00-0400	SM:SD37743	CN:UMCORE  
@PG	ID:bwa	VN:0.5.4  
@PG	ID:GATK TableRecalibration	VN:1.0.3471	CL:Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate, TileCovariate], default_read_group=null, default_platform=null, force_read_group=null, force_platform=null, solid_recal_mode=SET_Q_ZERO, window_size_nqs=5, homopolymer_nback=7, exception_if_no_tile=false, ignore_nocall_colorspace=false, pQ=5, maxQ=40, smoothing=1  

Here is another example of SAM alignment part, we separte each line for clear:

1:497:R:-272+13M17D24M	113	1	497	37	37M	15	100338662	0	CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAG	0;==-==9;>>>>>=>>>>>>>>>>>=>>>>>>>>>>	XT:A:U	NM:i:0	SM:i:37	AM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37  

19:20389:F:275+18M2D19M	99	1	17644	0	37M	=	17919	314	TATGACTGCTAATAATACCTACACATGTTAGAACCAT	>>>>>>>>>>>>>>>>>>>><<>>><<>>4::>>:<9	RG:Z:UM0098:1	XT:A:R	NM:i:0	SM:i: AM:i:0  X0:i:4	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:37  

19:20389:F:275+18M2D19M	147	1	17919	0	18M2D19M	=	17644	-314	GTAGTACCAACTGTAAGTCCTTATCTTCATACTTTGT	;44999;499<8<8<<<8<<><<<<><7<;<<<>><<	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:4	X1:i:0	XM:i:0	XO:i:1	XG:i:2	MD:Z:18^CA19  

9:21597+10M2I25M:R:-209	83	1	21678	0	8M2I27M	=	21469	-244	CACCACATCACATATACCAAGCCTGGCTGTGTCTTCT	<;9<<5><<<<><<<>><<><>><9>><>>>9>>><>	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:5	X1:i:0	XM:i:0	XO:i:1	XG:i:2	MD:Z:35  
  
The SAM is divided into two parts, the header section and the alignment section.

The comment information is not the focus of the SAM file, but is a record of the process of generating and processing the SAM file, and is specified to start with @, with different tags to indicate different information, mainly include:

@HD, a description of the standard-compliant version, the order in which the sequences are compared.  
@SQ, description of the reference sequence.  
@RG, description of the sequence (read) on which the comparison was made.  
@PG, description of the program used.  
@CO, arbitrary description information.  

Use command: samtools view -H <.bam> to see the header

In the comparison section we focus on the first 11 columns  

<b>The first column</b>: QNAME - Query template NAME is very straightforward, it is the name number of the fragment of the query, which corresponds to the name of the reads in the sequencing fastq data. For example, the above HWI-.... .2168

<b>The second column</b>: FLAG - bitwise FLAG tag that characterizes the result of the comparison in binary terms. If one or more of the following occur simultaneously, then the sum of the numbers in the corresponding first column is the number you see in the second column of the sam file. For example, if 83=64+16+2+1, it means that the reads are matched to the reference genome, both sequences are matched at the double end, and the current fragment is reversed to the reference genome.

1: means this sequence is sequenced by PE double-end sequencing  
2: means the sequence matches the reference sequence exactly, no insertion or deletion  
4: means the sequence is not mapped to the reference sequence  
8: means the other end of this sequence is not mapped to the reference sequence, for example, this sequence is R1, its corresponding R2 end sequence is not mapped to the reference sequence  
16: represents the sequence is compared to the negative chain of the reference sequence  
32: represents the sequence corresponding to the other end of the sequence is compared to the negative chain of the reference sequence  
64 : means this sequence is the R1 end sequence, read1  
128 : the sequence is the R2 end sequence, read2.  
256 : represents this sequence is not the primary comparison, a sequence may be compared to more than one position of the reference sequence, only one is the primary comparison position, the others are all secondary  
512: means this sequence failed in QC and was filtered out  
1024: means this sequence is a PCR duplicate  
2048: means this sequence has a part of chimeric PCR primers  

For example, if 99, 99=1+2+32+64, which means this sequence is R1 end sequence, which means this sequence corresponds to the other end sequence compared to the negative strand of the reference sequence, which means this sequence and the reference sequence match exactly, no insertion deletion, it is double end sequencing

<b>The third column</b>: RNAME-- Reference sequence NAME of the reference genome corresponding to the name of the chromosome or contig, scaffold. If * it means there is no result on the comparison. If @SQ header lines are present, RNAME must be present in one of the SQ-SN tag.

<b>The fourth column</b>: POS-- 1-based leftmost mapping POSition fragment leftmost mapping to the reference genome, counting from 1. If it is 0, it means there is no match.

<b>The fifth column</b>: MAPQ-- MAPping Quality mapping score, -10 log10 Pr{mapping position is wrong}- logarithmic transformation of the probability that the mapping position is wrong with a base of ten. According to this formula, we can also know that the larger the value, the lower the probability that the pairing is wrong.

<b>The sixth column</b>:  CIGAR - CIGAR string stands for Concise Idiosyncratic Gapped Alignment Report. First of all, there are the following letters that represent different meanings. Most of the meanings are easy to understand, such as sequence matching, insertions and deletions, and skipping the fragment. The document mentions that Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.  
<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/491ce3bf5a3dc72277520ada9c45236498586c4b/imgs/2.png">
</div>
 
For example, 251M, means the alignment of the full match 

30M3D126M3D58M37S, it means that the reads are 30bp match + 3bp missing + 126bp match + 3bp missing + 58bp match + 27bp soft clipping.  

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/491ce3bf5a3dc72277520ada9c45236498586c4b/imgs/3.png">
</div>

If there is a sequence that comes from a mature mRNA, if this sequence skips exactly 200bp of intron in the middle and 75bp mapped to exon before and after, what should the CIGAR value of this sequence be written?CIGAR:75M200N75M

<b>The seventh column</b>: read2's chromosome name on the reference sequence, if not available use "*", same as "="

<b>The eighth column</b>:: PNEXT: position of read2 on the reference sequence

<b>The ninth column</b>:: TLEN: length of inserted fragment  

TLEN - observed Template LENgth If all read segments are mapped to the corresponding reference sequence, the absolute value of TLEN is equal to the distance (end-start+1) between the mapped end of the template sequence and the mapped start of the template sequence (including both ends). It should be noted that the bases on the comparison do not include sof-clipped bases. If the read segment is compared to the start of the leftmost segment of the template, the TLEN field is positive, and if the comparison is to the start of the rightmost segment, it is actually an antisense strand and the TLEN field is negative. If the starting position of the comparison is the same at both ends, then any positive or negative number is assigned. If there is only a single chain, the value is 0. And the positive and negative numbers of any intermediate segments are undefined.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/c44f8e5480793e2a8932699303f4eb57640f47ca/imgs/4.png">
</div

Note that the length of the inserted fragment, is the entire length below, not between R1 and R2.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/c44f8e5480793e2a8932699303f4eb57640f47ca/imgs/5.png">
</div

Because the values given by sam are all at the 5 end, which is not the same as the actual in vivo DNA. So you need to add the later value, add the length of the sequence, and then subtract the starting value of the previous sequence.

<b>The 10th column</b>: SEQ: segment SEQuence. This field can be a * when the sequence is not stored. If not a * the length of the sequence must equal the sum of lengths of M/I/S/=/X operations in CIGAR. An = denotes the base is identical to the reference base. No assumptions can be made on the letter cases.

<b>The 11th column</b>: QUAL: ASCII of base QUALity plus 33 (same as the quality string in the Sanger FASTQ format). A base quality is the phred-scaled base error probability which equals −10 log10 Pr{base is wrong}. This field can be a * when quality is not stored. If not a ‘*’, SEQ must not be a ‘*’ and the length of the quality string ought to equal the length of SEQ.

<b>The tag column</b>: At the end of the SAM, you will see something like: XT:A:U NM:i:0 SM:i:37 AM:i:0 X0:i:1 X1:i:0 XM:i:0 XO:i:0 XG:i:0 MD:Z:37. TAGs are optional fields on a SAM/BAM Alignment. A TAG is comprised of a two character TAG key, they type of the value, and the value. A user can also use any additional tags to store any information they want. TAGs starting with X, Y, or Z are reserved to be user defined.

Example:  

XT:A:U  - user defined tag called XT.  It holds a character.  The value associated with this tag is 'U'.  
NM:i:2  - predefined tag NM means: Edit distance to the reference (number of changes necessary to make this equal the reference, excluding clipping)  
XS:A:+, XS:A- -XS to indicate the genomic strand that produced the RNA from which the read was sequenced. When your sequencing is unstranded, mappers can still add a strand to a read if it crosses a splice site, and that splice site has a canonical splice site sequence - it can use this sequence to work out which way the RNA is going. All spliced reads would contain the XS tag.  
SA: BWA uses SA tag for marking chimeric reads. 

## SA tag in SAM file
After aligning with bwa mem, chimeric reads will have an SA tag. Their format is: SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)
Each element in the list represents a part of the chimeric alignment. Conventionally, at a supplementary line, the first element points to the primary line. Strand is either ‘+’ or ‘-’, indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based coordinate.

Example:

HWI-ST387:139:C03WJABXX:5:2108:15315:193815 16 scf7180000067989 85156 60 60M41S * 0 0 TTGAAGTCAAGAAAGTGGTAAAGAGAGATTAATAGGGGTATCTCAGCTACAACAAATATTATATTAAATTAAATGGTTAATCTTGCTTTGCTCACCATAAA * NM:i:2 MD:Z:31G1C26 AS:i:50 XS:i:0 <b>SA:Z:scf7180000067989,85273,-,54S47M,60,1</b>;

HWI-ST387:139:C03WJABXX:5:2108:15315:193815 272 scf7180000067989 85273 60 54H47M * 0 0 AATATTATATTAAATTAAATGGTTAATCTTGCTTTGCTCACCATAAA * NM:i:1 MD:Z:11T35 AS:i:42 XS:i:22 <b>SA:Z:scf7180000067989,85156,-,60M41S,60,2</b>;

scf7180000067989 is the reference name, 85156 is the position, - is the strand, 60M41S is the cigar, 60 is the map quality, 2 is the NM. Each SA tag will store information of another alignments in a chimeric alignment, other itself.

The software TransIndel (Yang et al., 2018, BMC genomics) used the SA tag of BWA-MEM to detect Indel. The algorithm is as follows:
Sometimes it is hard to distinguish splicing and deletion, for a true deletion the length of left mapped + lengths of right mapped must be equal to the length of the read itself.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/bc9f0ec4b7ee9ac2b0b419c2d1204f596aec288a/imgs/transindel_Algorithm.png">
</div

.  
## Slippage during polymerase chain reaction amplification  
Slippage during PCR, also known as replication slippage, is a form of mutation that causes trinucleotide or dinucleotide amplification and sometimes even contraction during DNA replication, resulting in tandem repeat sequences.

## SNV (Single nucleotide variant)
A SNV can be rare in one population but common in a different population. Sometimes SNVs are known as single nucleotide polymorphisms (SNPs), although SNV and SNPs are not interchangeable. To qualify as a SNP, the variant must be present in at least 1% of the population.

## Soft clipping and hard clipping
About soft clipping and hard clipping, it means that when query matching, some sequences are not matched completely, but soft clipping will keep the unmatched part afterwards, and hard clipping can remove the matched part completely. For example, when there is only part of the fragment is compared to the reference sequence when comparing. The difference is that a Soft Clip will eventually retain the corresponding sequence in the sequence that follows, while a Hard Clip will delete the fragment directly in the sequence that follows. BWA-MEM uses soft clipping CIGAR operation for supplementary alignments. By default, BWA-MEM uses soft clipping for the primary alignment and hard clipping for supplementary alignments.

Hard Clip exists with the intention of reducing the redundancy of BAM file sequences, for example, there is a read which can be compared to two places A, B. In place A, it is 60M90S, and in place B it is 60H90M, at this time a read actually already has the complete sequence information in position A, and the information in position B is actually redundant. So a marker form like Hard Clip can be introduced at location B, and it will be able to mark the sequence at location B as secondary.

## Splicing junctions
Key to defining the complexity of alternative splicing within a gene is the identification of splice junctions (SJs), which occur at exon-exon boundaries and are typically characterized in pairs representing both the donor site (5’ intron boundary to 3’ upstream exon boundary) and acceptor site (3’ intron boundary to 5’ downstream exon boundary).
  
## Variant allele frequency (VAF)
The fraction of genome copies in the (tumor or control) sample affected by the variant, is either 0.0, 0.5, or 1.0 for germline variants, reflecting absence, heterozygosity, or homozygosity, respectively. In contrast, allele frequencies of somatic variants vary across the whole range from 0.0 to 1.0, depending on the clonal structure of the tumor sample and its impurity (the ratio of healthy genome copies in the tumor sample).
  
The variant allele frequency (VAF; also known as variant allele fraction) is used to infer whether a variant comes from somatic cells or inherited from parents when a matched normal sample is not provided. A variant is potentially a germline mutation if the VAF is approximately 50% or 100%.

## XS tag in SAM file (only for Tophat and HISAT, not BWA-MEM, BWA-MEM also has the XS tag, but with other meaning)
Where your sequencing is unstranded, mappers can still add a strand to a read if it crosses a splice site, and that splice site has a canonical splice site sequence - it can use this sequence to work out which way the RNA is going. Mappers such as TopHat (which cufflinks was originally designed to work with) only accepted splice sites with one of the canonical splice sites (i'm guessing), thus all spliced reads would contain the XS tag. This is definately true for HiSat.









