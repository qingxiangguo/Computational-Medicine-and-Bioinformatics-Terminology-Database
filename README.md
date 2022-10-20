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

## Copy number variation (CNV)
Copy number variation (CNV) is a phenomenon in which sections of the genome are repeated and the number of repeats in the genome varies between individuals. Such regions may or may not contain a gene(s).

## CRAM (file format)
Compressed Reference-oriented Alignment Map (CRAM) is a compressed columnar file format for storing biological sequences aligned to a reference sequence. CRAM was designed to be an efficient reference-based alternative to the Sequence Alignment Map (SAM) and Binary Alignment Map (BAM) file formats.

## Enhancer
An enhancer is a short (50–1500 bp) region of DNA that can be bound by proteins (activators) to increase the likelihood that transcription of a particular gene will occur. They can be located up to 1 Mbp (1,000,000 bp) away from the gene, upstream or downstream from the start site. Enhancers are found mostly in the intergenic and intronic regions, while a few enhancers have been found within exons.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a73102a47ec95cf8ca7fbb3f7e938ee279da01dc/imgs/enhancer.png">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a73102a47ec95cf8ca7fbb3f7e938ee279da01dc/imgs/enhancer2.png">
</div>

Here is an enhancer diagram. Within this DNA sequence, protein(s) known as transcription factor(s) bind to the enhancer and increases the activity of the promoter. 1. DNA 2. Enhancer 3. Promoter 4. Gene 5. Transcription Activator Protein 6. Mediator Protein 7. RNA Polymerase

## Epitope
It is capable of stimulating an immune response. This is usually one to six monosaccharides or five to eight amino acid residues on the surface of the antigen. Each antigen typically has many epitopes. 

## Exitron
Exitrons are exon-like introns located within protein-coding exons. Removal or retention of exitrons through alternative splicing increases proteome complexity and thus adds to phenotypic diversity.Exitrons are defined as introns within protein-coding exons that, when retained, maintain the protein-coding potential of the transcript. Marquez and colleagues argue that four features distinguish exitrons from other introns: high GC content, absence of stop codons, overrepresentation of a size corresponding to multiples of three nucleotides, and prevalence of synonymous substitutions (as usually observed for exonic sequences).

Transcripts with exitrons in their sequences can be distinguished from those with retained introns in three ways. First, transcripts containing exitrons are transported out of the nucleus to be translated, whereas those containing introns are identified as incompletely processed and are kept in the nucleus where they cannot be translated. Second, only transcripts with exitrons of lengths not divisible by three have the potential to incorporate premature termination sequences, while sequences with introns normally result in premature termination. Third, exitron transcripts are usually the major isoform, but those with introns are only present in small amounts.

## Extrachromosomal circular DNA （ecDNAs）
Extrachromosomal circular DNA (eccDNA) is a type of double-stranded circular DNA structure that is derived and free from chromosomes. It has a strong heterogeneity in sequence, length, and origin and has been identified in both normal and cancer cells. In contrast to previously identified circular DNA structures (e.g., bacterial plasmids, mitochondrial DNA, circular bacterial chromosomes, or chloroplast DNA), eccDNA are circular DNA found in the eukaryotic nuclei of plant and animal (including human) cells. Extrachromosomal circular DNA is derived from chromosomal DNA, can range in size from 50 base pairs to several mega-base pairs in length, and can encode regulatory elements and full-length genes. EcDNA, as the vehicles for oncogene and drug-resistance genes, enables them to be rapidly amplified, and lead to overexpression consequently.For instance, oncogenes EGFR and c-MYC were found in ecDNA and amplified in human cancer tissues than normal tissues. Lacking a high-order chromatin structure, suppressing histone modifications, and insulator shackle make ecDNAs more accessible than their genome counterparts which may facilitate promoter–enhancer interactions, transcription initialization, and achieving additional expression.

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

## Gene fusion - Anchor length
The length of the left and right ends of the read spanning the fusion site, as shown in the right panel in Figure B above.

## Gene fusion - Breakpoint  
The location on the genome where the fusion of two fused genes occurs, such as the site where Gene A (blue) and Gene B (green) are fused in Figure B above.

## Gene fusion - Broken exon (BE) type fusion  
Which refers to the fusion without retaining the original intact exon sequence. For example, in the above figure A, part of the sequence of Exon3 of Gene A is fused with Exon2 of Gene B. In the new gene after fusion, part of the sequence of Exon3 from Gene A is lost.

## Gene fusion - Circular RNAs (false-positive)
Circular RNAs are transcripts which are spliced in a non-canonical mannor such that the 5' and 3' ends are ligated together, yielding a closed ring. In RNA-Seq data they manifest as intragenic duplications with both breakpoints at splice-sites. It is very hard to distinguish genomic duplications from circular RNAs without whole-genome sequencing data. Similar to read-through fusions, these transcripts are abundant in healthy tissue and affect many genes. Arriba's blacklist effectively filters the majority of such events, but occasionally some pass the filter. For this reason, they are tagged as duplication/non-canonical_splicing to warn the user about a potential false positive.

## Gene fusion - Intact exon (IE) type fusion  
The original exon is retained intact after the fusion, and the original exon structure is not affected. For example, in the above figure A, Exon2 of Gene A and Exon1 of Gene B are fused and the sequences of the two exons are retained intact.

## Gene fusion - Long insert size
The longer distance between two reads in double-end sequencing mate-pair sequencing, generally several kilobases or even longer.

## Gene fusion - Read-through fusions (false-positive)
Even in healthy tissue it happens very frequently that transcripts are produced which are composed of exons from two neighboring genes on the same strand. These transcripts often contain all but the last exon of the 5' gene and all but the first exon of the 3' gene. They are caused by RNA polymerases missing the stop sign at the end of the 5' gene and "reading through" until reaching the next stop sign at the end of the 3' gene. The splicesome then removes the intergenic region between the two genes, yielding a chimeric transcript.

For some genes, this is a common phenomenon, which is not reflected in the gene annotation, however, and therefore appears like an aberrant transcript at first glance. The SLC45A3:ELK4 fusion, which has been discovered in prostate cancer and has also been described in benign prostate tissue, is one such example. It is risky to discard all read-through events, however, because they might be the result of a focal deletion, which fuses two neighboring genes together. For example, the GOPC:ROS1 fusion in glioblastoma multiforme is caused by a 240 kb deletion. Arriba uses a comprehensive blacklist trained on collections of RNA-Seq samples from healthy tissue to remove likely harmless transcript variants. Rare transcript variants still bypass this filter often enough. Read-through events are therefore tagged as deletion/read-through to make the user aware of a potential false positive. Further evidence should be sought to substantiate the prediction, such as:

(1) high expression of the transcript; (2) a deletion detected in DNA-Seq data as a structural variant or copy-number alteration; (3) multiple deletion/read-through events affecting the same region (it is unlikely that many events pass the blacklist).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/770f4c4c74ea13817d2474b553e3e75178a90e74/imgs/read-t.png">
</div>

## Gene fusion - Short insert size
The shorter distance between two reads in double-end sequencing paired-end sequencing, generally a few hundred bp.

## Gene fusion - Spanning read
a paired-end read that matches across the fusion site to two fusion genes, respectively, such as the pair of reads matching to Gene A (blue) and Gene B (green) in Figure B above.

## Gene fusion - Split read
A read that matches exactly to the fusion site, as shown in the right panel in Figure B above.

## Genotype phasing
Phasing is the process of inferring haplotypes from genotype data. This information is often important for understanding gene expression patterns for genetic disease research. The popular NGS sequencing technology is to mix the sequences together, and after sequencing, we cannot directly distinguish which of these sequences is the parental source and which is the maternal source. We usually only detect which variants are present in the genome and the base composition of these variants (pure and heterozygous), which is usually referred to as genotype. This distinction can only be achieved with Phasing.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/c15ad61ef3b8af4302985375a7cdec3f7fccdeaa/imgs/old_hap_lecture.jpg">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/c15ad61ef3b8af4302985375a7cdec3f7fccdeaa/imgs/phase.png">
</div>

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

## Illumina adapter portfolio
The first step of sequencing is to construct a library from DNA or RNA. A library contains DNA inserts flanked on each side by an adapter, as shown below:

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/34a9917508fc32c6c1e889f9ad86d50d50fa1e3e/imgs/Illumina%20adapter%20portfolio_1.png">
</div>

Schematic representation of a dual-indexed library fragment

Adapters contain:

Sequences that allow the library to bind and generate clusters on the flow cell (p5 and p7 sequences)  

Sequencing primer binding sites to initiate sequencing (Rd1 SP and Rd2 SP)  

Index sequences (Index 1 and, where applicable, Index 2), which are sample identifiers that allow multiplexing/pooling of multiple samples in a single sequencing run or flow cell lane.  

## Intergenic region
An intergenic region is a stretch of DNA sequences located between genes. Intergenic regions may contain functional elements and junk DNA. Intergenic regions should not be confused with intragenic regions (or introns), which are short, non-coding regions that are found within genes, especially within the genes of eukaryotic organisms.Intergenic regions contain a number of functional DNA sequences such as promoters and regulatory elements, enhancers, spacers, and (in eukaryotes) centromeres. They may also contain origins of replication, scaffold attachment regions, and transposons and viruses.[2] Non-functional DNA elements such as pseudogenes and repetitive DNA, both of which are types of junk DNA, can also be found in intergenic regions—although they may also be located within genes in introns.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/5efcaaf12e6dff7c54ad810d73c281f33d009c2b/imgs/intergenic.png">
</div>

## Intron retention 
An overview of the intron retention (IR) mechanism: different isoforms can be produced from a single gene through AS. (A), Isoforms with introns fully spliced are sent out of the nucleus for translation. Intron-retaining isoforms (IRIs) can be generated through IR (no intron retention): (B), In most cases, the IRIs are degraded by the nonsense-mediated decay (NMD) pathway, the reason being that retained introns often contain premature termination codons (PTCs) that can trigger NMD (with intron retention): (C), In some cases, the IRIs are detained in the nucleus, and in response to stimuli these IRIs can undergo further splicing to remove the retained intron, before being exported out of nucleus for translation (with intron retention): (D), In the case of cytoplasmic splicing, IRIs are shuttled to the cytoplasm for preservation and may be subject to further splicing (with intron retention): (E), In yet another case, IRIs escape from the NMD pathway and are translated into protein isoforms, which, compared with normal protein isoforms, are often truncated and may lose domains; however, it could also be that the alternative protein isoforms include extra domains formed by the amino acid sequences translated from retained introns (with intron retention).  

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/89c8759153191b960c4f4c064d2967a6965c0a87/imgs/fgene-11-00586-g001.jpg">
</div>

## Linear Amplification via Transposon Insertion (LIANTI)  
A linear whole genome amplification (WGA) method.LIANTI achieved linear amplification of the whole genome for the first time, enabling more uniform and accurate amplification.  

Genomic DNA is randomly fragmented and tagged by Tn5 transposon insertion containing T7 promoter sequence, and the resulting DNA fragments are linearly amplified into RNAs by T7 in vitro transcription. Following reverse transcription and second strand synthesis, double-stranded DNA amplicons are formed representing the linear amplification product of the original genomic DNA, which is suitable for DNA library preparation and sequencing.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/c73f61d8e0f1b6c68806bbba654f36e89de055bf/imgs/v2-8a02a512df398452595024dfbfc7d1a7_r.jpg">
</div>

<b>Experimental workflow</b>  

1. Cell lysis. Single cells are placed into PCR tubes containing lysis buffer by mouth pipetting, flow sorting, micromanipulator, laser dissection or microfluidic devices. Cells are subsequently lysed by Qiagen protease digestion. If the starting material is small amount of genomic DNA instead of single cells, cell lysis step can be skipped.

2. LIANTI transposome. LIANTI transposome is made by mixing equal molar of Tn5 transposase and LIANTI transposon DNA. The sequence of LIANTI transposon DNA is: 5'/Phos/CTGTCTCTTATACACATCTGAACAGAATTTAATACGACTCACTATAGGGAGATGTGTATAAGAGACAG-3' After self annealing, the LIANTI transposon DNA consists of a 19-bp double-stranded region for Tn5 transposase binding and dimerization, and a 30-nt single-stranded loop containing T7 promoter sequence.

3. Tn5 transposition. Genomic DNA from a single cell is randomly fragmented and tagged by LIANTI transposome insertion during transposition reaction.

4. Gap filling. After Tn5 transposition, both ends of each DNA fragment are gap filled and extended by DNA polymerase extension, converting single-stranded loops into double-stranded T7 promoters on both ends of each fragment. The residue Tn5 transposase and DNA polymerase are subsequently removed by protease digestion, followed by heat inactivation of the protease.

5. In vitro transcription (IVT) linear amplification. Still within the same PCR tube, overnight IVT reaction is assembled, including standard IVT buffer, NTPs, T7 RNA polymerase, RNase inhibitor, DMSO, etc.

6. Reverse transcription (RT). After overnight IVT, RNAs representing the linearly amplified products of the original genomic DNA are column purified, self primed on the 3′ end, and reverse transcribed. RNase digestion is carried out to convert the double-stranded DNA-RNA hybrids into single-stranded DNA.

7. Second strand synthesis (SSS). Taking advantage of the 19-bp specific sequence on the 3' end of each single-stranded DNA, SSS step is performed by specific priming and DNA polymerase extension. The resulting double-stranded DNA fragments are LIANTI amplicons linearly amplified from the original single-cell genomic DNA, with unique molecular barcodes attached on each amplicon.

8. Library prep and sequencing. Depending on transposome insertion density and specific applications, LIANTI amplicons can be subject to optional sonication, before proceeding to standard library prep pipelines (i.e. NEBNext Ultra II DNA Library Prep Kit for Illumina).

<b>Very important: sometimes it is hard to understand this process, just remeber two tips: (1) a Tn5 transposase will always recognize and combine the Tn5 transposon end sequence; (2) when a Tn5 transposase dimers (two Tnp together) each contains a sequence (whether it is adaptor P5, P7 or it is T7 protor), it looks like that the Tn5 transposase will insert two sequences at the 5 end and 3 end of the target DNA respectively. Infact that is not the truth. One Tnp will only insert the seq it carries into the DNA, this is because three kinds of sequences will be produced during the process </b>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/f77e783a52b182b3ed6fcf7e414a7d740c1aa204/imgs/lianti_2.png">
</div>

<b>Adavantages</b>  
1. Accurate detection of single-cell copy-number variations (CNVs) with kilobase resolution  
2. Accurate detection of single-cell single-nucleotide variations (SNVs)  

## LNCaP cells
LNCaP cells are a cell line of human cells commonly used in the field of oncology. LNCaP cells are androgen-sensitive human prostate adenocarcinoma cells derived from the left supraclavicular lymph node metastasis from a 50-year-old caucasian male in 1977. They are adherent epithelial cells growing in aggregates and as single cells.

## Melanoma
Melanoma, the most serious type of skin cancer, develops in the cells (melanocytes) that produce melanin — the pigment that gives your skin its color.

## MHC class I
MHC class I molecules are expressed in all nucleated cells and also in platelets—in essence all cells but red blood cells. It presents epitopes to killer T cells, also called cytotoxic T lymphocytes (CTLs). A CTL expresses CD8 receptors, in addition to T-cell receptors (TCR)s. When a CTL's CD8 receptor docks to a MHC class I molecule, if the CTL's TCR fits the epitope within the MHC class I molecule, the CTL triggers the cell to undergo programmed cell death by apoptosis. Thus, MHC class I helps mediate cellular immunity, a primary means to address intracellular pathogens, such as viruses and some bacteria, including bacterial L forms, bacterial genus Mycoplasma, and bacterial genus Rickettsia. In humans, MHC class I comprises HLA-A, HLA-B, and HLA-C molecules.

## P5 and P7 adaptors
Regardless of the library construction method, submitted libraries will consist of a sequence of interest flanked on either side by adapter constructs. On each end, these adapter constructs have flow cell binding sites, P5 and P7, which allow the library fragment to attach to the flow cell surface. All Paired-End Format sequencing on the HiSeq and All sequencing of any type on the MiSeq MUST HAVE FULL-LENGTH P5 and P7 sequences . (some of the small RNA libraries and alternative genomic library constructions use a partial P7, this is not supported by the HiSeq PE and MiSeq.)

## Precision and Recall
Precision is calculated by dividing the true positives by anything that was predicted as a positive. Recall (or True Positive Rate) is calculated by dividing the true positives by anything that should have been predicted as positive.

## Promoter
A promoter, as related to genomics, is a region of DNA upstream of a gene where relevant proteins (such as RNA polymerase and transcription factors) bind to initiate transcription of that gene. The resulting transcription produces an RNA molecule (such as mRNA).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/119dccf8f9f3b544668ba5945fb3206ac92180a5/imgs/promoter2.png">
</div>

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
  
## Structural variation  
Genomic structural variation is the variation in structure of an organisms chromosome. It consists of many kinds of variation in the genome of one species, and usually includes microscopic and submicroscopic types, such as deletions, duplications, copy-number variants, insertions, inversions and translocations. Originally, a structure variation affects a sequence length about 1kb to 3Mb, which is larger than SNPs and smaller than chromosome abnormality (though the definitions have some overlap). However, the operational range of structural variants has widened to include events > 50bp.The definition of structural variation does not imply anything about frequency or phenotypical effects. Many structural variants are associated with genetic diseases, however many are not.Recent research about SVs indicates that SVs are more difficult to detect than SNPs.

## T7 promoter and T7 RNA polymerase  
T7 is a promoter (T7 promoter), a strong promoter from T7 phage that can specifically respond to T7 RNA polymerase, and is a sequence that initiates transcription of the T7 phage gene.  

T7 RNA polymerase is highly promoter specific and only transcribes DNA or DNA copies downstream of the T7 promoter in T7 phages. t7 RNA polymerase selectively activates transcription of the T7 phage promoter and synthesizes mRNA at a rate about 5 times faster than that of the common E. coli RNA polymerase.  

In vitro transcription using T7 bacteriophage polymerase is widely used in molecular biology. T7 RNA polymerase is very selective and efficient, resulting both in a high frequency of transcription initiation and effective elongation. These features result in an RNA elongation that is approximately five-fold faster than for E. coli RNA polymerase. 

## Tn5 transposon  
Tn5 transposons were discovered in <i>Escherichia coli</i> and consist of a core sequence encoding three antibiotics (neomycin, bleomycin, and streptomycin) and two inverted IS50 sequences, IS50L and IS50R, which encode a Tn5 transposase (Tnp). IS50 has two pairs of 19-bp inverted ends that are outside ends (OEs) and inside ends (IEs). These inverted OEs are target sites of Tn5 transposase. When transposition occurs, transposases bind to the OEs of the Tn5 transposon, forming Tnp-OE complexes and then the two complexes join together. The C-terminus of Tnp interacts and dimerizes to form a synaptic complex that has the ability to cleave DNA. Tnps that bind to the right and left ends are responsible for catalyzing the hydrolysis of the phosphodiester bond at the left and right ends, respectively. Tnp activates water molecules that hydrolyze the DNA strand and forms a 3′-OH nucleophilic group at the 5′ ends of transposons, which in turn attacks the complementary strand to form a hairpin structure. Finally, the synaptic complex captures target DNA and finishes the strand transfer by nucleophilic attack on both strands of the target DNA with 3’ OH groups of the Tn5 transposon.

Outside end（OE）sequences composed of 19 bases that are recognized by Tnp and involved in the transposition of the entire Tn5.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/34cc72ed144af5775b596e14bd517c3ab844a8b7/imgs/tn52.png">
</div

.
## Tn5 transposase (Tnp)
Linear Amplification via Transposon Insertion (LIANTI) uses Tn5 transposase for linear amplification, haploid typing, and structural variation detection. The “cut and paste” function of Tn5 is widely used in genomics research. Subsequently, studies have shown that only OE sequences and Tn5 transposases are required for transposition in vitro. Tn5 transposases can randomly insert adaptors/barcodes into DNA, and the resulting DNA molecules are ready for PCR amplification and sequencing.

## Tn5 transposase (Tnp) for NGS sequencing library construction  
The key of transposase method library building technology is Tn5 transposon, a bacterial transposon, which is essentially a DNA fragment containing several resistance genes and edited transposase genes. The traditional library construction method requires DNA fragmentation, end repair, splice ligation, library amplification, and multi-step purification, etc. When Tn5 is used for library construction, the multi-step reactions of DNA fragmentation, end repair, and splice ligation can be transformed into a one-step reaction, shortening the library construction time.
  
By combining the P5 and P7 end partial junction sequences (Adapter 1/2) with transposon end sequences (OE), Tnp recognizes transposon ends to form a Tn5 transposon complex with P5 and P7 end partial junctions (containing two dimer-forming Tnp and two partial P5, P7 sequences). The complex recognizes the target sequence of the acceptor DNA, cuts the acceptor DNA (namely ,the DNA you want to sequence), and inserts the carried donor DNA to form DNA with a P5 part adapter Adapter 1 at one end and a P7 part adapter Adapter 2 at the other end, which is then added by PCR barcoding, and the rest of the linker form a DNA library with complete linkers at the P5 and P7 ends. The cut formed by the transposition is then patched using DNA polymerase. This product is amplified by two primer pairs, N5 (N5XX) and N7 (N7XX), and P5 and P7 (PCR Primer Mix), and the product is sorted and purified to form a sequencer-ready library, with N5 (N5XX) and N7 (N7XX) responsible for adding the barcode, and P5 and P7 (PCR Primer Mix) responsible for recovering the full-length sequences of P5 and P7.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/26d617ea1dda100a7ef4c3726faba077d01a7317/imgs/1656323644241570.png">
</div
  
The final structure:

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/26d617ea1dda100a7ef4c3726faba077d01a7317/imgs/library-3.png">
</div  

.

## Tn5 transposase (Tnp) for PacBio sequencing library construction (PCR-free)
A schematic overview of the SMRT-Tag approach. Hairpin adaptor-loaded triple mutant Tn5 transposase is loaded and used to fragment DNA into 2 – 10 kilobase (kb) fragments. After removing Tn5 transposase, an optimized gap repair is used to fill the resulting 9 bp gaps on either side of the molecule, and an exonuclease treatment is used to purify repaired covalently closed templates, which are sequenced on the PacBio Sequel II instrument. 
  
Crucially, combining tagmentation with optimized gap repair allowed the streamlined creation of PacBio libraries from 80 – 100 ng DNA (a minimum of 15,000 human cell equivalents) compared to current protocols that require > 0.5 – 5 μg DNA (a minimum of ∼200,000 human cell equivalents). 

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/3cc296e4e1b223f9e3d836fbd5f6fb0b54ca2c3e/imgs/pc.jpg">
</div

 .
## Transporter associated with antigen processing (TAP)
Transporter associated with antigen processing (TAP) protein complex belongs to the ATP-binding-cassette transporter family. It delivers cytosolic peptides into the endoplasmic reticulum (ER), where they bind to nascent MHC class I molecules.

## Untranslated region (UTR)
Refers to either of two sections, one on each side of a coding sequence on a strand of mRNA. If it is found on the 5' side, it is called the 5' UTR (or leader sequence), or if it is found on the 3' side, it is called the 3' UTR (or trailer sequence). mRNA is RNA that carries information from DNA to the ribosome, the site of protein synthesis (translation) within a cell. The mRNA is initially transcribed from the corresponding DNA sequence and then translated into protein. However, several regions of the mRNA are usually not translated into protein, including the 5' and 3' UTRs.

In protein-coding genes, the exons include both the protein-coding sequence and the 5′- and 3′-untranslated regions (UTR). 

Normal stop codons and the 3' UTR are usually located in the last exon of the sequence 

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/b572cb0bb2f36e30c3b64cae1fb1b19fecc10a3e/imgs/UTR2.png">
</div

.
## Variant allele frequency (VAF)
The fraction of genome copies in the (tumor or control) sample affected by the variant, is either 0.0, 0.5, or 1.0 for germline variants, reflecting absence, heterozygosity, or homozygosity, respectively. In contrast, allele frequencies of somatic variants vary across the whole range from 0.0 to 1.0, depending on the clonal structure of the tumor sample and its impurity (the ratio of healthy genome copies in the tumor sample).
  
The variant allele frequency (VAF; also known as variant allele fraction) is used to infer whether a variant comes from somatic cells or inherited from parents when a matched normal sample is not provided. A variant is potentially a germline mutation if the VAF is approximately 50% or 100%.
  
## VCF (file format)
VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.

```
##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3  
```

This example shows in order (1) a good simple SNP, (2) a possible SNP that has been filtered out because its quality is below 10, (3) a site at which two alternate alleles are called, with one of them (T) being ancestral (possibly a reference sequencing error), (4) a site that is called monomorphic reference (i.e. with no alternate alleles), and (5) a microsatellite with two alternative alleles, one a deletion of 3 bases (TCT), and the other an insertion of one base (A). Genotype data are given for three samples, two of which are phased and the third unphased, with per sample genotype quality, depth and haplotype qualities (the latter only for the phased samples) given as well as the genotypes. The microsatellite calls are unphased.

### 1. Meta-information lines
File meta-information is included after the ## string, often as key=value pairs.

The ‘fileformat’ field is always required and should detail the VCF format version number. For example, for VCF version 4.0, this line should read:

##fileformat=VCFv4.0

It is strongly encouraged that information lines describing the INFO, FILTER and FORMAT entries used in the body of the VCF file be included in the meta-information section. Although they are optional, if these lines are present then they must be completely well-formed.

INFO fields should be described as follows (all keys are required):

##INFO=<ID=ID,Number=number,Type=type,Description=”description”>

Possible Types for INFO fields are: Integer, Float, Flag, Character, and String.

The Number entry is an Integer that describes the number of values that can be included with the INFO field. For example, if the INFO field contains a single number, then this value should be 1. However, if the INFO field describes a pair of numbers, then this value should be 2 and so on. If the number of possible values varies, is unknown, or is unbounded, then this value should be ‘.’. Possible Types are: Integer, Float, Character, String and Flag. The ‘Flag’ type indicates that the INFO field does not contain a Value entry, and hence the Number should be 0 in this case. The Description value must be surrounded by double-quotes.

FILTERs that have been applied to the data should be described as follows:

##FILTER=<ID=ID,Description=”description”>

Likewise, Genotype fields specified in the FORMAT field should be described as follows:

##FORMAT=<ID=ID,Number=number,Type=type,Description=”description”>

Possible Types for FORMAT fields are: Integer, Float, Character, and String.

### 2. The header line syntax
The header line names the 8 fixed, mandatory columns. These columns are as follows:

#CHROM
POS
ID
REF
ALT
QUAL
FILTER
INFO
  
If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number of sample IDs. The header line is tab-delimited.

### 3. Data lines
### Fixed fields
There are 8 fixed fields per record. All data lines are tab-delimited. In all cases, missing values are specified with a dot (“.”). Fixed fields are:

CHROM chromosome: an identifier from the reference genome. All entries for a specific CHROM should form a contiguous block within the VCF file.(Alphanumeric String, Required)
  
POS position: The reference position, with the 1st base having position 1. Positions are sorted numerically, in increasing order, within each reference sequence CHROM. (Integer, Required)
  
ID semi-colon separated list of unique identifiers where available. If this is a dbSNP variant it is encouraged to use the rs number(s). No identifier should be present in more than one data record. If there is no identifier available, then the missing value should be used. (Alphanumeric String)
  
REF reference base(s): Each base must be one of A,C,G,T,N. Bases should be in uppercase. Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String. For InDels, the reference String must include the base before the event (which must be reflected in the POS field). (String, Required).
  
ALT comma separated list of alternate non-reference alleles called on at least one of the samples. Options are base Strings made up of the bases A,C,G,T,N, or an angle-bracketed ID String (”<ID>“). If there are no alternative alleles, then the missing value should be used. Bases should be in uppercase. (Alphanumeric String; no whitespace, commas, or angle-brackets are permitted in the ID String itself)
  
QUAL phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong). If ALT is ”.” (no variant) then this is -10log_10 p(variant), and if ALT is not ”.” this is -10log_10 p(no variant). High QUAL scores indicate high confidence calls. Although traditionally people use integer phred scores, this field is permitted to be a floating point to enable higher resolution for low confidence calls if desired. (Numeric)
  
FILTER filter: PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. “q10;s50” might indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total number of samples. “0” is reserved and should not be used as a filter String. If filters have not been applied, then this field should be set to the missing value. (Alphanumeric String)
  
INFO additional information: (Alphanumeric String) INFO fields are encoded as a semicolon-separated series of short keys with optional values in the format: <key>=<data>[,data]. Arbitrary keys are permitted, although the following sub-fields are reserved (albeit optional):
  
AA ancestral allele
  
AC allele count in genotypes, for each ALT allele, in the same order as listed
  
AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
  
AN total number of alleles in called genotypes
  
BQ RMS base quality at this position
  
CIGAR cigar string describing how to align an alternate allele to the reference allele
  
DB dbSNP membership
  
DP combined depth across samples, e.g. DP=154
  
END end position of the variant described in this record (esp. for CNVs)
  
H2 membership in hapmap2
  
MQ RMS mapping quality, e.g. MQ=52
  
MQ0 Number of MAPQ == 0 reads covering this record
  
NS Number of samples with data
  
SB strand bias at this position
  
SOMATIC indicates that the record is a somatic mutation, for cancer genomics
  
VALIDATED validated by follow-up experiment
  
etc. The exact format of each INFO sub-field should be specified in the meta-information (as described above).

Example for an INFO field: DP=154;MQ=52;H2. Keys without corresponding values are allowed in order to indicate group membership (e.g. H2 indicates the SNP is found in HapMap 2). It is not necessary to list all the properties that a site does NOT have, by e.g. H2=0.

### Genotype fields
If genotype information is present, then the same types of data must be present for all samples. First a FORMAT field is given specifying the data types and order. This is followed by one field per sample, with the colon-separated data in this field corresponding to the types specified in the format. The first sub-field must always be the genotype (GT).

As with the INFO field, there are several common, reserved keywords that are standards across the community:

GT genotype, encoded as alleles values separated by either of ”/” or “|”, e.g. The allele values are 0 for the reference allele (what is in the reference sequence), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X, mitochondrion, only one allele value should be given. All samples must have GT call information; if a call cannot be made for a sample at a given locus, ”.” must be specified for each missing allele in the GT field (for example ./. for a diploid). The meanings of the separators are:
/ : genotype unphased
| : genotype phased

DP read depth at this position for this sample (Integer)

FT sample genotype filter indicating if this genotype was “called” (similar in concept to the FILTER field). Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters that fail, or ”.” to indicate that filters have not been applied. These values should be described in the meta-information in the same way as FILTERs (Alphanumeric String)

GL : three floating point log10-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Numeric)

GQ genotype quality, encoded as a phred quality -10log_10p(genotype call is wrong) (Numeric)
 
HQ haplotype qualities, two phred qualities comma separated (Numeric)
  
If any of the fields is missing, it is replaced with the missing value. For example if the format is GT:GQ:DP:HQ then A|A:.:23:23,34 indicates that GQ is missing. Trailing fields can be dropped (with the exception of the GT field, which should always be present).

Additional Genotype fields can be defined in the meta-information. However, software support for such fields is not guaranteed.

## XS tag in SAM file (only for Tophat and HISAT, not BWA-MEM, BWA-MEM also has the XS tag, but with other meaning)
Where your sequencing is unstranded, mappers can still add a strand to a read if it crosses a splice site, and that splice site has a canonical splice site sequence - it can use this sequence to work out which way the RNA is going. Mappers such as TopHat (which cufflinks was originally designed to work with) only accepted splice sites with one of the canonical splice sites (i'm guessing), thus all spliced reads would contain the XS tag. This is definately true for HiSat.









