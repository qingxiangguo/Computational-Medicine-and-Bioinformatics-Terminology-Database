# Computational-Medicine-and-Bioinformatics-Terminology-Database

A collection of bioinformatics terms explained - because somewhere out there, someone is staring at 'epistasis' for the fifth time today.

# Contributors

Qingxiang (Allen) Guo  
Postdoctoral Fellow  
Northwestern University, Feinberg School of Medicine
Qingxiang.guo@northwestern.edu

# What's here:

Bioinformatics terms in plain English. Born from countless moments of quietly Googling during lab meetings.
Actively maintained. Contributions welcome - misery loves company.

# Terminology Content

## 3C-qPCR (one to one)

### What is 3C-qPCR?
3C-qPCR (Chromosome Conformation Capture combined with Quantitative PCR) is a technique used to analyze the interaction frequency between different genomic regions. It combines the spatial proximity detection capability of 3C technology with the quantitative analysis advantage of qPCR, allowing researchers to accurately determine the physical proximity of specific DNA fragments.

### How to Conduct a 3C-qPCR Experiment?
1. **Cell Nuclei Preparation and DNA Crosslinking**: Cells are fixed with formaldehyde, crosslinking spatially proximate DNA fragments.
2. **DNA Digestion**: Crosslinked DNA is digested with a restriction enzyme, producing sticky-ended DNA fragments.
3. **DNA Fragment Ligation**: The digested DNA fragments are randomly ligated by DNA ligase, forming ligation products.
4. **DNA Purification and Secondary Digestion** (if necessary): Purify the ligation products and perform a secondary digestion to remove unwanted DNA forms.
5. **Primer Design and qPCR Analysis**: Design specific primers for the anticipated ligation products and analyze them using qPCR.

### Restriction Enzyme Selection in 3C-qPCR Experiments

- **First Enzyme**: Choose enzymes like HindIII or EcoRI for their even distribution and high cutting efficiency across the target genomic areas. Optimizing the first enzyme's digestion efficiency is crucial for comprehensive DNA digestion.

- **Second Enzyme**: Select high-frequency cutters such as CviQI, DpnII, or StyI to minimize background noise. Ensure the second enzyme's sites are not near the first enzyme's to avoid digestion interference.

### Primer Design for qPCR
Primers should be designed at the ends of the expected ligation products, ensuring one primer is on the 5' side of the enhancer's restriction site and the other on the 5' side of the promoter's restriction site. This design ensures that primers can specifically amplify the ligation product between the enhancer and promoter.

### CT Value Comparison and ΔΔCT Calculation Example
Suppose you conducted a 3C-qPCR experiment and obtained the following CT values:

- **Control Group** (without IGF1R knock-out):
  - GAPDH CT value: 15
  - Enhancer-Promoter Interaction Region CT value: 23
- **Experimental Group** (with IGF1R knock-out):
  - GAPDH CT value: 15
  - Enhancer-Promoter Interaction Region CT value: 25

**Calculation Steps**:

1. **Calculate the ΔCT for each sample**:
   - ΔCT (Control Group) = 23 (Interaction Region CT value) - 15 (GAPDH CT value) = 8
   - ΔCT (Experimental Group) = 25 (Interaction Region CT value) - 15 (GAPDH CT value) = 10

2. **Calculate the ΔΔCT**:
   - ΔΔCT = ΔCT (Experimental Group) - ΔCT (Control Group) = 10 - 8 = 2

3. **Calculate the Relative Expression**:
   - Relative Expression = 2^(-ΔΔCT) = 2^(-2) = 0.25

This means that, compared to the control group, the interaction frequency between the enhancer and MYC promoter in the experimental group decreased, specifically to 25% of the original frequency. This result indicates that knocking out IGF1R might reduce the physical proximity between the enhancer and MYC promoter, potentially affecting MYC gene expression.

### Incorporating Positive Control into the Simulation
Let's assume we have a positive control region (Control Interaction Region) whose interaction frequency remains unchanged under all conditions, providing its CT values as an example:

- **Positive Control Region CT Values**:
  - Control Group CT value: 20
  - Experimental Group CT value: 20

Since this positive control region is not affected by IGF1R, its CT values should remain consistent across the control and experimental groups, providing us with a baseline to correct experimental variations.

**Correction Steps**:
- **Calculate the ΔCT for the positive control** (using the example values directly as the positive control CT values are the same in both groups):
  - ΔCT (Positive Control, Control Group) = 20 (Control Interaction Region CT value) - 15 (GAPDH CT value) = 5
  - ΔCT (Positive Control, Experimental Group) = 20 (Control Interaction Region CT value) - 15 (GAPDH CT value) = 5

The consistency of the ΔCT values for the positive control across the control and experimental groups confirms the appropriate selection of the reference

<div align=center>
<img src="imgs/3c.png">
</div>

## 4C-seq (Circular Chromosome Conformation Capture followed by sequencing) (one vs all)

### Overview
4C-seq is a robust method used to investigate the spatial organization of chromosomes. It is based on Inverse PCR, amplifying DNA sequences in close proximity to a known sequence to reveal interactions with unknown genomic sequences.

### Enzyme Digestion
The method requires two rounds of enzyme digestion producing **sticky ends** essential for subsequent ligation.

- **First Restriction Enzyme Digestion**: Utilizes a primary enzyme like DpnII, NlaIII, Csp6I, or MboI, selected based on recognition site distribution near the region of interest.
  
- **Second Restriction Enzyme Digestion**: Employs a secondary enzyme post-ligation to further digest the DNA, ensuring smaller fragment sizes for efficient PCR amplification.

### PCR Primer Design
Primers for 4C-seq PCR include Illumina P5 and P7 adapters for seamless sequencing.

- **Reading Primer**: Features an overhang incorporating a portion of the Illumina sequencing primer hybridization sites, **strategically positioned adjacent to the primary restriction site (RE1)** to maximize sequence capture from the viewpoint.

- **Non-Reading Primer**: Located within 50 bp of the secondary restriction site (RE2) and is **designed with an index** to minimize PCR product size and facilitate multiplex sequencing.

### Two-Step PCR Strategy
A two-phase PCR method amplifies ligated fragments while attaching sequencing adapters:

1. **First PCR Step**: Inverse PCR is conducted with VP-specific primers that have overhangs complementary to adapter sequences.
   
2. **Second PCR Step**: Uses universal primers binding to the first PCR step's overhangs, integrating the full Illumina adapter sequences for direct sequencing.

### Key Points
- **Inverse PCR**: Relies on inverse PCR to infer spatial genomic structure from a known sequence.
- **Sticky Ends Production**: Both enzymatic digestions create sticky ends to assist DNA ligation.
- **PCR Primers**: Primers contain P5 and P7 adapters, vital for sequencing setup.
  
### Library Preparation for Sequencing
The PCR-generated product, now ready for Illumina sequencing.

<div align=center>
<img src="imgs/4c1.jpg">
</div>

<div align=center>
<img src="imgs/4c2.png">
</div>

<div align=center>
<img src="imgs/4c3.png">
</div>

## 5mC, 5hmC, CpG Sites, and CpG Islands

#### 5-Methylcytosine (5mC)
5mC is a form of DNA methylation occurring at the 5th carbon of the cytosine ring in DNA. It's often associated with gene silencing, particularly when present in CpG sites. 5mC is a common epigenetic mark and plays a crucial role in regulating gene expression.

#### 5-Hydroxymethylcytosine (5hmC)
5hmC is formed from 5mC through the addition of a hydroxyl group, marking an intermediate step in DNA demethylation. Though it's a methylated form, it often indicates a transition towards gene activation, showcasing a key role in dynamic gene regulation.

#### CpG Sites
CpG sites are DNA regions where a cytosine nucleotide is followed by a guanine nucleotide. Methylation at these sites, primarily involving 5mC, is vital for controlling gene expression. These sites are significant in epigenetic studies due to their influence on gene activity.

#### CpG Islands
CpG islands are dense CpG site regions typically located near gene promoters. Unlike isolated CpG sites, these islands are usually unmethylated, facilitating active gene expression. However, methylation can occur under certain conditions (like in diseases), leading to gene silencing. CpG islands can undergo demethylation, reactivating previously silenced genes.

#### Interplay
The interplay among 5mC, 5hmC, CpG sites, and CpG islands is central to DNA methylation dynamics and gene regulation:

1. **CpG Sites and 5mC**: In CpG sites, cytosine is often methylated to form 5mC, which is associated with gene suppression. This methylation plays a pivotal role in silencing gene functions.

2. **CpG Islands and Transcriptional Activity**: CpG islands, rich in CpG sites, are usually unmethylated, maintaining transcriptional activity in the associated gene regions.

3. **Methylation of CpG Islands**: When CpG islands undergo methylation, it leads to the suppression of downstream gene expression.

4. **5hmC and Gene Reactivation**: The formation of 5hmC is sometimes linked to the demethylation of CpG islands, potentially leading to the reactivation of previously silenced genes. This shift from 5mC to 5hmC signals a tendency towards demethylation and gene activation.

## Allelic balance
The proportion of reads covering a variant’s location that support the variant. For example, if a variant’s location is covered by 100 reads, of which 25 support the variant and 75 do not, then the variant would have an allelic balance of 25/100 = 0.25.

## Amplification Efficiency Differences Post-Bisulfite Treatment in scCOOL-seq, scNanoCOOL-seq, and scBS-seq

### Key Factors Influencing Efficiency Differences

1. **Chemical Transformation Impact**:
   - Bisulfite treatment chemically transforms unmethylated cytosines (C) into uracil (U), affecting the stability and replication of the DNA strand.
   - Transformed U, acting like thymine (T), pairs with adenine (A) instead of guanine (G), altering the original base pairing pattern of the DNA sequence.

2. **Pairing Stability**:
   - Methylated C-G pairings are more stable than U-A pairings due to the presence of three hydrogen bonds in C-G, compared to two in U-A.
   - Stability of base pairing is crucial for PCR amplification efficiency. Instability may lead to errors or reduced efficiency during amplification.

3. **Primer Specificity and Hybridization Efficiency**:
   - Specific primers are required to hybridize with bases on the DNA template during amplification. Primers designed to hybridize more easily with methylated C may lead to higher amplification efficiency in these regions.
   - Transformed U sites might not be ideal targets for all primers, leading to reduced amplification efficiency in these regions.

4. **Polymerase Selectivity**:
   - DNA polymerases used for DNA amplification may exhibit different efficiencies in processing different base pairings. Some enzymes might be less efficient in incorporating bases opposite transformed U (read as T) compared to C.

### Core Reasons for Efficiency Differences

- Methylated C remains a standard nucleotide, capable of standard base pairing with primers and enzymes.
- Modified U, though capable of pairing with A, may have slightly inferior interactions with primers and polymerases.
- In PCR cycles, these minor efficiency differences accumulate, leading to significant differences in amplification product abundance.
- The key reason lies in the subtle but crucial differences in interactions between the modified non-natural base U and the standard components of PCR, leading to overall amplification efficiency differences.

### Conclusion

- Amplification efficiency differences post-bisulfite treatment between methylated C and unmethylated U (originally C) are primarily due to the chemical transformation's impact on base pairing stability, primer specificity, and DNA polymerase selectivity.
- These factors collectively affect the amplification efficiency and outcomes of bisulfite-treated DNA.

## Androgen Deprivation Therapy (ADT)

Androgen Deprivation Therapy (ADT) reduces androgen levels or blocks their effects using surgical castration, LHRH agonists/antagonists, and anti-androgens. It is typically the first-line treatment for advanced or metastatic prostate cancer.

## ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing)

It is a powerful technique for profiling chromatin accessibility and identifying regions of open chromatin in the genome.

ATAC-seq works by using a hyperactive Tn5 transposase enzyme to tag regions of open chromatin with sequencing adapters. The transposase inserts sequencing adapters into the chromatin accessible regions, which can then be PCR-amplified and sequenced using next-generation sequencing technology.

The resulting sequence data can be used to identify regions of open chromatin and to map transcription factor binding sites, enhancers, promoters, and other regulatory elements. By comparing ATAC-seq data from different cell types or under different conditions, researchers can gain insights into the regulation of gene expression and cellular differentiation.

So what are the common chromatin open regions?

The common chromatin open regions are mainly the promoter upstream of the gene and the distal regulatory elements such as enhancers and silencers. The promoter is the DNA region near the transcription start site (TSS), which contains the transcription factor binding site (TFBS), so the transcription factor can bind to the TFBS on the promoter and recruit RNA polymerase to transcribe the gene. The TSS contains a transcription factor binding site (TFBS), so the transcription factor can bind to the TFBS on the promoter and recruit RNA polymerase to transcribe genes. Enhancers are generally located in the 1 Mb DNA region downstream or upstream of the promoter, and when transcription factors bind to enhancers and make contact with the promoter region, they can promote gene transcription. Conversely, silencers reduce or repress the expression of genes.

<div align=center>
<img src="imgs/promoter3.png">
</div>

Workflow

<div align=center>
<img src="imgs/atac22.jpg">
</div>

Workflow

<div align=center>
<img src="imgs/atac33.png">
</div>

Therefore, ATAC-seq can help identify promoter regions, potential enhancers or silencers, that is, the peaks in ATAC-seq, are often sites where promoters, enhancer sequences, and some trans-regulators bind.

Isn't the body of that gene a chromatin open region? Yes ~ but the chromatin opening in the body region of a gene is not stable; when transcription is performed, the body region is transcribed with each step forward, after opening, and quickly returns to its original state after transcription.

So let's summarize, what can ATAC-seq be used for in the end?

1. Identification of important transcription factors: According to the principle, the chromatin open region captured by ATAC is generally upstream and downstream of the part of the DNA sequence being transcribed, so that we can combine the enriched sequences with motif analysis to identify which transcription factors are involved in the regulation of gene expression. (particularly effective for transcription factors with poor antibody quality)

2. Generate features of transcription factor binding regions (footprinting): After a transcription factor binds to DNA, it occupies a space that prevents the transposase Tn5 enzyme from cutting in other nucleosome-free regions, which leaves a small region, called a footprint, in which reads drop abruptly from high coverage peaks. So ATAC-seq footprints can help us to see the state of transcription factors binding on the whole genome, and are mainly applied to study cellular reprogramming mechanisms, chromatin remodeling factors, domains of epigenetic modifications on diseases, T cell depletion, etc.

3、Generate the epigenome map

4、Get the corresponding accessibility regions in different tissues or under different conditions.

5、Get the nucleosome location

## Avidin 
Avidin is a protein derived from both avians and amphibians that shows considerable affinity for biotin, a co-factor that plays a role in multiple eukaryotic biological processes. Avidin and other biotin-binding proteins, including Streptavidin and NeutrAvidin protein, have the ability to bind up to four biotin molecules.

Avidin is a biotin-binding protein that is believed to function as an antibiotic in the eggs of birds, reptiles and amphibians. Chicken Avidin has a mass of 67,000-68,000 daltons and is formed from four 128 amino acid-subunits, each binding one molecule of biotin. 

Avidin has a very high affinity for up to four biotin molecules and is stable and functional over a wide range of pH and temperature. Avidin is amenable to extensive chemical modification with little to no effect on function, making it useful for the detection and protein purification of biotinylated molecules in a variety of conditions. Because Avidin is easily purified from chicken egg whites, it is very economical to produce (much more so than streptavidin).

## Avidin-biotin interaction
Avidin and other biotin-binding proteins, including Streptavidin and NeutrAvidin protein, have the ability to bind up to four biotin molecules. The Avidin-biotin complex is the strongest known non-covalent interaction (Kd = 10-15M) between a protein and ligand. The bond formation between biotin and Avidin is very rapid, and once formed, is unaffected by extremes of pH, temperature, organic solvents and other denaturing agents. These features of biotin and Avidin – features that are shared by Streptavidin and NeutrAvidin Protein – are useful for purifying or detecting proteins conjugated to either component of the interaction.

## Biotin

Biotin is a vitamin (Vitamin H, Vitamin B7, Coenzyme R) that is present in small amounts in all living cells and is critical for a number of biological processes including cell growth and the citric acid cycle. The valeric acid side chain of the biotin molecule can be derivatized in order to incorporate various reactive groups that facilitate the addition of a biotin tag to other molecules.  Because biotin is relatively small (244.3 Daltons), it can be conjugated to many proteins and other molecules without significantly altering their biological activity. The highly specific interaction of biotin-binding proteins with biotin makes it a useful tool in assay systems designed to detect and target biological analytes.

## Biotinylated dCas9

Biotinylated dCas9, a non-cutting version of Cas9 with a biotin tag, targets DNA precisely without damage, used for genomic research.

### Features

- Retains DNA binding, loses cutting ability.
- Biotin enables easy molecular tagging.

### Uses

Target genome sequence.

## Blunt ends

In blunt ends, both strands are of equal length – i.e. they end at the same base position, leaving no unpaired bases on either strand. 

## Blunt-end cloning overview

Blunt-end cloning is a fundamental molecular biology technique for inserting blunt-ended DNA fragments into blunt-ended vectors. This method allows for the direct cloning of DNA without the need for complementary sticky ends, making it a versatile tool for various genetic manipulations.

### Generating Blunt Ends
Blunt ends are created when double-stranded DNA has nucleotides that pair perfectly across from each other. Methods to generate blunt ends include:
- **DNA Amplification**: Using proofreading DNA polymerases that ensure amplicons have blunt ends.
- **Restriction Enzyme Digestion**: Employing enzymes that cut DNA to produce blunt ends directly.
- **End Modification**: Converting sticky ends to blunt ends through end repair mechanisms.

### Cloning Process
- **Vector and Insert Preparation**: Both need blunt ends for successful ligation. Vectors may undergo dephosphorylation to prevent self-ligation, whereas inserts are prepared via PCR, restriction digestion, or end repair.
- **Ligation**: DNA ligase catalyzes the formation of a phosphodiester bond between the 3’ hydroxyl and 5’ phosphate groups, joining the vector and insert.
- **Transformation and Selection**: Following ligation, the recombinant plasmid is introduced into E. coli, and transformants are selected.

### Advantages
- **Simplicity**: No need for specific sequences at the primer ends.
- **Universality**: Can clone any DNA fragment without complementary sequences.
- **Versatility**: Suitable for sub-cloning, sequencing, and library construction.

### Limitations
- **Orientation**: Non-directional cloning leads to inserts in both orientations.
- **Efficiency**: Lower recombination efficiency compared to sticky-end cloning.
- **Self-Ligation**: Vectors can self-ligate, necessitating additional steps to mitigate.

### Key Steps in Procedure
1. **Insert Preparation**: Via PCR with proofreading polymerases, restriction enzyme digestion or modification of existing ends.
2. **Vector Preparation**: Amplifying with proofreading polymerases, using blunt-end generating enzymes, or modifying sticky ends.
3. **Ligation**: Employing DNA ligase to join insert and vector.
4. **Transformation**: Introducing the ligated plasmid into bacteria for propagation.
5. **Selection**: Identifying and selecting the correct transformants.

<div align=center>
<img src="imgs/bl1.png">
</div>


<div align=center>
<img src="imgs/bl2.png">
</div>

Strategies to generate an insert for blunt-end cloning. 3.1 One of the common ways to generate a blunt-ended insert is to amplify the desired piece of the template DNA using a proofreading polymerase. 3.2 Amplicons generated by a non-proofreading polymerase have overhangs (extra nucleotides on one strand of dsDNA) that need to be removed by PCR polishing. 3.3 Generation of blunt-end insert via restriction digestion with blunt-end generating enzymes. 3.4 If sticky-ends are generated in the workflow, they can be converted to blunt-ends by end repair (Fill in 3.4A and chew back 3.4B). 3.5 Alternative to end repair, mung bean nuclease generates blunt-ends via digestion of 3’ and 5’ overhangs.

<div align=center>
<img src="imgs/bl3.png">
</div>

Strategies to generate a plasmid for blunt-end cloning. 4.1 One of the common ways to generate a blunt-ended vector is to amplify the desired piece of the template DNA using a proofreading polymerase. 4.2 Amplicons generated by a non-proofreading polymerase have overhangs (extra nucleotides on one strand of dsDNA) that need to be removed by PCR polishing. 4.3 If the vector has sticky-ends, they can be converted to blunt-ends either by end repair or mung bean nuclease.

## Blunt-end cloning: Phosphatase Reaction of the Vector Before Ligation

In the context of molecular cloning, treating vectors with phosphatase before ligation is crucial to prevent the vector's self-ligation and to ensure successful cloning. 

### Key Insight

- **T4 DNA Ligase Requirement**: A fundamental requirement for the ligation process, catalyzed by T4 DNA ligase, is the presence of at least one phosphate group at the DNA ends being joined. This enzyme facilitates the formation of a phosphodiester bond between the 3' hydroxyl (OH) group of one nucleotide and the 5' phosphate group of another.

### Phosphatase Treatment Rationale

- By removing the 5' phosphate groups from the vector, phosphatase treatment prevents the vector from self-ligating. This significantly reduces the background of clones without the insert, which are otherwise a common outcome due to vector self-ligation.
- **Critical Balance**: While this treatment is effective in reducing self-ligation, it is essential to remember that for successful ligation to occur, the insert must retain its 5' phosphate group, or the vector must be re-phosphorylated after phosphatase treatment if both were dephosphorylated.

### Enzymes and Strategies

- **Enzymes Used**: Common phosphatases such as Shrimp Alkaline Phosphatase (rSAP), Calf Intestinal Phosphatase (CIP), and Antarctic Phosphatase (AnP) are employed for this purpose.
- **Post-Treatment Considerations**: If both vector and insert end up dephosphorylated, adding a phosphate group back to one of them is necessary for ligation. T4 polynucleotide kinase is typically used for this re-phosphorylation step.

<div align=center>
<img src="imgs/pho.png">
</div>

## Blunt end ligation

Usually, a straight cut creates blunt ends or non-overhanging ends. These ends can be joined using a DNA ligase enzyme. The joining of two blunt ends is called blunt end ligation. This does not need matching or complementary ends for ligation.

## Breakend symbol in BND translocation - square brackets

These 3 elements are combined in 4 possible ways to create the ALT. In each of the 4 cases, the assertion is that s
is replaced with t, and then some piece starting at position p is joined to t. The cases are:
REF ALT Meaning

s t[p[ piece extending to the right of p is joined after t

s t]p] reverse complementary sequence piece extending left of p is joined after t

s ]p]t piece extending to the left of p is joined before t

s [p[t reverse complementary sequence piece extending right of p is joined before t

So what does "piece extending left of p" mean? It is very hard to understand without figures.

There are two key points to understand:

One is that the breakpoints are connected. You have to make sure the connected reads are from 5 to 3 end, so it can be translated properly. In this way, sometimes the positive strand might be connected to the minus strand of another chromosome to make sure that the result sequence is from 5 to 3 end. The other point is that when we describe this connection in the VCF file, we represent the case of a connection on a positive chain (and what is the status of the other strand). In fact, in biological phenomena, both chains are connected.

<div align=center>
<img src="imgs/bnd6.png">
</div>

<div align=center>
<img src="imgs/bnd1.jpeg">
</div>

<div align=center>
<img src="imgs/bnd2.jpeg">
</div>

<div align=center>
<img src="imgs/bnd3.jpeg">
</div>

<div align=center>
<img src="imgs/bnd4.jpeg">
</div>

<div align=center>
<img src="imgs/bnd5.jpeg">
</div>



## Cancer immune Evasion through loss of MHC Class I antigen presentation

Major histocompatibility class I (MHC I) molecules bind peptides derived from a cell's expressed genes and then transport and display this antigenic information on the cell surface. This allows CD8 T cells to identify pathological cells that are synthesizing abnormal proteins, such as cancers that are expressing mutated proteins. In order for many cancers to arise and progress, they need to evolve mechanisms to avoid elimination by CD8 T cells. MHC I molecules are not essential for cell survival and therefore one mechanism by which cancers can evade immune control is by losing MHC I antigen presentation machinery (APM). Not only will this impair the ability of natural immune responses to control cancers, but also frustrate immunotherapies that work by re-invigorating anti-tumor CD8 T cells, such as checkpoint blockade. 

## Castration-resistant prostate cancer (CRPC) 

A form of advanced prostate cancer. With CRPC, the cancer no longer completely responds to treatments that lower testosterone. It shows signs of growth, like a rising PSA (prostate-specific antigen), even with low levels of testosterone.

## CCAT1 Expression in Prostate Cancer

CCAT1, a long non-coding RNA, is upregulated in prostate cancer and is present in multiple prostate cancer cell lines. It is expressed in LNCaP, DU145, PC3, and 22RV1 cells, while notably absent in VCaP cells.

### Impact of CCAT1 Knockdown
- **MTT Assay (Figure B)**: si-CCAT1 treatment in LnCaP and PC3 cells leads to reduced cell viability over time, implying that CCAT1 supports cell survival and proliferation.
- **EdU Staining (Figure C)**: si-CCAT1 decreases DNA synthesis in LnCaP and PC3 cells, as evidenced by fewer EdU-positive cells.
- **Flow Cytometry (Figure D)**: Higher apoptosis rates are observed in si-CCAT1 treated LnCaP and PC3 cells, highlighting CCAT1's role in apoptosis inhibition.
- **Transwell Assays (Figures E and F)**: si-CCAT1 reduces migration and invasion in LnCaP and PC3 cells, indicating CCAT1's involvement in metastatic processes.

<div align=center>
<img src="imgs/ccat1.png">
</div>

<div align=center>
<img src="imgs/ccat2.png">
</div>


## CCK-8 Assay (Cell Counting Kit-8)
The CCK-8 assay is a cell viability and proliferation detection method based on the reduction of the WST-8 tetrazolium salt, assessing cellular health by measuring absorbance.

## cDNA second strand synthesis vs SMART-seq2

Normal cDNA second strand sythesis is not a PCR process. But Takara smart-seq2 containing PCR of cDNA. The first strand cDNA synthesis is typically carried out using reverse transcription (RT) enzymes, which use a reverse transcriptase to synthesize cDNA from an RNA template. The resulting single-stranded cDNA is then converted to double-stranded cDNA through a second strand synthesis step that typically uses DNA polymerase and RNase H to remove the RNA template.

Smart-seq2 is an advanced RNA-Seq methodology that integrates template switching with PCR amplification to create comprehensive, full-length cDNA libraries from individual cells. Post reverse transcription, a template switching oligonucleotide (TSO) is introduced, which contains an rGrGrG sequence at its 3' end for hybridization with the cDNA's complementary sequence. This TSO facilitates template switching by the reverse transcriptase when it reaches the 5' end of the mRNA, allowing the addition of a known sequence to the 3' end of the first-strand cDNA, enhancing cDNA completeness and providing a uniform priming site for PCR amplification. 

After the second strand synthesis, PCR amplification is carried out to amplify the cDNA libraries. The PCR step is designed to specifically amplify only full-length cDNAs, which helps to reduce bias and improve the quality of the final library. The amplified cDNA libraries can then be used for high-throughput sequencing.

So while the second strand synthesis in Smart-seq2 does involve PCR amplification, it is not a standard cDNA second strand synthesis process.

## CAPTURE Technology Summary (one DNA to many RNA, protein)

## CAPTURE Technique Overview

CAPTURE technique integrates CRISPR-Cas9 system with biotin-streptavidin affinity purification to study protein-DNA interactions and chromatin architecture. 

1. **Vector Construction**: 

**Biotinylated dCas9 Vector**: This vector contains the gene sequence for dCas9 fused with a biotin acceptor peptide (BAP) sequence, which can be biotinylated by BirA ligase.
- **BirA Ligase Expression Vector**: BirA ligase, derived from E. coli, specifically biotinylates proteins with a BAP sequence. A separate vector is constructed for expressing BirA ligase in cells.
- **sgRNA Expression Vector**: Contains the specific sgRNA gene sequence for guiding dCas9 to the target genomic location.

2. **Cell Culture**: Grow suitable cell lines for transfection.
3. **Transfection**: Co-transfect cells with vectors for biotinylated dCas9, BirA ligase, and sgRNAs (lentiviruses method by 293T).
4. **Crosslinking**: Fix cells with formaldehyde to stabilize protein-DNA interactions.
5. **Chromatin Preparation and Shearing**: Lyse cells, extract chromatin, and shear it into fragments using sonication.
6. **Affinity Purification**: Capture biotinylated dCas9-bound DNA-protein complexes with streptavidin-coated magnetic beads.
7. **Reverse Crosslinking and DNA Purification**: Reverse crosslinks to free DNA and proteins, followed by DNA purification.
8. **Analysis**: Perform proteomics and genomics analysis to identify interacting proteins and DNA regions.

The CAPTURE technique provides a powerful tool for dissecting complex regulatory networks and understanding chromatin dynamics by enabling the study of specific genomic loci and their interacting factors.

<div align=center>
<img src="/imgs/capture.png">
</div>

## CHAR-seq (Chromatin-Associated RNA sequencing, DNA linker)

ChAR-seq (Chromatin-Associated RNA sequencing) is a comprehensive **all-to-all RNA-DNA interaction** technique designed to detect interactions between RNA and DNA within the genome. At the core of this method is a biotinylated bridge oligonucleotide that facilitates the capture of RNA-DNA contacts.

Here’s a brief overview of the ChAR-seq process:

- **Biotin Bridge**: A key feature of ChAR-seq is the use of a biotin-labeled bridge molecule. This bridge is a double-stranded DNA sequence with one strand containing a biotin modification.
  
- **Reverse Transcription**: RNA molecules in close proximity to the DNA are reverse transcribed to cDNA using the bridge as a primer, effectively linking the RNA sequence to the bridge.

- **Enzymatic Digestion**: The genomic DNA and bridge are both then digested, and the bridge molecule with its attached cDNA is ligated to nearby DNA fragments.

- **Purification and Sequencing**: Following ligation, a second enzymatic digestion step is employed to remove any bridges that have not properly ligated to the genomic DNA but have instead ligated to sequencing adaptors. This is achieved by utilizing a secondary restriction enzyme site engineered into the bridge, which ensures that only the desired chimeric cDNA-DNA molecules are amplified and sequenced.

<div align=center>
<img src="/imgs/ch1.jpeg">
</div>

<div align=center>
<img src="/imgs/ch2.jpeg">
</div>

<div align=center>
<img src="/imgs/ch3.png">
</div>

## Checkpoint blockade immunotherapy

Patients are treated with antibodies that block negative regulatory molecules, such as PD-1/PD-L1 or CTLA4, which normally restrain T cell responses. This kind of therapy can reinvigorate a patient's anti-tumor T cell responses, which then can cause tumors to shrink and even lead to cures in some patients

## ChIA-PET: (All to all DNA-DNA and CHIP-seq) 

Chromatin Interaction Analysis by Paired-End Tag Sequencing (ChIA-PET) is a technique combining chromatin immunoprecipitation (ChIP) with high-throughput sequencing to map long-range chromatin interactions and the DNA-binding sites of proteins of interest.

### Key Steps in ChIA-PET

- **Cross-linking:** Fix chromatin interactions in place using formaldehyde.
- **Chromatin Shearing:** Breakdown the chromatin into manageable pieces.
- **ChIP Enrichment:** Isolate DNA segments that interact with the protein of interest using specific antibodies.
- **Ligation of Linkers:** Add distinct linkers to two separate aliquots of the enriched chromatin.
- **Proximity Ligation:** After linkers are attached, recombine the aliquots to allow ligation between physically interacting DNA fragments.
- **Restriction Digest:** Cut back the DNA at linkers to create short, manageable tags.
- **Sequencing:** Perform high-throughput sequencing on the ligated DNA.
- **Data Analysis:** Identify and characterize the chromatin interactions from the sequencing data.

### Significance of Dual Linker Strategy

Splitting chromatin into two aliquots and ligating different linkers (Linker A and Linker B), then mixing them does not alter the inherent spatial interactions between chromatin fragments. The purpose of this process is to distinguish self-ligation events, where a chromatin fragment loops back on itself, from inter-ligation events, where different chromatin fragments that are proximal in the nucleus are joined. This distinction is crucial for accurate mapping of chromatin interactions.

<div align=center>
<img src="/imgs/chia-pet.png">
</div>

## Chimera rates in Whole Genome Amplification (WGA)  

Defined as the number of reads that are improperly connected (including abnormal fragment size and interchromosomal connection) divided by the total number of mappable reads. Here are the chimera rates for some commercial kits: 

the SigmaAldrich GenomePlex Single Cell Whole Genome Amplification Kit (DOP-PCR) 15%

the Qiagen REPLI-g Single Cell Kit (MDA) 2%

the General Electric (GE) illustra Single Cell GenomiPhi DNA Amplification Kit (MDA)  3%

the Yikon Genomics Single Cell Whole Genome Amplification Kit (MALBAC)  5%

the Rubicon Genomics PicoPLEX WGA Kit (MALBAC-like)  13%

the Bioskryb PTA  About 15% to 7.5%

## Chimeric/split/supplementary reads

They are the same thing.

Chimeric reads occur when one sequencing read aligns to two distinct portions of the genome with little or no overlap. This could be like sequence A mapped to 85156-85257 bp of genome, while part of sequence A mapped to 85273-85320 bp of genome. Then, sequence A is a chimeric read. Chimeric reads are indicative of structural variation. Chimeric reads are also called split reads.

An alignment of a read that cannot be represented as a linear alignment. A chimeric alignment is represented as a set of linear alignments that do not have large overlaps. Typically, one of the linear alignments in a chimeric alignment is considered the “representative” alignment, and the others are called “supplementary” and are distinguished by the supplementary alignment flag. All the SAM records in a chimeric alignment have the same QNAME.

Split alignments, or chimeric alignments, are alignments where part of the read maps to one place, and another part to another. For example, part of a long read may map to chr1 and part of it maps to chr4.

One record is marked as "representative", sometimes also called the "primary" record, while the other components of the split read are marked "supplementary", given the 2048 flag. The "primary" record generally has a SEQ field that represents the entirety of the original read's sequence (with CIGAR soft clipping operators saying which part of that sequence aligned), and the "supplementary alignments" will have SEQ field but sometimes just segments of the original read's sequence with CIGAR hard clipping operators indicating that it is partial.

Supplementary alignments are especially common with long reads, and it can be a signal for structural variants e.g. where two chromosomes are fused together, and parts of the read align to multiple chromosomes, or the split alignment may align to either side of a large deletion, or they may be split to align through an inversion (part of it aligns to the forward strand, part of it to the reverse strand, and again the forward strand)

There is no limitation on how many splits might occur so the split can align to 3, 4, or more different places. Each part of the split puts a new line in the SAM file, and note that all the records also have the same read name, or QNAME (first column of SAM).

<div align=center>
<img src="imgs/326013ae-35ce-4f20-8992-09c89110.png">
</div>

The figure is from: https://cmdcolin.github.io/posts/2022-02-06-sv-sam

## Chimeric read in PacBio

PacBio generates two types of chimeric reads:

1. unsplit subreads, i.e. siamaera
Reads structure looks like this: ----R1--->--A--<--R1.rc--. This happens quite frequently, however, it strongly depends on chemistry and particularly the quality of the libary prep. For siamaera detection, Siamaeric reads usually are separated by short joint sequences (corrupted adapter). Detection is based on blasting and identifying reads with reverse complement self hits. Reads are trimmed to the longest non-chimaeric subsequence without joint sequence.

```
----R--->--J--<--R.rc--  siamaeric read
----R--->                trimmed read

--R->-J-<----R.rc-----   siameric read2
        <----R.rc-----   trimmed read2
```

2. random fusion chimeras

I’m not exactly sure, when or how this happens, but there is a fraction of reads, where random sequences seem to be fused together. Probably some blunt end ligations during library prep, or similar effect … This seems to happen quite rarely, and it is hard to quantify exactly, as there are other effects, that can cause reads to look like chimeras, although they aren’t.

SMRT Bells are blunt end ligated to the DNA fragment being sequenced, sometimes a DNA fragment ligating to another DNA fragment before having a SMRT Bell adapter added. This would generate the classic chimera, a sequenced read being from two random parts of a genome. Note this is random and will not happen at exactly the same location more than once, so it is easily dealt with at the analysis stage. The "siamaeras" or missing adapter look like a sequence followed by the reverse complement of the sequence, and results from the SMRT Bell being missing/not-detected on one end of the insert i.e. you read a sequence forward then backwards without an adapter so the software does not know that the read needs to be split. This can happen due to sample prep - a long overhang forms a hairpin, mimicking a SMRT Bell, or less frequently due to a real SMRT bell being missed in software.

## ChIP-Seq Technology Overview

TFBS, or Transcription Factor Binding Sites, can be located inside promoters or enhancers, meaning transcription factors can bind to both promoters and enhancers.

Promoters are generally upstream of the TSS (Transcription Start Site).

In ChIP-seq, peaks identified using an antibody against a specific transcription factor actually represent TFBS. These peaks might be inside known promoter regions (where gene annotations indicate promoters), or if they are upstream of promoter regions, they could indicate enhancers.

### ChIP-Seq Workflow

ChIP-Seq (Chromatin Immunoprecipitation Sequencing) is a powerful technique used to study protein-DNA interactions within the cell. The process involves several key steps:

1. **Cross-linking**: Proteins are cross-linked to DNA using a fixative like formaldehyde.
2. **Fragmentation**: The chromatin is then fragmented, typically using sonication or enzymatic digestion.
3. **Immunoprecipitation**: Fragmented chromatin is incubated with a specific antibody (primary antibody) that binds to the protein of interest.
4. **Capture of Antibody-Chromatin Complexes**: After incubation, the antibody-chromatin complexes are pulled down using magnetic beads or agarose beads coated with protein A/G.
5. **Washing**: The beads are washed to remove non-specifically bound material.
6. **Reverse Cross-linking**: The protein-DNA cross-links are reversed, usually by heating, to free the DNA.
7. **DNA Purification**: The DNA is then purified from the solution.

### Q&A Related to ChIP-Seq

**Q1: Why is ChIP-Seq also used to study histones?**

- **A1**: ChIP-Seq is not limited to studying transcription factors or DNA-binding proteins. It can also be used to study histones and their modifications. Histones are proteins around which DNA winds, and their modifications can influence gene expression. By using antibodies specific to certain histone modifications, ChIP-Seq can provide insights into the role of these modifications in gene regulation.

**Q2: What is the significance of using H3 or Pol II as positive controls in ChIP-Seq?**

- **A2**: Histone H3 and RNA Polymerase II (Pol II) are commonly used as positive controls in ChIP-Seq. H3 is a core histone protein, and its presence across the genome provides a uniform control signal. Pol II, involved in mRNA synthesis, is enriched at active gene promoters, providing a distinct peak pattern in these regions. These controls help validate the effectiveness of the ChIP process.

**Q3: Why is the negative control IgG selected based on the host species and subtype of the primary antibody?**

- **A3**: The negative control, typically normal IgG, is chosen to match the host species and subtype of the primary antibody to ensure that the behavior of the control mimics that of the primary antibody, except for the specific binding to the target protein. This matching helps in accurately determining the non-specific background signal in the assay.

**Q4: What is the role of protein A/G in ChIP-Seq?**

- **A4**: Protein A and protein G are bacterial proteins that bind to the Fc region of antibodies. In ChIP-Seq, they are used to capture the antibody-chromatin complexes. Beads coated with protein A/G are used to pull down these complexes from the solution. The choice between protein A and G depends on the species and subtype of the antibody used, as they differ in their binding affinities to different antibodies.

**Q5: Why is chromatin fragmented and a portion used as an 'input' control in ChIP-Seq?**

- **A5**: Fragmentation of chromatin is necessary to expose the DNA-binding sites and to create manageable lengths of DNA for sequencing. A portion of this fragmented chromatin, known as 'input,' is set aside before immunoprecipitation. The input sample represents the baseline DNA content and is used in data analysis as a control to account for differences in DNA accessibility and fragmentation efficiency, allowing for normalization and accurate interpretation of ChIP-Seq results.

**Q6: What are the differences between agarose beads and magnetic beads in ChIP, and how do you choose between them?**

- **A6**: Agarose beads and magnetic beads are both used in ChIP assays to pull down the antibody-chromatin complexes, but they have different properties and practical implications:
    - **Agarose Beads**: These beads are a traditional choice for immunoprecipitation. They are usually cheaper and can bind a large amount of antibody. However, the process of separation by centrifugation can be time-consuming and may lead to sample loss or contamination.
    - **Magnetic Beads**: Magnetic beads offer a more modern approach. They allow for a faster and cleaner separation process using a magnet, reducing sample loss and increasing efficiency. They are particularly useful when working with small sample volumes or when a rapid turnaround is needed.
    - **Choosing Between Them**: The choice between agarose and magnetic beads depends on several factors, including the scale of the experiment, available resources, and personal preference. Magnetic beads are generally favored for their ease of use and efficiency, especially in high-throughput settings, but agarose beads can be a cost-effective alternative for larger volumes or when budget constraints are significant.

**Q7: How many biological replicates are generally needed in ChIP-Seq, and how many replicates are required for the input control?**

- **A7**: The number of biological replicates in ChIP-Seq is crucial for ensuring the reliability and reproducibility of the results:

    - **Biological Replicates**: A minimum of two to three biological replicates for each ChIP-Seq experiment is typically recommended to account for biological variability. However, the exact number can vary based on the sample type, the robustness of the ChIP protocol, and the statistical power required for the study. Some high-precision studies may need more replicates to achieve reliable results.
    
    - **Input Control Replicates**: For the input control, which serves as a baseline reference for normalizing ChIP-Seq data, it is generally considered sufficient to have one input control, even when there are multiple biological replicates of the ChIP sample. While it's ideal to match the number of input controls to the number of biological replicates, many studies opt for a single input control due to its inherent stability and lower variability. This single input sample is processed alongside the ChIP samples to control for factors like DNA accessibility and fragmentation efficiency.

In practice, the use of one input control is often a balance between resource optimization and the empirical observation that input DNA tends to have less variability between different biological conditions. Therefore, while having an input control for each biological replicate might provide a slight increase in accuracy, it is often deemed unnecessary for the normalization process in many ChIP-Seq analyses.

<div align=center>
<img src="imgs/chip.png">
</div>

## ChIP-seq and RNA-seq co analysis

In genomic research, ChIP-seq (Chromatin Immunoprecipitation Sequencing) and RNA-seq (RNA Sequencing) are commonly used experimental techniques. When presented side by side in the same chart, they are often used to illustrate the relationship between the state of epigenetic modifications in chromatin (such as methylation, acetylation, etc.) and gene expression.

Specifically, ChIP-seq is primarily used to detect the distribution of specific proteins or epigenetic markers on chromatin, such as histone modifications or transcription factor binding sites. On the other hand, RNA-seq is used to quantitatively measure gene expression levels.

In a chart with aligned coordinates, the x-axis represents the position on the genome, and the y-axis represents the signal intensity at that position. The peaks in ChIP-seq indicate the presence of specific epigenetic markers or protein binding at that location, while the peaks in RNA-seq represent gene expression at that location. By jointly analyzing these two types of data, the relationship between gene expression and chromatin status can be derived. For example, a specific histone modification might be enriched in areas with high gene expression, suggesting that this modification might facilitate gene expression.

Additionally, this analysis can also help to reveal the mechanisms of gene regulation. For example, if the binding sites of a transcription factor (detected by ChIP-seq) highly correlate with the expression level of a gene (detected by RNA-seq), it can be inferred that this transcription factor may play a crucial role in regulating the expression of this gene.

## CHIRP-seq (RNA-chromatin interaction, one RNA to many DNA, can also detect RNA and proteins)

CHIRP-seq (Chromatin Isolation by RNA Purification followed by sequencing) is a sophisticated method utilized to elucidate the interaction sites between chromatin-associated RNAs and the genomic DNA. This technique is pivotal for understanding the complex roles of long non-coding RNAs (lncRNAs) in gene regulation and chromatin organization.

### Controls and primers in CHIRP-seq

- For RNA, we can carry out ChIRP-RT-qPCR, where typically a target-specific probe and a LacZ negative control probe are employed. The amplified regions are then detected using primers specific for the RNA of interest for all probes.

- In the case of DNA, ChIRP-qPCR is performed using both the target-specific and LacZ negative control probes. Primers designed to amplify the DNA regions where RNA is thought to bind are used.

### Experimental Probes

For the target lncRNA under investigation, probes must be designed in even and odd sets for cross-validation purposes.

- **Odd Probes**: Bind to the odd-numbered regions of the target lncRNA.
- **Even Probes**: Bind to the even-numbered regions of the target lncRNA.

This design provides a cross-validation mechanism


### Implementation Steps

1. Perform cross-linking on cell lysates to stabilize RNA-DNA interactions.
2. Hybridize biotinylated even and odd probe sets to the target lncRNA.
3. Use magnetic streptavidin beads (e.g., `PureProteome™ Streptavidin magnetic beads`) to isolate the RNA-probe complexes.
4. Perform stringent washes to remove non-specific interactions.
5. Elute and purify the bound RNA and any associated genomic DNA.
6. Use qPCR or sequencing to analyze the DNA regions enriched by the CHIRP.

<div align=center>
<img src="/imgs/chirp-0.png">
</div>

### Data Analysis

- Significant enrichment over the LacZ control observed via qPCR indicates specific probe-target interactions.
- Peak identification in IGV against the input control confirms the genomic locations of RNA-DNA interactions with high confidence"


<div align=center>
<img src="/imgs/chirp-1.png">
</div>

<div align=center>
<img src="/imgs/chirp-2.png">
</div>

<div align=center>
<img src="/imgs/chirp-3.png">
</div>

### ChIRP-Seq Probe Binding Note

ChIRP-seq employs single-stranded DNA probes targeting specific lncRNAs. Theoretically, these probes might also bind to the lncRNA's encoding genomic DNA due to sequence complementarity. However, since ChIRP-seq lacks a denaturation step for DNA, the actual likelihood of probes binding to double-stranded DNA is very low. This ensures the experiment predominantly focuses on identifying RNA-protein interactions without interference from unintended genomic DNA bindings.

Uses 20-mer biotin-conjugated DNA oligonucleotide, with one probe every 100 base pair of RNA.

The program has a maximum target size of 8 kb, so if the target RNA is greater than 8 kb, the sequence should be broken into segments.

## CPU vs Core vs Thread vs Node

CPU (Central Processing Unit): The CPU is the core part of the computer and is responsible for executing most of the instructions of the computer program. It is the core part of the computer hardware system and is responsible for interpreting and executing instructions from the operating system, applications, and controlling other hardware. A CPU can have one or more cores.

Core: A core is a component of a CPU that can execute instructions independently. A CPU may contain one or more cores (e.g. dual-core, quad-core, octa-core, etc.). Each core can execute a program or part of a program individually, so multiple cores can execute multiple programs or multiple parts of a program at the same time, which can greatly improve computing performance. For example, you mentioned the Intel Xeon Gold 6338, which is a CPU with 32 cores.

Thread: A thread is the smallest unit in program execution and is the basic unit of processor scheduling (task execution). It is contained within the process and is the actual unit of operation within the process. A thread refers to a single sequential flow of control within a process. A process with multiple threads can perform concurrent execution of multiple tasks (functions).

When we talk about a CPU having multi-threading capabilities, we usually mean that each core of the CPU is capable of handling multiple threads of execution simultaneously.

Node: In computing, a node can mean many things, but usually we refer to a device or server running a certain computational task. In a computing network or cluster, a node is usually a single computer that contains its own CPU, memory, and storage devices, and can run an operating system and various software to perform specific computing tasks independently. When multiple nodes join together to form a network or cluster, they can work together to handle larger, more complex computing tasks.

<b>For example, a node with one CPU (32 cores, 64 threads) and 2 A100 GPUs</b>

## Chromatin

Chromatin refers to the complex of DNA, RNA, and proteins that make up the chromosomes in eukaryotic cells. Chromatin is composed of nucleosomes and other proteins that interact with DNA and regulate its function. chromatin can contain RNA, although the amount and type of RNA can vary depending on the specific cell type, developmental stage, and physiological conditions.

There are different types of RNA that can associate with chromatin, including messenger RNA (mRNA), long non-coding RNA (lncRNA), small nucleolar RNA (snoRNA), and others. These RNA molecules can interact with chromatin-associated proteins, including histones and other transcriptional regulators, to modulate gene expression and chromatin structure.

It consists of genomic DNA together with all directly or indirectly associated protein and RNA molecules.

## Chromatin accessibility

Chromatin accessibility refers to the degree to which DNA is accessible to transcriptional machinery, such as RNA polymerase and other regulatory proteins. When the chromatin structure is open and accessible, transcription factors and other regulatory proteins can easily bind to specific DNA sequences and initiate gene expression. Conversely, when the chromatin structure is tightly packed and inaccessible, these proteins may not be able to access the DNA, and gene expression may be suppressed.

<div align=center>
<img src="imgs/chromatin_access.jpg">
</div>

## Chromium X Series - 10x Genomics

A commom strategy for single cell RNA library prep. The figure explains better. Please see below.
<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/1368b03b35308002fcf2890f0f2820a7e16312e4/imgs/10x1.png">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/1368b03b35308002fcf2890f0f2820a7e16312e4/imgs/10x2.png">
</div>

## Chromium X Series - 10x Genomics single cell library construction

A detailed description of the whole library preparation workflow.

<div align=center>
<img src="/imgs/10x-2.png">
</div>

10x Genomics technology is based on the partitioning of samples and reagents into droplets, each called a Gel
Bead in Emulsion (GEM). Once partitioned, the Gel Bead dissolves and its oligo primers are released into the
aqueous environment of the GEM. The cell captured in the GEM is also lysed. The contents of the GEM (oligos,
lysed cell components and Master Mix) are incubated in an RT reaction to generate full-length, barcoded cDNA
from the poly A-tailed mRNA transcripts. The reverse transcription reaction is primed by the barcoded Gel Bead
oligo and the reverse transcriptase incorporates the template switch oligo via a template switching reaction at
the 5’ end of the transcript. The GEMs are then “broken”, pooling single-stranded, barcoded cDNA molecules
from every cell. A bulk PCR-amplification and Enzymatic Fragmentation follow. Size selection is then used to
optimize the insert size of the double-stranded cDNA prior to library construction. During library construction
Read 2 is added by Adapter ligation. Illumina P5 and P7 sequences and sample index sequences are added
during the Sample Index PCR. The final library fragments contain the P5, P7, Read 1 and Read 2 sequences used
in Illumina bridge amplification and sequencing. Additionally, each fragment contains the 10x Barcode, UMI and
cDNA insert sequence used in data analysis.

<div align=center>
<img src="/imgs/10x-3.png">
</div>

What does the final seq looks like?

<div align=center>
<img src="/imgs/10x-40.png">
</div>

## Chromium X Series - Why Single Cell RNA-seq Favors the 3' End

In 10x Genomics single-cell RNA-seq, even with the TSO method that could capture the whole RNA, there's still a 3' end bias. Here's why:

- During library prep, RNA is randomly cut.
- Only pieces with the full Read1 sequence get sequenced.
- That's because only these pieces can be spotted by primers to add the P5 adapter during the library making.
- The Read1 part is on the RNA's 3' end, so that's the part that mostly gets read.

So, in the end, you get more reads from the RNA's 3' end because that's where Read1 is, and that's necessary for the sequencing to work.

<div align=center>
<img src="/imgs/sc1.png">
</div>

<div align=center>
<img src="/imgs/sc2.png">
</div>

<div align=center>
<img src="/imgs/sc3.png">
</div>

<div align=center>
<img src="/imgs/sc4.png">
</div>



## Chromium X Series - 3' vs 5' Libraries

The main difference between Single Cell 3' and Single Cell 5’ Gene Expression libraries lies in which end of the RNA transcript they target. The 3' library captures sequences starting from the RNA's 3' end, utilizing a polydT sequence on a gel bead oligo. Conversely, the 5' library begins sequencing from the TSO end after reverse transcription, using a polydT supplied as an RT primer. This distinction dictates whether the 3' or 5' end of the RNA is represented in the sequencing data.

<div align=center>
<img src="/imgs/sc5.png">
</div>

## Chromium X Series - Nanopore + 10X Sequencing Prep

For the combination of Nanopore and 10X sequencing, the library is prepared without enzymatic fragmentation. This method preserves the full-length cDNA and avoids the introduction of 3'/5' end bias in the sequencing data.

## Chromoplexy and chromothripsis

Chromoplexy:
Chromoplexy is like a complicated puzzle. It takes pieces from different chromosomes (DNA sequences) and mixes them up, creating a new configuration. In terms of sequences, it might look like parts of Chromosome 1, 3, and 5 getting jumbled together.

Chromothripsis:
Chromothripsis is like shattering a glass. It typically affects one or two chromosomes that get broken into many pieces. These pieces are then put back together in a haphazard manner, often missing some fragments (like when you can't find all pieces of the shattered glass).

Their difference:

1. Scope: Chromoplexy usually involves multiple chromosomes, while chromothripsis typically affects one or two chromosomes.
Outcome: Chromoplexy results in a more mixed or shuffled sequence with pieces from different chromosomes. Chromothripsis results in a more shattered and reassembled sequence from the same chromosome.

2. Missing DNA: Chromothripsis often has missing pieces, while chromoplexy usually doesn't.

Difference from normal SVs (Structural Variants):

Normal structural variants are usually simpler. They can involve deletions (missing pieces), duplications (extra copies), inversions (flipped sequences), and translocations (pieces switching places). Chromoplexy and chromothripsis are more complex, involving multiple changes happening all at once.

Think of normal structural variants like moving furniture in your house (a single couch or table), while chromoplexy and chromothripsis are like rearranging or remodeling multiple rooms at once.

## Chromosome

A chromosome is a linear strand of DNA that is compacted and organized by proteins, including histones, into a highly condensed structure. Chromosomes carry genetic information in the form of genes and are passed down from parent to offspring during cell division.

## Cis vs. Trans Regulatory Elements

- **Cis and Trans**: Terms indicating physical position; "cis" means "on the same side," and "trans" means "across from."

- **Cis-acting Elements**: Regulatory DNA sequences (e.g., promoters, enhancers, silencers) located on the same DNA molecule as the gene they regulate, affecting gene expression on the same chromosome. Enhancers, even if physically distant, regulate through folding and are considered cis-acting because they control genes on the same chromosome.

- **Trans-regulatory Factors**: Molecules like transcription factors and lncRNAs that are distinct from the gene they regulate (different molecule types). They can regulate genes located on different chromosomes or distant locations within the same chromosome.

In summary, cis-regulation occurs within the same DNA molecule, affecting genes on the same chromosome, while trans-regulation involves different molecules (proteins or RNAs) acting across DNA molecules or chromosomes.

- **Cis and Trans in Red-C**: In the context of Red-C technology, "cis" refers to the RNA interacting with DNA regions on the same chromosome, including those physically distant but brought into proximity through chromosomal folding. "Trans" describes RNA interactions with DNA regions on different chromosomes. This usage underscores the spatial organization and dynamics within the nucleus, which are crucial for understanding gene regulation and expression patterns.

## Complex structural variation (SV)  
While SV is typically defined by its canonical forms – duplication, deletion, insertion, inversion and translocation – recent breakpoint mapping studies have revealed a surprising number of “complex” variants that evade simple classification. Complex SVs are defined by clustered breakpoints that arose through a single mutation but cannot be explained by one simple end-joining or recombination event.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/e3f4b7bc9721bb073c43a3d64f89c0f0c49fe7e3/imgs/1-s2.0-S1672022921001431-gr1_lrg.jpg">
</div>

## Copy number variation (CNV)

Copy number variation (CNV) is a phenomenon in which sections of the genome are repeated and the number of repeats in the genome varies between individuals. Such regions may or may not contain a gene(s).

## CpG site and CpG islands

1. **CpG Sites**:
   - CpG sites refer to specific DNA sequences where a cytosine (C) is followed by a guanine (G) directly connected through a phosphate bond, denoted as "CpG" (the 'p' indicates the phosphate link).
   - These sites are widely distributed across the entire genome but are generally less frequent in number.
   - In most regions of the genome, CpG sites tend to be methylated.
   - Methylation of CpG sites can influence gene expression, but this effect is usually more localized compared to CpG islands.

2. **CpG Islands**:
   - CpG islands are longer regions in DNA with a high density of CpG sites, typically extending over at least 200 base pairs and characterized by a higher CpG density and G+C content compared to the surrounding regions.
   - They are often located near gene promoters, especially around the transcription start sites (TSS).
   - Unlike CpG sites in other genomic regions, CpG islands are usually unmethylated and associated with active gene expression.
   - The methylation status of CpG islands is closely linked to gene expression regulation, cell-type-specific expression patterns, and the development of diseases, particularly cancer.

Therefore, while both CpG sites and CpG islands involve sequences of cytosine and guanine, they differ significantly in terms of function, distribution, and their impact on gene expression. Research on CpG islands is particularly important in epigenetics and disease studies due to their key role in regulating gene expression.

## CRAM (file format)
Compressed Reference-oriented Alignment Map (CRAM) is a compressed columnar file format for storing biological sequences aligned to a reference sequence. CRAM was designed to be an efficient reference-based alternative to the Sequence Alignment Map (SAM) and Binary Alignment Map (BAM) file formats.

## CRISPR/Cas9 Gene Editing

- **gRNA Equivalence**: gRNA = sgRNA (in lab) = crRNA (in nature) + tracrRNA (in nature).

- **Editing Principle**: The PAM sequence is adjacent to the target DNA sequence. The sgRNA binding sequence is typically 17-20 base pairs long, known as the "seed sequence."

- **Complex Formation and Action**: sgRNA forms a complex with Cas9, targets the specific DNA sequence, unwinds it, and Cas9 cuts the target gene.

- **Gene Repair**: Subsequently, the cell's repair mechanism is leveraged to accomplish gene editing following the Cas9-induced cut.

<div align=center>
<img src="/imgs/cas9.png">
</div>

There are several applications of CRISPR/Cas9 technology, including but not limited to:

- The simplest form, where cutting a gene and then allowing it to repair tends to introduce errors, thereby silencing the gene.

- Inactivating Cas9's cutting domain and then fusing it with a deaminase to achieve precise point mutations.

- Completely inactivating Cas9 and then linking it to a transcription factor (or linking the gRNA to a transcription factor) to enhance gene expression through targeted binding.

- Completely inactivating Cas9 and then attaching it to the KRAB domain, which recruits and physically silences the gene.

- Completely inactivating Cas9 and then attaching it to GFP for structural observation.

- Large fragment deletions require two sets of sgRNA, each guiding Cas9 to cut one end of the deletion segment.

- Cutting open a gene and then introducing a new sequence for insertion can achieve sequence insertion.


## CRISPRi Technology Overview

- **sgRNA (Single-Guide RNA)**: Directs CRISPR-associated proteins (Cas) to specific DNA sequences, targeting the dCas9-MeCP2-KRAB complex to genomic regions like enhancers for gene expression modulation without DNA cleavage.

- **MeCP2-BSD Vector**: Incorporates dCas9 linked with MeCP2 and KRAB repression domains that bind targeted DNA and remodel chromatin into a repressive state, reducing gene expression. It includes a blasticidin resistance gene for cell selection after transduction.

- **dCas9 (Dead Cas9)**: A modified version of Cas9 that binds to DNA without cutting, used in CRISPRi to deliver repressor domains to the target site, effectively silencing gene activity without altering the DNA sequence.

- **CRISPRi Mechanism**: Uses dCas9 fused with transcriptional repressors to inhibit transcription at targeted sites including gene promoters and enhancers, offering a precise method to study gene function and regulatory elements without making permanent changes to the DNA.


## Differential transcript usage (DTU) analysis

A differential transcript usage (DTU) analysis is testing for proportional differences in a gene's transcript composition, and has been of rising interest for many research questions, such as analysis of differential splicing or cell-type identification.

## Directional cloning using restriction enzymes

The use of two restriction enzymes with different recognition sites and incompatible sticky ends is the most efficient and reliable method for directional cloning. 

<div align=center>
<img src="/imgs/direct.png">
</div>

## Discordant reads and cordant reads

1) Concordant reads (Properly aligned reads)
2) Discordant reads (improperly aligned reads: important to identify genome alteration events)

Concordant reads: have span size within the range of expected fragment size and consistent orientation of read pairs with respect to reference.

Discordant reads: have unexpected span size/inconsistent orientation of read pairs.

Example:

Concordant Reads** - R1----->Expected Mapping distance and orientation<-----R2

Discordant reads: Discordant reads have different categories

A) Based on mapping Distance

R1-----> Unexpected mapping distance <-----R2
B) Based on read orientation (Expected read orientation for Paired-end data should be R1 (Forward) R2 (Reverse): FR orientation, but in case of discordant reads, orientations are either FF or RR)

R1-----> R2-----> [FF orientation]

R1<----- R2<---- [RR orientation]

## Double-Strand cDNA

Double-strand cDNA (ds-cDNA) is created by reverse transcription of RNA followed by a second-strand synthesis reaction that uses the single-strand cDNA as a template to form a complementary strand. ds-cDNA is typically used in the preparation of libraries for next-generation sequencing (NGS) applications. The double-stranded nature of the cDNA allows for more stable storage and provides a format that is compatible with various sequencing technologies. The ds-cDNA libraries can be used to perform in-depth analysis of transcriptomes, gene expression profiling, and discovery of novel transcripts or alternative splicing events.

## DNA dA-Tailing
This is a next step of DNA end-repair. Adding a non-template dAMP (dA) to the 3’ end of a blunt-ended DNA fragment. This incorporated 3’-dA prevents concatemer formation and prepares the DNA fragment for subsequent ligation of adaptors or cloning vectors that have complementary 3’-dT overhangs. 

In other words, add "A" base to the 3´ end of a blunt phosphorylated DNA fragment. This treatment creates compatible overhangs for the next step of DNA sample preparation.T

## Dr-Seq single cell sequencing technology

Dr-Seq technology separates the sequencing of gDNA and mRNA from a single cell. The workflow is as follows:

### Cell lysis and reverse transcription
- Cells are lysed and reverse transcription is performed using primer **Ad-1x** to obtain cDNA with cell-specific barcodes  
- The cDNA contains **Ad-1x** sequence on 5' end
- Quasilinear amplification is then performed to amplify cDNA, which introduces **Ad-2** sequence
  - The amplified cDNA can have two forms:  
    1. **3' end with Ad-2**  
    2. **5' end with Ad-1x and 3' end with Ad-2** (due to linear amplification nature)

### Sequencing of gDNA and mRNA
The amplified product is split into two parts:

#### mRNA Sequencing
- In vitro transcription using T7 promoter in **Ad-1x** to specifically amplify cDNA
- Build library and sequence

#### gDNA Sequencing 
- PCR using primers complementary to **Ad-2** to amplify gDNA
- Remove Ad-2 sequence
- Build library and sequence 

### Bioinformatic analysis
- Align reads to reference transcriptome and genome
- Distinguish cDNA and gDNA reads by alignment
  - However, has limitations due to sequence ambiguities

In summary, Dr-Seq lyses single cells, reverse transcribes mRNA, and amplifies both gDNA and cDNA with different barcoded primers. It then separates the two and processes them independently for sequencing.

<div align=center>
<img src="/imgs/DR-seq.png">
</div>

## Enhancer
An enhancer is a short (50–1500 bp) region of DNA that can be bound by proteins (activators) to increase the likelihood that transcription of a particular gene will occur. They can be located up to 1 Mbp (1,000,000 bp) away from the gene, upstream or downstream from the start site. Enhancers are found mostly in the intergenic and intronic regions, while a few enhancers have been found within exons.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a73102a47ec95cf8ca7fbb3f7e938ee279da01dc/imgs/enhancer.png">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a73102a47ec95cf8ca7fbb3f7e938ee279da01dc/imgs/enhancer2.png">
</div>

Here is an enhancer diagram. Within this DNA sequence, protein(s) known as transcription factor(s) bind to the enhancer and increases the activity of the promoter. 1. DNA 2. Enhancer 3. Promoter 4. Gene 5. Transcription Activator Protein 6. Mediator Protein 7. RNA Polymerase

## Enzymatic Methyl-seq (EM-seq)

EM-seq is a bisulfite-free, enzymatic method for detecting 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC) in DNA. This technique is an alternative to bisulfite sequencing, offering a less damaging approach to detect DNA methylation.

### Detailed Workflow
1. **DNA Shearing**: Genomic DNA is fragmented into smaller pieces.
2. **End Repair/A-Tailing**: The DNA fragments are polished and adenine-tailed, preparing them for ligation.
3. **Adaptor Ligation**: Adaptors are ligated to the DNA fragments, preparing them for amplification and sequencing.
4. **Oxidation and Protection**: An oxidation step protects 5mC and 5hmC from subsequent deamination, thus preserving the methylation information.
5. **Deamination**: Unprotected cytosines are deaminated to uracil (U), differentiating them from methylated cytosines.
6. **PCR Amplification**: The NEBNext® Q5U Master Mix amplifies the DNA, enriching for fragments that have undergone the conversion.
7. **Sequencing**: High-throughput sequencing, like Illumina sequencing, is used to read the DNA sequences and identify methylation patterns.

### Advantages
- Preserves DNA integrity by avoiding harsh chemical treatments.
- Enables detection of both 5mC and 5hmC.
- Compatible with low input DNA and single-cell analysis.

The EM-seq approach provides greater accuracy and less DNA damage, enhancing the reliability of methylation data for epigenetic studies.

## Enhancer

- **Definition**: Enhancers are cis-regulatory DNA elements that upregulate gene transcription, with no fixed nomenclature.
- **Function**: They can be located far from the gene they regulate and work by looping to interact with the gene's promoter.
- **Variability**: Their sequences and locations vary between cell types, requiring specific identification and validation.
- **Study**: Their activity is studied through techniques like ChIP-seq and Hi-C.

### Characteristics of Enhancer Regions

- Enhancers are identified by specific histone modifications, notably H3K4me1 (monomethylation of histone H3 at lysine 4) signifying enhancer presence, and H3K27ac (acetylation of histone H3 at lysine 27) indicating active enhancers.
- Proteins like Cohesin, Mediator, and CTCF bind to enhancers, playing a crucial role in the formation of chromatin loops that bring enhancers in proximity to promoters for gene regulation.
- The binding of key transcription factors to enhancers dictates their regulatory activity.
- Enhanced chromatin accessibility is a hallmark of enhancers, reflecting a more open chromatin state that facilitates transcriptional machinery access.

## Enhancer-Promoter Loops

Enhancer-promoter loops play a pivotal role in the regulation of gene expression in eukaryotic cells. These structures are crucial for facilitating the appropriate levels of gene transcription and are key to understanding the complexities of gene regulation.

### Understanding Enhancers and Promoters

- **Enhancers**: DNA sequences that may be located far from the genes they regulate. They contain specific binding sites for transcription factors and are essential in upregulating the expression of genes.

- **Promoters**: Sequences located at the beginning of genes, serving as the primary binding sites for RNA polymerase and other regulatory proteins for the initiation of transcription.

### Sequence of Interactions in Loop Formation

1. **Initial Binding to Enhancers**: Transcription factors typically bind to enhancers first. These enhancer-bound factors recruit additional proteins, creating a conducive environment for transcriptional activation.

2. **Chromatin Remodeling and Looping**: The chromatin structure is remodeled, facilitating the physical proximity of the enhancer to the promoter region. This looping is mediated by various proteins and changes in the chromatin landscape.

3. **Subsequent Promoter Interaction**: Following the formation of the loop, the enhancer, now in close spatial relation to the promoter, influences the activity at the promoter. This can involve direct interactions between the enhancer-bound transcription factors and the promoter region or indirect effects mediated by chromatin structure.

### Functional Significance

- **Gene Expression Regulation**: These loops enable enhancers to influence gene transcription effectively, even when located at a significant distance from the target gene.

- **Cell Type Specificity**: The pattern of enhancer-promoter loops varies across different cell types, contributing to cell-specific gene expression profiles.

- **Implications in Development and Disease**: Abnormalities in these interactions can lead to misregulation of gene expression, with implications for various diseases, including cancer.

### Research Techniques

Methods like Chromosome Conformation Capture (3C) and its derivatives, along with ChIP-seq, are employed to study enhancer-promoter loops, offering insights into the three-dimensional organization of chromatin in gene regulation.

## Enhancer RNAs (eRNAs)

Enhancer RNAs (eRNAs) are non-coding RNA molecules transcribed from enhancer regions, playing a key role in gene regulation. They enhance gene expression by:

- **Modifying Chromatin Structure**: eRNAs facilitate changes in chromatin that bring enhancers and gene promoters into closer spatial proximity.
- **Recruiting Transcription Factors**: By attracting transcription factors to enhancers and promoters, eRNAs increase gene transcription activity.
- **Bridging Enhancers and Promoters**: eRNAs may physically bridge enhancer-promoter interactions, stabilizing the transcriptional machinery.

**Action Scope**:
- Primarily, eRNAs act locally, boosting the activity of the enhancer from which they are transcribed and its target gene.
- Some eRNAs might also have distal effects, potentially affecting the activity of other genes or enhancers indirectly through nuclear mobility or chromatin remodeling.

## Epigenome

An epigenome consists of a record of the chemical changes to the DNA and histone proteins of an organism; these changes can be passed down to an organism's offspring via transgenerational stranded epigenetic inheritance. Changes to the epigenome can result in changes to the structure of chromatin and changes to the function of the genome.

The epigenome refers to a collection of chemical modifications that occur on the DNA molecule and its associated proteins, which can alter gene expression without changing the underlying DNA sequence. These modifications include the addition or removal of small chemical tags, such as methyl groups, to the DNA molecule or its associated histone proteins, as well as the presence of non-coding RNA molecules. These modifications can influence the way that the DNA is packaged and accessed within the cell, allowing genes to be turned on or off in response to environmental and developmental cues.

The epigenome can be inherited from one generation to the next. In some cases, the modifications to the epigenome that occur during an individual's lifetime can be passed on to their offspring through a process called epigenetic inheritance. This occurs when the modifications to the epigenome are present in the sperm or egg cells that will form the embryo.

ATAC-seq and ChIP-seq are two methods used to study the epigenome.

## Epitope

It is capable of stimulating an immune response. This is usually one to six monosaccharides or five to eight amino acid residues on the surface of the antigen. Each antigen typically has many epitopes. 

## Exitron

Exitrons are exon-like introns located within protein-coding exons. Removal or retention of exitrons through alternative splicing increases proteome complexity and thus adds to phenotypic diversity.Exitrons are defined as introns within protein-coding exons that, when retained, maintain the protein-coding potential of the transcript. Marquez and colleagues argue that four features distinguish exitrons from other introns: high GC content, absence of stop codons, overrepresentation of a size corresponding to multiples of three nucleotides, and prevalence of synonymous substitutions (as usually observed for exonic sequences).

Transcripts with exitrons in their sequences can be distinguished from those with retained introns in three ways. First, transcripts containing exitrons are transported out of the nucleus to be translated, whereas those containing introns are identified as incompletely processed and are kept in the nucleus where they cannot be translated. Second, only transcripts with exitrons of lengths not divisible by three have the potential to incorporate premature termination sequences, while sequences with introns normally result in premature termination. Third, exitron transcripts are usually the major isoform, but those with introns are only present in small amounts.

## Extrachromosomal circular DNA （ecDNAs）
Extrachromosomal circular DNA (eccDNA) is a type of double-stranded circular DNA structure that is derived and free from chromosomes. It has a strong heterogeneity in sequence, length, and origin and has been identified in both normal and cancer cells. In contrast to previously identified circular DNA structures (e.g., bacterial plasmids, mitochondrial DNA, circular bacterial chromosomes, or chloroplast DNA), eccDNA are circular DNA found in the eukaryotic nuclei of plant and animal (including human) cells. Extrachromosomal circular DNA is derived from chromosomal DNA, can range in size from 50 base pairs to several mega-base pairs in length, and can encode regulatory elements and full-length genes. EcDNA, as the vehicles for oncogene and drug-resistance genes, enables them to be rapidly amplified, and lead to overexpression consequently.For instance, oncogenes EGFR and c-MYC were found in ecDNA and amplified in human cancer tissues than normal tissues. Lacking a high-order chromatin structure, suppressing histone modifications, and insulator shackle make ecDNAs more accessible than their genome counterparts which may facilitate promoter–enhancer interactions, transcription initialization, and achieving additional expression.

## F1-score
The F1-score combines the precision and recall of a classifier into a single metric by taking their harmonic mean. It is primarily used to compare the performance of two classifiers. Suppose that classifier A has a higher recall, and classifier B has higher precision. In this case, the F1-scores for both the classifiers can be used to determine which one produces better results.

## FISH Technology

Fluorescence in situ hybridization (FISH) is a cytogenetic technique that uses fluorescent probes to detect and localize specific DNA or RNA sequences within cells.

### Steps in FISH:
1. **Probe Design**: Fluorescent probes target specific nucleic acid sequences.
2. **Hybridization**: Probes bind to their targets within fixed cells.
3. **Imaging**: Fluorescence microscopy reveals the location of the probes, indicating target sequence distribution.

### Example: SNHG6 in Osteosarcoma Cells
- **Cy3 Channel**: Shows SNHG6 lncRNA signals in red, localized mainly in the cytoplasm.
- **DAPI Channel**: Stains nuclei in blue to provide a reference point.
- **Merge**: Combines Cy3 and DAPI images, illustrating the spatial relationship between SNHG6 lncRNA and the nucleus. 

This visualization highlights SNHG6's cytoplasmic distribution, suggesting its role in cytoplasmic functions rather than nuclear processes.

<div align=center>
<img src="/imgs/fish.jpg">
</div>

## Flongle adaptor and flow cell
Flongle is an adapter for MinION or GridION that enables direct, real-time DNA or RNA sequencing on smaller, single-use flow cells. It delivers up to 1.8 Gb of data. Providing immediate access to sequence data, Flongle is designed to be the quickest, most accessible and cost-efficient sequencing system for smaller or more frequently performed tests and experiments. Flongle flow cells should be used within 4 weeks of receipt.

## Flow Cytometry for Apoptosis

In flow cytometry, live cells do not take up dyes like Propidium Iodide (PI) and are negative for Annexin V, indicating intact cell membranes. This technique clearly distinguishes live cells from those in early or late stages of apoptosis or necrosis, making it essential for cell health assessments.

<div align=center>
<img src="/imgs/snnexin.png">
</div>

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

## Germline mutation and somatic mutation
A genetic change in a germ cell (egg or sperm) that becomes the DNA of each cell in the offspring's body. A variant (or mutation) contained in the germline can be passed from parent to offspring and is therefore inherited. They are also called germline mutations. INDELs are identified by removing the germline mutations from INDELs.

<b>Example</b>

Every person is made up of lots of tiny things, which we call cells, just like little factories. Each little factory has a set of instructions that tells it how to work. These instructions actually consist of two parts, like you have a storybook, half of the stories come from your dad, and the other half come from your mom. This is what we call "diploid".

Now, there are two situations that could make these instructions wrong:

The first situation is when you just begin as a new life, that is, when your dad's sperm and mom's egg combine, you get your own instructions. These instructions are copied and assembled from your dad and mom's instructions. But sometimes, there may be errors in this copying process, just like when printing a new book, some words get printed wrong. This kind of error that appears at the very beginning of your life is called a germline mutation. Because this error appeared at the very beginning of your life, the instructions in every little factory (each cell) in your body will have this error.

The second situation is when you start to grow up, your body needs more little factories. These new little factories are replicated from existing little factories. This is like you already have a book, but you need to copy more of this book. But in the process of photocopying, sometimes there will be errors, just like when the photocopier is copying, some words get copied wrong. This kind of error that appears after your life has been going for a while, in the process of little factories replicating, is called a somatic mutation. This kind of error only appears in a part of the little factories because they occur during the replication of the little factories.

So, simply put, a germline mutation is an error that existed when your life just began, and a somatic mutation is an error that occurred after your life has been going for a while, in the process of little factories replicating.

Let's introduce the concept of VAF (Variant Allele Frequency). This concept is similar to wanting to find out the proportion of wrong words in the entire book in the instructions.

Suppose we find an error in the instructions (genome), we can see the proportion of this error in all little factories (cells). This is what we call VAF.

For germline mutations, because this error appeared at the very beginning of your life, the instructions in every little factory (each cell) in your body will have this error. If this error appears in all instructions (two haplotypes), then its VAF is 100%. If this error only appears in half of the instructions (one haplotype), then its VAF is 50%. If this error does not appear in any instructions, then its VAF is 0%.

For somatic mutations, because this error appeared after your life has been going for a while, in the process of a part of the little factories (cells) replicating, this error will only appear in a part of the little factories' instructions. Therefore, its VAF could be any value, depending on the proportion of this error appearing in all little factories.

So, simply put, the VAF of germline mutations is usually 0%, 50% or 100%, while the VAF of somatic mutations could be any value, because they only appear in a part of the little factories (cells).

## GRID-seq (All to all RNA-DNA interaction)

GRID-seq is an innovative technique for mapping genome-wide RNA-chromatin interactions, leveraging a unique bridge linker and the special cutting properties of the MmeI enzyme.

### Key Components

- **Bridge Linker**: Combines a double-stranded DNA (dsDNA) segment for DNA attachment and a single-stranded RNA (ssRNA) overhang for RNA ligation.
- **MmeI Enzyme**: Cuts about 20 nucleotides **upstream** from its recognition site, enabling precise generation of 85 bp DNA fragments.

<div align=center>
<img src="imgs/grid.png">
</div>

## GridION Mk1
The GridION Mk1 provides users with five sequencing ports where MinION flow cells or Flongle adapters with flow cells can be connected, as well as a high performance integrated computer and basecalling accelerator. The device can basecall, in real-time, the data generated by five flow cells/Flongles. The current chemistry and software enables generation of up to 150 Gbases of data during a GridION Mk1 run. Up to 250 Gb (all 5 flow cells sequencing).

## HiChIP-seq

### Overview
- **Immunoprecipitation Required**: HiChIP requires immunoprecipitation (IP)
- **End Ligation**: Involves blunt-end ligation of DNA fragments.
- **Purpose**: Detects all-to-all DNA-DNA interactions for specific proteins.

### Loop calling (determine enhancer and promoter) Using H3K27ac Antibody
- **Marker for Active Enhancers/Promoters**: H3K27ac antibody marks active enhancers and promoters, facilitating loop calling.
- **Enrichment**: HiChIP enriches for interactions between regions marked with H3K27ac, primarily enhancers and promoters.
- **Anchors Classification**:
  1. **Promoter Anchor (P)**: Contains a transcription start site (TSS), with or without H3K27ac mark.
  2. **Enhancer Anchor (E)**: Contains H3K27ac mark but no TSS.
  3. **Non-enhancer/Non-promoter Anchor (X)**: Lacks both H3K27ac mark and TSS.
- **Loop Types**:
  1. **Enhancer-Enhancer (E-E)**
  2. **Enhancer-Promoter (E-P)**
  3. **Promoter-Promoter (P-P)**
  4. **Enhancer-Non-regulatory (E-X)**
  5. **Promoter-Non-regulatory (P-X)**
- **Loop Calling Principle**: Based on whether anchors contain H3K27ac and TSS, loops are classified by combinations of different anchor types.

### Loop calling Using Transcription Factor Antibodies
- **TFBS Identification**: To perform loop calling with transcription factor antibodies, knowledge of transcription factor binding sites (TFBS) is crucial.
- **Anchor Classification with TF Antibodies**:
  - **TFBS + TSS**: Promoter anchor.
  - **TFBS + No TSS**: Enhancer anchor.
  - **No TFBS**: Non-regulatory anchor.
- **Loop Classification**: Similar to H3K27ac-based classification, but based on the presence of TFBS and TSS in anchors.
- **Obtaining TFBS Information**: Through ChIP-seq or bioinformatics motif prediction.

<div align=center>
<img src="imgs/hichip.png">
</div>


## Hi-C Technology Overview (All to all DNA-DNA interaction)

### Workflow Steps:

1. **Cross-Linking**: Cells are treated with formaldehyde, cross-linking proteins to DNA, capturing the chromatin interactions.
2. **Digestion**: Chromatin is digested with a restriction enzyme that cuts DNA to create 5’ overhangs.
3. **Biotinylation**: The 5’ overhangs are filled in with biotinylated nucleotides, allowing for the blunt-ended fragments to be ligated together.
4. **Ligation**: DNA ends are ligated, favoring joined ends that are spatially close in the nucleus.
5. **Purification**: Ligation products with biotin are purified using streptavidin beads.
6. **Library Preparation**: DNA is sheared, size-selected, and prepared into a sequencing library.
7. **Sequencing**: High-throughput sequencing is used to read the ligation junctions, revealing which DNA regions are interacting.

<div align=center>
<img src="imgs/hiC_prot.png">
</div>

## Histone proteins

Histone proteins are a family of small, positively charged proteins that are rich in arginine and lysine residues. They are the major protein component of nucleosomes and are involved in the packaging of DNA into chromatin. There are five main types of histones: H1, H2A, H2B, H3, and H4.

<div align=center>
<img src="imgs/Histone.jpg">
</div>

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

Here is a concise Markdown note on the key points of the iMARGI technology:

## iMARGI (All-to-all RNA-DNA interactions, similar to CHAR-seq)

### Overview

- Converts RNA-DNA interactions into unique DNA sequences
- Sequencing reads mapped to infer original RNA-DNA interaction sites

### Key Steps

- Designed linker sequence ligated between RNA and DNA fragments
- Linker contains restriction site for later linearization
- cDNA synthesized from linker-ligated RNA-DNA chimera 
- Bottom DNA strand released, circularized 
- Circularized DNA cut at restriction site to linearize
- Linearized DNA now contains linker halves at 5' and 3' ends
- These linkers contain primer sequences for paired-end sequencing

### Only ssDNA

- Only single round of reverse transcription, as both circularization and linearization use single-stranded DNA
- Restriction site in linker enables controlled linearization

<div align=center>
<img src="/imgs/ima.png ">
</div>

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

## Inverse PCR (I don't like this name, doesn't make any sense)

Inverse PCR is a molecular technique used to amplify unknown sequences adjacent to known regions of DNA. It involves cutting the target DNA, self-ligating into a circular template, and then using outward-facing primers to perform PCR amplification. This allows for the extension from known sequences to identify their neighboring unknown regions.

<div align=center>
<img src="/imgs/inverse.jpg">
</div>

## Lambda control DNA  

Lambda DNA is 48 kb long and serves as a good model system to evaluate the sequencing workflow. The DNA molecule of 48502 basepairs is linear and except for the extreme ends double-stranded. At each end the 5' strand overhangs the 3' strand by 12 bases. The sequences of the ends are complementary.

## Lentivirus transfection

Lentivirus transfection is a way to put genes into cells. It works like this:

1. **Prepare** three plasmids: one with your gene, and two others (PAX2 and VSVG) to help make the virus.
2. **Co-transfect** these plasmids into 293T cells to produce lentivirus particles.
3. **Collect** the virus from the culture medium after a few days.
4. **Infect** your target cells with this virus.
5. The virus enters the cells, and the gene you sent is added to the cells' DNA.
6. **Express** the gene in the target cells to either turn off or increase the activity of certain genes.

<div align=center>
<img src="imgs/lenti.png">
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

## Lipofectamine

Transfection reagents developed by Thermo Fisher Scientific, delivering nucleic acids (such as DNA, siRNA, and miRNA) into mammalian cells, facilitating gene expression studies and RNA interference experiments. Lipofectamine RNAiMAX for siRNA.

## LNCaP cells
LNCaP cells are a cell line of human cells commonly used in the field of oncology. LNCaP cells are androgen-sensitive human prostate adenocarcinoma cells derived from the left supraclavicular lymph node metastasis from a 50-year-old caucasian male in 1977. They are adherent epithelial cells growing in aggregates and as single cells.

## (MALBAC) Multiple Annealing and Looping Based Amplification Cycles

A quasilinear whole genome amplification method. Unlike conventional DNA amplification methods that are non-linear or exponential (in each cycle, DNA copied can serve as template for subsequent cycles), MALBAC utilizes special primers that allow amplicons to have complementary ends and therefore to loop, preventing DNA from being copied exponentially. This results in amplification of only the original genomic DNA and therefore reduces amplification bias.The key to MALBAC is to not make copies of copies, and instead make copies only of the original genomic DNA by protecting the amplification products. MALBAC uses a thermally stable DNA polymerase with strand displacement activity. This important preamplification stage is then followed by exponential amplification of the full amplicons by PCR, generating the amount of DNA required for next-generation sequencing. MALBAC’s sequence-dependent bias is reproducible along the genome from cell to cell. Therefore, signal normalization for CNV noise reduction can be carried out. As shown below, after signal normalization with a reference cell, MALBAC offers the best CNV accuracy.

<div align=center>
<img src="/imgs/malbac2.png">
</div>

## MAS-seq from Pacbio

The MAS-Seq for 10x Single Cell 3’ kit takes as input single-cell cDNA and outputs a sequencing-ready library that results in a 16-fold throughput increase compared to regular single-cell Iso-Seq® libraries.

However, MAS-Seq cannot detect RNA modifications.

Regular 10x single-cell RNA sequencing does not directly detect RNA modifications. The standard protocol for 10x single-cell RNA sequencing uses oligo-dT primers to capture the poly-A tails of mRNA molecules, which are then reverse transcribed to generate cDNA libraries for sequencing. During this process, RNA modifications are not specifically labeled or detected, and therefore are not directly captured in the resulting sequencing data.

Steps include: (1) PCR with biotin-primer and TSO primer, to label full length cDNA; (2) Beads pull down to remove TSO artifact; (3) another PCR to amplify; (4) Add adaptor.

<div align=center>
<img src="/imgs/mas-seq.png">
</div>

## MDA - inverted chimera

The displaced 3'-terminus would be free to reanneal, preferentially at randomly occurring complementary segments on nearby 5'-strands. The outcome will be the joining of two sequences in inverted orientation. The finding that 85% of chimeras had this inverted form supports this as the likely mechanism. The chimeric junctions also reveal the site where the displaced 3'-end annealed to the second template and continued elongation. In the example , from one of the sequenced chimeric junctions, priming was initiated on the new template where the sequence on the 3'-end had annealed to the sequence on the 5'-strand. In 76.8% of cases, the 3'-ends initiated priming at complimentary sequences of ≥ 2 bp on the new templates. These ranged up to 21 base pairs of complimentarity (Fig 4). The complimentarity occurred in 93.5% and 90.9% of cases for inverted and directly joined segments, respectively, when <10 kb apart, but only 34.8% and 27.6%, respectively, when >10 kb apart.

Inverted chimera is characterized by the fact that the original sequence (A), after a small sequence (green, usually 3-31bp), the new sequence is an inverted sequence (B) at another position. This inverted sequence is compared to the negative strand of the genome. The sequence that corresponds to the positive strand is the reverse complementary sequence of B, and a small segment of sequence downstream of this new reverse complementary sequence in plus strand (towards the 3 end of the genome) should be reverse complementary to a segment of the original sequence (A). The segment is taken at breakingpoint forward (toward the 5 end of the genome).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/7bfec898f0dbdf331add143ac431db4384ae59d3/imgs/inveerted_chimera.png">
</div>

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/7bfec898f0dbdf331add143ac431db4384ae59d3/imgs/inverted_chimera.png">
</div>

## Melanoma

Melanoma, the most serious type of skin cancer, develops in the cells (melanocytes) that produce melanin — the pigment that gives your skin its color.

## MHC class I

MHC class I molecules are expressed in all nucleated cells and also in platelets—in essence all cells but red blood cells. It presents epitopes to killer T cells, also called cytotoxic T lymphocytes (CTLs). A CTL expresses CD8 receptors, in addition to T-cell receptors (TCR)s. When a CTL's CD8 receptor docks to a MHC class I molecule, if the CTL's TCR fits the epitope within the MHC class I molecule, the CTL triggers the cell to undergo programmed cell death by apoptosis. Thus, MHC class I helps mediate cellular immunity, a primary means to address intracellular pathogens, such as viruses and some bacteria, including bacterial L forms, bacterial genus Mycoplasma, and bacterial genus Rickettsia. In humans, MHC class I comprises HLA-A, HLA-B, and HLA-C molecules.

## Microhomology in structural variations

Besides the micro-insertions, short stretches (usually 2–6 nucleotides) of identical sequence, known as overlapping microhomology31, were frequently spotted at SV breakpoint junctions.

<div align=center>
<img src="imgs/micro.png">
</div>

## Micro-insertion in structural variantions

Micro-insertions refer to the addition of small segments of DNA (usually less than 500 bp in length) at the breakpoints of structural variations (SVs) in the genome. These inserted sequences are often non-aligned, meaning that they do not match up with the reference genome, and are frequently novel with no homology to known sequences. Micro-insertions are thought to be the result of non-templated DNA synthesis at the rearrangement junctions during the repair of DNA double-strand breaks through the non-homologous end-joining (NHEJ) mechanism.

## MinION

Similar to MinION Mk1C, but come without screen. Starter Packs from $1,000 including 1 flow cell andd 1 kit. The enhanced has 4 flow cells and 1 kit with $3,250.

## MinION Mk1C

Nanopores read the length of DNA or RNA presented to them — from short to ultra-long (longest >4 Mb). Up to 50 Gb per MinION Flow Cell / 2.8 Gb per Flongle Flow Cell. Theoretical max output when system is run for 72 hours (or 16 hours for Flongle) at 420 bases / second. Outputs may vary according to library type, run conditions, etc. Starter Packs from $4,900, including consumables Compatible with Flongle Flow Cells for smaller tests and analyses. The starter pack includes six flow cells and 1 library kit. It is for filed-use and come with a screen.

## Multiple displacement amplification (MDA)  

MDA offers much higher genome coverage than DOP-PCR. However, like DOP-PCR, MDA is an exponential amplification process. This results in sequence-dependent bias, causing overamplification in certain genomic regions and underamplification in other regions. However, such uneven bias of MDA is not reproducible along the genome from cell to cell, causing CNV measurements noisy and normalization ineffective. Nevertheless, MDA has been widely applied since its invention. <b> It is PCR-free, and under isothermal conditions, The DNA fragments are 50–100 kb long. </b>

## Multiplex sequencing 

Allows large numbers of libraries to be pooled and sequenced simultaneously during a single run on Illumina instruments. Sample multiplexing is useful when targeting specific genomic regions or working with smaller genomes. Pooling samples exponentially increases the number of samples analyzed in a single run, without drastically increasing cost or time.

With multiplex sequencing, individual "barcode" sequences are added to each DNA fragment during next-generation sequencing (NGS) library preparation so that each read can be identified and sorted before the final data analysis. These barcodes, or index adapters, can follow one of two major indexing strategies depending on your library prep kit and application.

## Nanopore sequencing
Using nanopore sequencing, a single molecule of DNA or RNA can be sequenced without the need for PCR amplification or chemical labeling of the sample. A MinION flow cell contains 512 channels with 4 nanopores in each channel, for a total of 2,048 nanopores used to sequence DNA or RNA. The wells are inserted into an electrically resistant polymer membrane supported by an array of microscaffolds connected to a sensor chip. Each channel associates with a separate electrode in the sensor chip and is controlled and measured individually by the application-specific integration circuit (ASIC). Ionic current passes through the nanopore because a constant voltage is applied across the membrane, where the trans side is positively charged. Under the control of a motor protein, a double-stranded DNA (dsDNA) molecule (or an RNA–DNA hybrid duplex) is first unwound, then single-stranded DNA or RNA with negative charge is ratcheted through the nanopore, driven by the voltage. As nucleotides pass through the nanopore, a characteristic current change is measured and is used to determine the corresponding nucleotide type at ~450 bases per s (R9.4 nanopore).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/b0129ac21a844b3a0252b247a2dc1086a95bcc8a/imgs/nanopore.png">
</div>

## Nanopore sequencing - duplex basecalling  

Oxford Nanopore announced a new method, “Duplex”  which enables nanopore devices to sequence a template and complement strand of a single molecule of DNA, in succession, in order to achieve very high accuracy sequencing results.  This system can currently achieve modal accuracy trending towards Q30 for a single double stranded molecule of DNA.  The goal of Duplex experiments is to achieve a high proportion of instances where the complement follows the template strand through the nanopore and gives the system ‘two looks at the same sequence pairs’.

## Nanopore sequencing - Targeted sequencing - enrichment strategy - amplicon

## Nanopore sequencing - Targeted sequencing - enrichment strategy - hybrid-capture

## Nanopore sequencing - Targeted sequencing - enrichment strategy - CRISPR/Cas9

## Nanopore sequencing - Targeted sequencing - enrichment strategy - adaptive sampling

## Nanopore sequencing - whole genome DNA sequencing kits  

A wide range of library preparation kits are available to suit all whole genome sequencing requirements. Amplification-free kits allow direct, long-read sequencing of native DNA, eliminating the potential for PCR bias and allowing the detection of base modifications alongside nucleotide sequence. Amplification-based kits are also available, enabling whole genome sequencing from low input amounts or poor quality DNA (e.g. FFPE).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/8611e16d9176ca1a927e9020fd63350694b6c812/imgs/WGS.png">
</div>

## Nanopore sequencing - whole genome DNA sequencing kits - Ligation Sequencing Kit - SQK-LSK109

An old version of whole genome DNA sequencing kits. Recommend using our latest version: SQK-LSK114.

## Nanopore sequencing - whole genome DNA sequencing kits - Ligation Sequencing Kit - SQK-LSK110  

Processing singleplex samples. $599.00. Utilise upstream processes such as size selection or whole genome amplification. ~60 minutes. Input requirement	1000 ng gDNA. No PCR. Shipped at 2–8°C. Long-term storage -20°C. The library preparation method is straightforward: DNA ends are repaired and dA-tailed using the NEBNext End Repair/dA-tailing module, and then sequencing adapters, supplied in the kit, are ligated onto the prepared ends.

<b>The kit is optimised to generate maximum output and read length. The kit contains an updated sequencing adapter (Adapter Mix F, or AMX-F), which is designed to turn over significantly less fuel during a sequencing run compared to previous kit versions.</b> If your experiment requires long reads, it is recommended to start with full-length gDNA, and fragmentation/shearing is neither advised nor required. To determine purity, we suggest using the Nanodrop to measure the A260/280 and A260/230 ratios and we recommend that the sample should meet the following criteria: A260/280 = 1.8, A260/230 = 2.0-2.2. The Ligation Sequencing Kit is compatible with upstream whole genome amplification (for applications where under 1 ng of sample is available).

## Nanopore sequencing - whole genome DNA sequencing kits - Ligation Sequencing Kit - SQK-LSK112  

Legacy product, no longer support.

## Nanopore sequencing - whole genome DNA sequencing kits - Ligation Sequencing Kit - SQK-LSK114  

The kit is optimised to achieve sequencing accuracies of over 99% (Q20+) with high output on our latest nanopore: R10.4.1. Our flow cell priming and sequencing reagents have been reformulated to be compatible with this improved Kit 14 adapter and R10.4.1 nanopore.Kit 14 also includes previous updates such as the higher capture rate of DNA to enable lower flow cell loading amounts, and fuel fix technology, allowing users to run longer experiments without the need for fuel addition during the run. Users can either start with 1 μg of gDNA or 100-200 fmol of shorter-fragment input such as amplicons or cDNA. If your experiment requires long reads, it is recommended to start with full-length gDNA, and fragmentation/shearing is not advised. 

## Nanopore sequencing - whole genome DNA sequencing kits - Rapid Sequencing Kit - SQK-RAD004  

generates sequencing libraries from extracted gDNA in 10 minutes using a simple two-step protocol. At the heart of the kit is a transposase which simultaneously cleaves template molecules and attaches tags to the cleaved ends. Rapid Sequencing Adapters are then added to the tagged ends.  Input requirement 400 ng gDNA (<30 kb). No PCR. Shipped at 2–8°C. Long-term storage -20°C.  Due to the simple nature of the workflow and the fact that little sample manipulation is required (e.g. minimal pipetting steps and no clean-ups) some very long reads can be achieved with this kit, despite the required transposase fragmentation. However, in order for long reads to be observed in sequencing, long fragments need to be present in the sample in the first place. Typical throughput is lower than Ligation Sequencing Kit.

## Nascent RNA Characteristics

Nascent RNA, freshly transcribed within the cell. 

- **Composition**: Includes enhancer RNAs (eRNAs) and can lack the polyadenylation (poly(A) tail) typical of mature mRNAs.
- **Regulation**: Capable of both cis-regulation, affecting nearby genes at the site of synthesis, and trans-regulation, influencing distant genomic loci.
- **Sequence Features**: Often capped but not always polyadenylated, nascent RNAs have unique sequence characteristics that facilitate their regulatory functions.

## NEBNext Quick Ligation Module (E6056)

Another name, NEBNext Quick T4 DNA Ligase. This module is also compatible with some Oxford Nanopore Technologies workflows. 

## Nucleosome

 A nucleosome is the basic unit of chromatin structure and consists of a segment of DNA wrapped around a core of eight histone proteins (two copies of each of the histones H2A, H2B, H3, and H4). Nucleosomes compact DNA and help regulate access to the DNA by various proteins.

 <div align=center>
<img src="imgs/Chromosomes-are-made-of-DNA-histone-protein-complexes-Chromosomal-DNA-is-packaged.png">
</div>

## Nucleosome-Depleted Regions (NDRs) and Cis-Regulatory Elements

Open chromation = Nucleosome-Depleted Regions (NDRs) 

The identification of open chromatin, i.e., Nucleosome-Depleted Regions (NDRs), is indeed a critical early step in the process of distinguishing cell-type-specific NDRs and cis-regulatory elements. 

The analysis includes:

### Identification of Open Chromatin
Initially, researchers use a sliding window chi-square test to identify NDRs through the analysis of single-cell-level methylation data. This step involves searching for regions that show a methylation pattern at GCH sites, areas that have a higher level of GCH than the genomic background and meet specific criteria for size and coverage of GCH sites.

### Identification of Cell-Type-Specific NDRs
Once open chromatin regions have been identified, researchers look for specific NDRs within each major cell type group (such as spermatogonia, spermatocytes, and spermatozoa). These are regions that display a methylation level at least 10% higher than the maximum level found in other cell types within the group.

### Identification of Shared NDRs Across Consecutive Stages
In addition, NDRs shared between consecutive stages are sought after. These NDRs display a methylation level difference of less than 10% across consecutive stages but more than 10% compared to the remaining stages within each group.

### Correlation Calculation
For the identified specific NDRs, researchers calculate the correlation between their methylation levels and the RNA expression levels of their neighboring genes. "Neighboring genes" refer to the closest gene within 2 kb of the TSS (proximal NDRs) or within 100 kb of the TSS (distal NDRs). Only NDRs with a correlation coefficient greater than 0.6 are defined as cis-regulatory elements.

The process is iterative and multi-staged, potentially involving extensive data analysis and statistical validation. Through this method, researchers can pinpoint those NDRs and cis-regulatory elements that are crucial for gene expression regulation in specific cell types or developmental stages.

 <div align=center>
<img src="imgs/ndr1.png">
</div>

So the NDRs for the red color of these three genes are NDRs that are shared between the three periods, except that the chromatin level changes and gets lower and lower.

The red NDRs (Nucleosome-Depleted Regions) of the three genes (UTF1, MEIOB, TNP1) shown in the figure are regions shared by these genes in sperm cells at different developmental stages. These areas are marked in red because they are open chromatin regions at these stages, indicating their crucial role in gene expression regulation during spermatogenesis.

The changes in chromatin accessibility shown in the figure suggest that as development progresses, the chromatin accessibility near these genes decreases, which might be related to a reduction in gene expression levels. For instance, for the MEIOB gene, the chromatin is more accessible in the early stages, which may be related to its function in the early stages of spermatogenesis. As development progresses, the chromatin accessibility decreases, possibly reflecting that this gene no longer needs to be actively expressed in later developmental stages.

## Nucleosome-Depleted Regions (NDRs) and transcription factor motif

<div align=center>
<img src="imgs/ndr2.png">
</div>

### What is a Motif?

A motif in genetics refers to a specific, short sequence of nucleotides within DNA that is recognized and bound by a particular transcription factor. These motifs are crucial for gene expression regulation, as they direct the binding of transcription factors to specific DNA locations, impacting the gene's transcription process.

### Insights from Figure 3e

Figure 3e focuses on transcription factor motif enrichment in distal Nucleosome-Depleted Regions (NDRs).

- **Transcription Factor Binding**: This figure analyzes which transcription factors' binding motifs are enriched in distal NDRs. The presence of specific motifs indicates potential activity of certain transcription factors.

- **Circle Representation**: The size of each circle represents the level of motif enrichment (-log10 P-value), and the color indicates the RNA expression level of the corresponding transcription factor.

- **Key Findings**: 
    - Motifs of transcription factors like KLF family and FOXP1 were highly enriched in specific cell populations (e.g., Undiff.SPG-1).
    - The enrichment suggests roles these transcription factors might play in spermatogenesis stages.

### Methodology for Motif Identification

The article outlines a comprehensive process for motif identification:

1. **Identification of Open Chromatin Regions (NDRs)**:
    - Using scCOOL-seq, open chromatin regions are identified in each cell type.

2. **Searching for Transcription Factor Motifs**:
    - HOMER software's `findMotifsGenome.pl` command is used to search for transcription factor binding motifs in these open regions.
    - Specific parameters: Searching for motifs with lengths of 8-12bp (`-size given -len 8,10,12`).

3. **Selection Based on Enrichment and Expression Levels**:
    - Motifs are reported based on their enrichment (P-value ≤ 10^-10) and expression levels (TPM ≥ 10 in at least one cell type).

### Conclusion

This methodology, combining single-cell chromatin accessibility data and motif analysis, helps predict key transcription regulatory factors by utilizing open chromatin regions as prior information and correlating motif enrichment with transcription factor expression levels. The study of motifs, especially in NDRs, provides significant insights into regulatory mechanisms controlling gene expression during different cell development stages like spermatogenesis.


## P5 and P7 adaptors

Regardless of the library construction method, submitted libraries will consist of a sequence of interest flanked on either side by adapter constructs. On each end, these adapter constructs have flow cell binding sites, P5 and P7, which allow the library fragment to attach to the flow cell surface. All Paired-End Format sequencing on the HiSeq and All sequencing of any type on the MiSeq MUST HAVE FULL-LENGTH P5 and P7 sequences . (some of the small RNA libraries and alternative genomic library constructions use a partial P7, this is not supported by the HiSeq PE and MiSeq.)

## Pore-C Technology Overview

### Introduction
Pore-C is an innovative genomic technique that enhances our understanding of the 3D genome structure. While sharing initial steps with Hi-C, Pore-C extends capabilities to capture complex multi-way DNA interactions, offering a more comprehensive view of chromatin architecture.

### Comparison with Hi-C
- **Shared Initial Steps**: Both Pore-C and Hi-C begin with similar wet lab procedures, including cell culture, cross-linking, DNA digestion, and ligation.
- **Divergence in Sequencing and Analysis**:
  - **Sequencing Method**: Pore-C employs long-read nanopore sequencing, contrasting Hi-C's short-read Illumina sequencing.
  - **Data Complexity**: Pore-C captures multi-way DNA contacts, providing richer spatial chromatin structure information, unlike Hi-C which focuses on pairwise interactions.

### Experimental Workflow
1. **Cell Culture and Cross-Linking**: Similar to Hi-C, cells are grown and treated with formaldehyde for DNA-protein cross-linking.
2. **DNA Digestion and Ligation**: DNA is digested with restriction enzymes and ligated, mirroring Hi-C's approach.
3. **Key Divergences Post-Ligation**:
   - **DNA Purification and Sequencing**: Post-ligation steps involve DNA purification followed by long-read nanopore sequencing.
   - **Bioinformatics Analysis**: Pore-C requires specialized analysis pipelines to handle long-read data and multi-way contact mapping.

### Bioinformatics Workflow (Pore-C Specific)
1. **Local Alignments to Reference(s)**: DNA reads are aligned to the reference genome using BWA.
2. **Optimize Alignment Path Through Read**: Identifying the most probable DNA segment combinations covering each read.
3. **Assign Alignments to Restriction Fragments**: Alignments are mapped to specific restriction fragments.
4. **Assign Restriction Fragments to Bins**: Fragments are categorized into larger genomic bins.
5. **Tabulate Support for Bin-Bin Contacts**: Calculating interaction frequencies between genomic bins.
6. **Generate Contact Map**: Creating a comprehensive chromatin contact map from the interaction data.

### Advantages of Pore-C
- **Detailed Chromatin Interactions**: Captures complex, high-order DNA interactions.
- **Long-Read Sequencing**: Provides more extensive coverage and insight into chromatin structure.
- **Enhanced Analytical Depth**: Offers deeper insights into genome organization and gene regulation.

### Conclusion
Pore-C represents a significant advancement in chromatin conformation studies. By building on the foundational steps of Hi-C and incorporating advanced sequencing and analysis techniques, Pore-C provides a more nuanced and comprehensive understanding of the genome's 3D structure.

<div align=center>
<img src="imgs/porec.png">
</div>


## Post-Bisulfite Adapter Tagging (PBAT) - core technique of scBS-seq, scNanoCOOL-seq, scCOOL-seq

PBAT is a technique used in DNA methylation studies, specifically designed for bisulfite-treated DNA (BS-DNA), which allows for the analysis of methylation patterns at a single-nucleotide resolution. The figure provided illustrates the PBAT method, which is adapted for sequencing on Illumina platforms. 

 <div align=center>
<img src="imgs/PBAT1.png">
</div>

 <div align=center>
<img src="imgs/PBAT2.png">
</div>

### Key Steps in the PBAT Process:

1. **Random Priming**: Initially, a random primer with a sequence of NNNN (where N represents any nucleotide) is annealed to the BS-DNA at random locations across the genome.

2. **First Strand Synthesis**: The 5' end of the random primer is extended using a DNA polymerase, creating the first DNA strand complementary to the bisulfite-converted DNA template.

3. **Second Random Priming**: After the synthesis of the first strand, a second random priming event occurs, where a new random primer with an NNNN sequence anneals to the newly synthesized first strand DNA.

4. **Second Strand Synthesis**: This step extends the second primer, creating a double-stranded DNA with adapters at both ends.

### Adapter Design for Illumina Sequencing:

In the context of Illumina sequencing, the oligo sequences are designed to include the P5 and P7 sequences which are standard Illumina adapters. These are crucial for the DNA library to properly bind to the flow cell and for bridge amplification to occur.

- The **P5+N6** sequence corresponds to the first adapter, which is bound to the first strand synthesized from the BS-DNA.
- The **P7+N6** sequence corresponds to the second adapter, which is attached to the complementary strand generated from the first cDNA.

### Importance for Illumina Sequencing:

The inclusion of P5 and P7 sequences in the PBAT process is essential to align the DNA fragments correctly onto the Illumina flow cell. These sequences are recognized by the Illumina sequencer, allowing the bound DNA fragments to undergo bridge amplification, which is a hallmark of the Illumina sequencing process.

By integrating PBAT with Illumina’s sequencing technology, researchers can efficiently analyze genome-wide methylation patterns at single-base resolution in single cells or low-input DNA samples, which is critical for understanding epigenetic regulation and its impact on gene function and disease.

## Percent spliced in (PSI) for short reads

The percent spliced in index (PSI) indicates the efficiency of splicing a specific exon into the transcript population of a gene. It shows the frequency/degree of a specific alternative splicing event. It is calculated by IR/(IR + ER)%, where IR is inclusion reads, ER is exclusion reads. IR means reads that support the event. ER means reads that do not support the event.

However, without normalization, the PSI of long exons would always approximate 100%. We need to normalize them by dividing the possible mapping length of a certian reads. For short reads NGS, the read length is fixed, like 100 bp, 150 bp.

So normalized IR = IR/(exon length + (read length -1)), normalized ER = ER/(read length -1).

How to understand this? For a certain read, the possible mapping read range is the exon length + read length - 1. The read can be mapped to either the upstream or downstream of the exon. For ER, the maximum range is the read length - 1. It has nothing to do with the exon length.

<b>The traditional PSI can not be calculated for long read sequencing</b>, since the read length is not fixed. You can calulate just the number of IR and ER without normalization. This may not be that strict, but it can work.

<div align=center>
<img src="imgs/PSI1.png">
</div>

<div align=center>
<img src="imgs/PSI2.png">
</div>

## Precision and Recall (Sensitivity)

Precision is calculated by dividing the true positives by anything that was predicted as a positive. 

Recall (or True Positive Rate) is calculated by dividing the true positives by anything that should have been predicted as positive.

Recall and Sensitivity are one and the same.

## Primary template-directed amplification (PTA)

An isothermal WGA method that reproducibly captures >95% of the genomes of single cells in a more uniform and accurate manner than existing approaches, resulting in significantly improved variant calling sensitivity and precision.

1. PTA reads are short (they are generated intentionally to avoid uneven amplification) and PTA is not compatible with long-read sequencing.  
2. PTA does create fewer artifacts (including chimera) than MDA, and this is experimentally verified in papers (https://www.nature.com/articles/s41588-022-01180-2). However, they focus on artifactual SNVs and small INDELs induced by ingle-strand dropout (SSD). 
3. Accompany with the PTA, the authors developed a tool SCAN-SNV (https://www.nature.com/articles/s41467-019-11857-8) to detect somatic SNV and small INDELs aware of MDA artifact and a tool SCAN2 (https://www.biorxiv.org/content/10.1101/2021.04.30.442032v1.abstract) to detect somatic SNV and small INDELs aware of PTA artifact. 
4. The PTA is provided as the ResolveDNA Kit in BioSkryb company, with price unlisted. They are provided as 24 and 96 Reaction Kits.  

PTA is provided as  ResolveDNA® Whole Genome Sequencing Workflow, with low-input (>4 pg to 10 ng) DNA samples.For the whole pack you need ResolveDNA Whole Genome Amplification Kit - 96 Reactions, ResolveDNA Bead Purification Kit - 96 Reactions, ResolveDNA Library Preparation Kit - 96 Reactions, ResolveDNA Multi-Use Library Adapters - 960 Reactions, ResolveDNA Cell Buffer Pack, 25 Low Bind 96 Well PCR Plates, PCR Plate Sealing Film - Pack of 100

## P2
Component	Specicification
Processor	Latest Generation CPU processor
GPU	Latest generation Ampere NVIDIA GPU for basecalling acceleration
Storage	16 TB SSD
RAM	128 GB

Self-contained benchtop device with compute inside
Powerful GPU and built-in screen
Can run up to two PromethION flow cells
Simple, preconfigured device - plug in and start running

## P 2 Solo
Up to 580 Gb (both flow cells sequencing). 290 Gb per flow cell. PromethION 2 Solo is a modular nanopore sequencing device using the same technology found in the MinION and GridION devices. It allows up to two sequencing experiments to be run concurrently or individually. The device has no integrated compute, but can be plugged into a GridION Mk1 or a stand-alone computer that meets the minimum spec.

## PromethION P24 and P48 
Up to 7,000 Gb (P24) or 14,000 Gb (P48) (all 24 or 48 flow cells sequencing respectively). 290 Gb per flow cell.

## Promoter

A promoter, as related to genomics, is a region of DNA upstream of a gene where relevant proteins (such as RNA polymerase and transcription factors) bind to initiate transcription of that gene. The resulting transcription produces an RNA molecule (such as mRNA).

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/119dccf8f9f3b544668ba5945fb3206ac92180a5/imgs/promoter2.png">
</div>

## Prostate Cancer Cell Lines & AR Status

**AR-Positive Cell Lines:**
- **LNCaP:** Sensitive to androgens. Good for anti-androgen therapy research.
- **VCaP:** Has ERG gene changes. Helps study AR's role.
- **22Rv1:** Grows despite low androgens. Important for understanding treatment resistance.

**AR-Negative Cell Lines:**
- **DU145:** Does not respond to androgens. Used for other cancer growth mechanisms.
- **PC3:** Similar to DU145, important for studying cancer spread and invasion.

**Special Mention:**
- **RWPE-1:** Normal prostate cells for control studies.

### Importance
- **AR-Positive:** Crucial for androgen-related research and therapy development.
- **AR-Negative:** Offers insights into androgen-independent cancer progression.
- **RWPE-1:** Provides a baseline for comparison.

## REPLI-g Single Cell Kit

Yields up to 40 µg/reaction, average product length >10 kb.

## REPLI-g Advanced DNA Single Cell Kit

Yield from a single cell is 25–35 µg of amplified DNA. Saves at least 1 hour of time versus the first-generation REPLI-g Single Cell Kit (cat. nos. 150343 and 150345).

## REPLI-g UltraFast Mini Kit

Resulting in typical DNA yields of 7 μg per 20 μl reaction. Sufficient product is available for downstream genetic analysis after just 45 minutes. Input, 10 ng DNA, 0.5 µl whole blood, ~300 cells/µl.

Certainly, I can add that for you:

## Restriction Enzymes in Molecular Cloning

Restriction enzymes, also known as restriction endonucleases, are key tools in molecular biology that facilitate the process of cloning by cutting DNA at specific sequences. Here's a concise overview of their role and considerations in molecular cloning.

### Characteristics of Restriction Enzymes
- **Specificity**: Each enzyme recognizes specific nucleotide sequences, typically 4 to 8 bases in length, often palindromic, and cuts the DNA at or near these sites.

### Cleavage Types:
- **Blunt Ends**: Enzymes make a straight cut across the DNA, leaving no overhangs.
- **Sticky Ends**: Enzymes make staggered cuts resulting in overhangs that can be either 5' or 3', facilitating the ligation of complementary sequences.

### Usage Frequency
- **Sticky End Enzymes**: Preferred for their precise and efficient ligation capabilities, hence widely used in cloning and recombinant DNA techniques. They are also the preferred choice for 3C (Chromosome Conformation Capture) techniques due to their efficient ligation of DNA ends that are in close proximity in the three-dimensional nucleus.
- **Blunt End Enzymes**: Useful in certain cloning strategies, such as TA cloning or when using ligases like T4 DNA Ligase that can join any type of ends.

### Selection in Cloning
- **Directional Cloning**: Two enzymes with different sticky ends are typically chosen to ensure correct orientation of the insert.
- **Optimal Conditions**: Each enzyme requires specific conditions like pH, ionic strength, and temperature for optimal activity.

## Reverse ChIP-seq Workflow

Reverse ChIP-seq is a technique to identify the DNA-binding proteins associated with specific genomic regions. The workflow includes:

1. **Crosslinking**: Fixing DNA-protein interactions with formaldehyde.
2. **Chromatin Shearing**: Fragmenting chromatin using sonication.
3. **Probe Hybridization**: Binding biotinylated DNA probes to target sequences.
4. **Capture**: Isolating probe-bound chromatin with streptavidin-coated beads.
5. **Identification**: Analyzing captured proteins via mass spectrometry.

<div align=center>
<img src="/imgs/reverse-chip.png">
</div>

The distribution of probes is usually uniform and covers the entire region, for example, the researchers designed 12 biotin-labeled probes to capture the 1281 bp AtCAT3 promoter region, with 6 probes for each DNA strand evenly distributed.

## RNA Immunoprecipitation-Sequencing (RIP-Seq) Protocol

### Overview
RNA Immunoprecipitation-Sequencing (RIP-Seq) is an essential technique for studying RNA-protein interactions and understanding gene regulation. It involves isolating RNA-protein complexes and analyzing associated RNA molecules. This note summarizes the protocol, focusing on whether to use cross-linking and the differences between cross-linking methods.

### Cross-Linking in RIP-Seq
- **Optional Use**: Cross-linking in RIP-Seq is optional. The absence of cross-linking simplifies the procedure but may miss weak RNA-protein interactions.
- **Formaldehyde and UV Cross-Linking**:
  - **Formaldehyde**: Creates reversible covalent bonds between proteins and RNA. Requires reverse cross-linking to break these bonds.
  - **UV**: Forms non-reversible covalent bonds, thus, no reverse cross-linking is needed.
- **Proteinase K Treatment**: Regardless of the cross-linking method used, treatment with Proteinase K is essential for digesting and removing proteins before RNA extraction.

### Procedure
1. **Cell Preparation**: Grow cells to 80-90% confluence.
2. **Cross-Linking (Optional)**:
   - **Formaldehyde Method**: Treat cells with formaldehyde, followed by quenching with glycine.
   - **UV Method**: Expose cells to UV light to fix RNA-protein interactions.
3. **Cell Lysis**: Use appropriate buffers for cell and nuclear lysis.
4. **Immunoprecipitation**: Incubate cell lysates with antibody-coated magnetic beads to pull down RNA-protein complexes.
5. **Washing**: Remove non-specifically bound material.
6. **Proteinase K Treatment**: Digest proteins to release RNA.
7. **RNA Extraction**: Use phenol:chloroform extraction and ethanol precipitation to isolate RNA.
8. **RNA Analysis**: Perform reverse transcription followed by PCR, microarray, or sequencing.

### Notes
- **Quality of Antibody**: The success of RIP-Seq greatly depends on the quality of the antibody used.
- **Native vs. Cross-Linked RIP**: Choose based on the nature of RNA-protein interactions and experimental requirements.
- **Bioinformatic Analysis**: Essential for mapping reads to transcripts and identifying binding sites.

<div align=center>
<img src="/imgs/ripsew.png">
</div>

For the first time doing RIP-Seq, in addition to the experimental group, a positive control can also be done, which is an antibody interacting with known RNA. The input should definitely be saved, and IgG can also be saved, but IgG is mainly used for qPCR validation (PCR various RNAs in the IgG pull-down samples), and the positive control is mainly used for qPCR validation as well. For qPCR validation, target primers need to be designed, such as the target RNA of your antibody-binding protein, the known interacting RNA with your antibody-binding protein, and the RNA interacting with the antibody of the positive control.

When sending for sequencing, typically only the experimental group and input are sequenced.

Additionally, peak calling and differential peak analysis in RIP-Seq are two different things. Peak calling is comparing the experimental group and input samples, identifying specific protein-DNA interaction regions within each group, which is a question of presence or absence.

On the other hand, differential peak analysis is a biological question, typically involving three experimental groups vs. three treatment groups, each compared with their own input after which peaks are extracted. Then, peaks from the experimental group and the treatment group are compared to see if there is differential expression.

In the qPCR graph, the y-axis is the ratio of RIP/Input, usually around 2%. Then, generally speaking, if the experiment is higher than 2%, and higher than the qPCR positive control (the same antibody), and five times higher than IgG, it is considered to have interaction. Here we are looking at presence or absence. It is not necessarily much higher than input, depending on the antibody, and it is even possible to be lower than input.


## How to Calculate and Plot RIP-qPCR Bar Charts

There are two main methods to analyze RIP-qPCR data:

### Percentage of Input
This method quantifies the enrichment of RNA as a percentage of the total input. To calculate the percentage of input, you need the Ct values from the 10% input, the antibody of interest, and the IgG control.

### Fold Enrichment
Alternatively, you can calculate the fold enrichment of your protein of interest over the IgG control. This method only requires the Ct values from the antibody of interest and the IgG control, without the need for the 10% input.

For the detailed calculation methods, please refer to: [Top Tip Bio's guide on analyzing ChIP-qPCR data](https://toptipbio.com/analyse-chip-qpcr-data/).


## RIP-qPCR vs RT-qPCR Analysis

### RIP-qPCR Analysis

RIP-qPCR is a technique used to determine whether a specific RNA is associated with a particular protein, indicating an RNA-protein interaction. The analysis typically follows these steps:

- **RNA Extraction**: Extract RNA from the cell lysate before immunoprecipitation (IP). A fraction of this lysate is taken as the 'input' sample.
  
- **cDNA Preparation**: Synthesize cDNA from RNA extracted from both IP and input samples.

- **qPCR**: Perform qPCR using specific primers to quantify the levels of RNA of interest in both the IP and input samples.

- **Data Calculation**: Calculate the enrichment of the RNA-protein interaction using either the 'Percentage of Input' or 'Fold Enrichment' method based on the Ct values from the IP, IgG, and input samples.

The purpose of RIP-qPCR is to assess the presence and level of interaction between the target RNA and protein, rather than quantifying the absolute expression levels of RNA.

### RT-qPCR Analysis

RT-qPCR quantifies the expression levels of specific genes from total RNA. The typical steps involved are:

- **RNA Extraction**: Extract total RNA using reagents like TRIzol.

- **cDNA Synthesis**: Synthesize cDNA from the total RNA.

- **qPCR**: Quantify the levels of gene expression using qPCR with specific primers.

- **Data Calculation**: The Ct value for each gene is normalized to a reference or housekeeping gene (e.g., GAPDH) to account for sample-to-sample variation.

For relative quantification in RT-qPCR, calculating amplification efficiency and generating a standard curve are not always necessary, particularly if the comparison is made using the ΔΔCt method with a suitable housekeeping gene as the reference.

### Summary

- **RIP-qPCR** focuses on the relative enrichment of target RNA associated with a specific protein.
  
- **RT-qPCR** measures the relative expression levels of genes using a reference gene for normalization, with a simple calculation method like the ΔΔCt being sufficient for relative quantification.

Normalization to a housekeeping gene like GAPDH is common in RT-qPCR for gene expression studies, while in RIP-qPCR, the enrichment calculation is based on the comparison of target RNA in IP versus the input sample.

## RIP-seq Peak Calling Methods

In the analysis of RNA Immunoprecipitation sequencing (RIP-seq) data, two primary methods of peak calling are employed to identify RNA fragments that are bound to specific proteins. These methods are crucial for understanding the interactions between RNA and proteins and their roles in gene regulation.

### Using Peak Calling Software

This approach involves using software designed to identify enriched RNA fragments in RIP-seq experiments. These tools analyze IP (Immunoprecipitation) samples and control samples (such as input or negative control) to identify regions of RNA that are specifically bound to proteins, known as "peaks". Commonly used software packages include MACS (Model-based Analysis of ChIP-Seq), HOMER, and CLIPper. Originally designed for ChIP-seq or CLIP-seq data, these tools can also be adapted for RIP-seq data by adjusting parameters.

### Treating RIP-seq Data as Differential Expression Analysis

This method is analogous to processing RNA-seq data. Here, researchers compare IP samples with control samples (such as input) and use methods of differential expression analysis to identify significantly enriched or depleted RNA fragments. This approach may involve standard RNA-seq data analysis workflows, using software like DESeq2 or edgeR, to identify regions with significant differential expression. A key aspect of this method is the processing and normalization of data to ensure valid comparisons.

DESeq2 [[29]] is used to identify differentially bound genes (DBGs). The read count for each gene across all samples is input to compare IP with input samples. A false discovery rate (FDR) of < 0.05 and a fold change of ≥ 2 (enriched) or ≤ 0.5 (depleted) are set as the cutoff criteria for identifying DBGs.

Programs for de novo transcript assembly followed by differential expression (DE) analysis, such as the Cufflinks/Cuffdiff suite [[15, 16]], and for DE on a set of known transcripts, such as DESeq [[17]], may seem applicable to RIP-seq analysis. However, unlike peak-calling strategies, these transcript-based methods assume the entire transcriptome is being sequenced at fairly deep coverage (as is usually the case in RNA-seq) and thus may be sensitive to background noise typical of IP-based protocols.

## Key Points in RIP-seq Experimental Design

### Q: What is the source of samples in a RIP-seq experiment?
A: In the ideal RIP-seq experimental setup, each type of sample (RIP, IgG control, or Input) should originate from the same source to ensure they are directly comparable. This means that all samples for a given replicate come from the same pool of cells.

### Q: How is Input extracted for RIP-qPCR?
A: The Input for RIP-qPCR is not calculated as a percentage of the total volume but rather in absolute terms. If 10 μL is designated as the Input, the corresponding experimental or IgG sample would be 100 μL.

### Q: Can RIP-qPCR be performed on each replicate, and is the use of IgG necessary?
A: Yes, RIP-qPCR can be performed on each replicate. The use of IgG is optional; it can be included for quality control purposes or omitted entirely. If included, IgG can serve as a control to demonstrate the specificity of the RNA-protein interaction.

### Q: Is it mandatory to measure IgG and Input in sequencing lab groups?
A: While measurement of the Input is mandatory for sequencing experiments, the measurement of IgG is optional. If IgG is measured, it is typically used only for illustrative purposes and not for differential expression analysis.

### Additional Notes:
- Each replicate can undergo a single round of RIP-qPCR, where the use of IgG is optional.
- IgG, if measured, is typically used only for graphical representation and is not included in the RIP-seq peak calling or differential expression analysis.
- For peak calling, it is common practice to compare the IP samples directly with their corresponding Input samples to identify areas of enrichment.
- The same quantity of RNA from input, RIP, IgG were used for reverse transcription in RIP-qPCR.

<div align=center>
<img src="/imgs/rip-seq.png">
</div>

## Understanding Long-Read RIP-Seq Analysis

### Long-Read vs Short-Read Sequencing in RIP-Seq

RIP-seq (RNA Immunoprecipitation Sequencing) experiments typically involve the pulldown of RNA-binding proteins (RBPs) along with their associated RNAs, with the subsequent identification of these RNAs via sequencing. While short-read sequencing has been a common approach, long-read sequencing provides a distinct advantage for RIP-seq analysis.

- **Long-Read Sequencing**: Offers continuous reads that cover entire RNA transcripts, capturing the full landscape of RNA-protein interactions. It is not dependent on the identification of peaky binding sites, as it reveals the entire RNA molecule.
- **Short-Read Sequencing**: Generates shorter fragments of RNA, traditionally requiring peak calling to identify regions of protein-RNA interaction.

### Analysis Differences: RNA-Seq vs Peak Calling

- **RNA-Seq Analysis**: When applying RNA-seq analysis tools such as EdgeR or DESeq2 to long-read RIP-seq data, we shift our focus from peak calling to differential expression analysis. This allows us to compare the abundance of RNAs between input (baseline RNA) and IP (pulldown) samples, aiming to identify RNAs that are preferentially bound by the protein of interest.
- **Peak Calling**: Traditional peak calling, used in short-read sequencing and techniques like ChIP-seq, is not suited for long-read RIP-seq data due to the continuous nature of the reads. Tools designed for punctate or peaky signals (identifying precise binding sites) are not appropriate for full-transcript capture.

### Practical Considerations for Long-Read RIP-Seq Analysis

- **Experimental Design**: It's crucial to understand the nature of the experiment. Whether the RNA is fragmented before pulldown (similar to CLIP techniques) or captured in full length (true RIP-seq) will determine the analysis approach.
- **Background Noise Sensitivity**: Transcript-based methods such as DESeq and EdgeR may be sensitive to the background noise typical in IP-based protocols. They assume deep coverage sequencing of the full transcriptome, which might not be the case in all RIP-seq protocols.
- **RIPSeeker**: Among the discussed methods, only RIPSeeker is set up to handle full-transcript RIP-seq data. It's advisable to consider the suitability of this tool for your specific data type.

### Analysis Strategy for Long-Read RIP-Seq

When dealing with long-read RIP-seq data, the goal is to identify RNAs that show significant enrichment in the IP samples compared to the input. This is indicative of RNA-protein interactions rather than the identification of localized peaks.

- **Differential Expression Analysis**: DESeq2 can be utilized for comparing groups of input and IP samples (e.g., 3 vs 3 comparisons) to identify RNAs that are significantly enriched in the IP samples.
- **Terminology Clarification**: While the term "peak calling" is commonly used, in the context of long-read RIP-seq, we are not identifying traditional peaks but rather enriched RNA transcripts indicative of binding events.

Below is RIP-seq for long read
<div align=center>
<img src="/imgs/rip-long.png">
</div>


## RNA-Binding Proteins (RBPs)

RNA-binding proteins (RBPs) play a crucial role in various aspects of RNA metabolism and gene regulation. They interact with different types of RNA molecules, influencing their function, localization, stability, and processing.

## RNA quantity conversions from tissue and cells

<div align=center>
<img src="/imgs/rnaq.png">
</div>

A 80% confluency VCaP 10 cm dish can get 30-40 ug total RNA.


### Functions of RBPs

#### RNA Processing and Splicing
- RBPs are involved in the processing of pre-mRNA into mature mRNA.
- They participate in RNA splicing, editing, and modification processes.

#### RNA Transport and Localization
- RBPs facilitate the transport of RNA molecules to specific locations within the cell.
- This localization is crucial for localized protein synthesis and cellular function.

#### RNA Stability and Degradation
- RBPs can increase or decrease the stability of RNA molecules.
- They play a role in determining the half-life and degradation of RNAs.

#### Protein Synthesis Regulation
- By binding to mRNA, RBPs regulate protein synthesis.
- They modulate gene expression through promotion or inhibition of translation.

#### Involvement in Non-Coding RNA Functions
- RBPs interact with non-coding RNAs like miRNAs and lncRNAs.
- They are involved in gene silencing, chromatin remodeling, and other regulatory processes.

### Interaction with Various RNA Types

#### Messenger RNA (mRNA)
- RBPs regulate the translation, stability, and processing of mRNA.
- They are key players in post-transcriptional gene regulation.

#### Ribosomal RNA (rRNA) and Transfer RNA (tRNA)
- RBPs are involved in the structural formation and function of rRNA and tRNA in protein synthesis.

#### Small Nuclear RNA (snRNA) and Small Nucleolar RNA (snoRNA)
- These RNAs, in association with RBPs, are involved in RNA splicing and modification.

#### Long Non-Coding RNA (lncRNA)
- RBPs interact with lncRNAs in gene expression regulation, chromatin structure organization, and other cellular processes.

#### MicroRNA (miRNA)
- RBPs play a role in the maturation and function of miRNAs, impacting gene silencing and post-transcriptional regulation.

### Importance of RBPs

RBPs are integral to cellular functioning and gene expression regulation. Their study helps in understanding RNA roles and regulatory mechanisms in biology and disease. Alterations in RBP interactions with RNA can lead to various diseases, making them a significant focus in biomedical research.

## RT-qPCR Analysis Example for Tissue Samples

When conducting RT-qPCR analysis on tissue samples for relative expression comparison between treated and control groups, a stable reference gene, such as GAPDH or ACTB, is typically sufficient. The reference gene serves to normalize the data, compensating for technical variations and differences in RNA amounts across samples.

Here is a specific example:

Suppose we have tissue samples from both treated and control groups and wish to measure the relative expression levels of a particular gene, Gene A. Each group has three replicates.

**Treated Group (Treated) Ct Values**:

- Gene A: 24.0, 23.8, 24.2
- Reference GAPDH: 18.0, 18.1, 18.1

**Control Group (Control) Ct Values**:

- Gene A: 26.0, 26.2, 26.1
- Reference GAPDH: 18.0, 18.2, 18.0

### Calculation Steps:

1. **Calculate the ΔCt for each sample** (Gene A Ct value - GAPDH Ct value):

   **Treated Group**:
   
   - ΔCt = 24.0 - 18.0 = 6.0
   - ΔCt = 23.8 - 18.1 = 5.7
   - ΔCt = 24.2 - 18.1 = 6.1

   **Control Group**:
   
   - ΔCt = 26.0 - 18.0 = 8.0
   - ΔCt = 26.2 - 18.2 = 8.0
   - ΔCt = 26.1 - 18.0 = 8.1

2. **Calculate the average ΔCt value** for each group:

   **Treated Group** average ΔCt: (6.0 + 5.7 + 6.1) / 3 = 5.93
   
   **Control Group** average ΔCt: (8.0 + 8.0 + 8.1) / 3 = 8.03

3. **Calculate ΔΔCt** (average ΔCt of the Control group - average ΔCt of the Treated group):

   ΔΔCt = 8.03 - 5.93 = 2.1

4. **Calculate the relative expression level** (using the negative exponent of 2^(-ΔΔCt)):

   Relative expression level = 2^(-ΔΔCt) = 2^(-2.1) ≈ 0.25

This means that Gene A's expression level in the treated group is approximately 25% of that in the control group, indicating a significant decrease in Gene A expression due to the treatment.

The calculations assume that the reference gene GAPDH is consistently expressed across all samples, which should be validated in the experiment. This example demonstrates how to use RT-qPCR data and a reference gene to assess relative changes in gene expression.

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
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a39ab01d88f7f5f998c7bd3558a055148faeb8d8/imgs/transindel_Algorithm2.png">
</div

.  

## Single-cell Bisulfite Sequencing (scBS-seq)

<div align=center>
<img src="/imgs/scBS-seq.png">
</div

<div align=center>
<img src="/imgs/scbs.png">
</div

### **Key Insight**

**During the initial oligo1 random priming in scBS-seq, various product combinations are generated. Oligo1 acts not only on the template but also on the newly synthesized products, leading to a diverse range of amplification products.**

### Introduction
Single-cell bisulfite sequencing (scBS-seq) is a specialized protocol for analyzing DNA methylation at single-cell resolution. It is a modification of traditional bisulfite sequencing and incorporates unique steps suitable for single cells. The process involves treating isolated single-cell genomic DNA with sodium bisulfite, followed by random priming and PCR amplification for sequencing.

### Detailed Workflow

#### Cell Lysis and Bisulfite Conversion

- **Cell Lysis**: Cells are lysed using a protein lysis buffer, and genomic DNA is released.
- **Bisulfite Conversion**: The Imprint DNA Modification Kit (Sigma, Cat# D5044) is used for bisulfite conversion, modifying unmethylated cytosines to uracil while leaving methylated cytosines unchanged. This process involves denaturation and incubation steps at specific temperatures.

#### DNA Priming and Amplification

- **Random Priming**: Post-conversion, DNA undergoes multiple rounds of random priming. Key primers include:

  - **Oligo1 Primer**: `[Btn]CTACACGACGCTCTTCCGATCTNNNNNNNNN`
  - **Oligo2 Primer**: `TGCTGAACCGCTCTTCCGATCTNNNNNNNNN`

- These primers facilitate the synthesis of new DNA strands from the bisulfite-treated DNA.

#### PCR and Library Preparation

- **PCR Amplification**: Using primers such as the PE1.0 forward primer and indexed iPCRTag reverse primer, the DNA is amplified.
  - **PE1.0 Forward Primer**: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- Amplified libraries are purified and quantified, ensuring only high-quality DNA is used for sequencing.

#### Sequencing and Analysis

- The prepared libraries are then sequenced using deep sequencing technology, providing a high-resolution profile of methylated cytosines at a single-cell level.

### Applications

scBS-seq is invaluable in epigenetics for understanding methylation patterns in individual cells, crucial for studying:

- Cellular heterogeneity and epigenetic variability.
- Developmental processes and disease progression at the single-cell level.

By providing detailed methylation profiles of individual cells, scBS-seq offers insights into the complex regulation of gene expression and contributes significantly to the understanding of epigenetic mechanisms in health and disease.

## Single-Strand cDNA

Single-strand cDNA (ss-cDNA) is synthesized from an RNA template by the process of reverse transcription. The enzyme reverse transcriptase catalyzes the formation of ss-cDNA by using RNA as a template and incorporating deoxyribonucleotides complementary to the RNA template. This form of cDNA is primarily used in quantitative polymerase chain reaction (qPCR) assays, a technique to measure the abundance of a particular mRNA or to compare mRNA levels in different sample sets in a quantitative and sensitive manner. The advantage of using ss-cDNA for qPCR is that it represents a precise copy of the mRNA, allowing for accurate amplification and quantification of the target sequence.

## siRNA

siRNA, 21-23 nucleotides in length and double-stranded, enters cells, may be processed by Dicer, then forms a RISC complex; it becomes single-stranded within RISC, then targets and degrades complementary mRNA.

## scCOOL-seq (Chromatin Overall Omic-scale Landscape Sequencing) 

<div align=center>
<img src="/imgs/sccool-seq.png">
</div

DNA sequences are made up of four bases: A, T, C, and G. CG sites refer to the occurrence of CG base pairs in DNA sequences. For example, in the sequence AGCTCGAT, there is one CG site.

In mammalian genomes, DNA methylation primarily occurs at CG sites (CpG dinucleotides), meaning methylation mainly occurs on these CG base pairs. This is a key characteristic of DNA methylation in mammals.

In the study, the enzyme M.CviPI is used, which can introduce artificial methylation at GCH sites (where H is A, C, or T). This means it adds CH3 groups to DNA sequence sites like GCA, GCC, GCT, etc., but does not affect other sequences, such as CpG (WCG) sites.

In single-cell COOL-seq technology, the specificity of M.CviPI enzyme is utilized to distinguish between open and closed chromatin regions. In open chromatin regions, GCH sequences are accessible to M.CviPI enzyme, so these sites are methylated. In closed chromatin regions, where DNA is tightly packed, M.CviPI enzyme cannot access these sites, so they remain unmethylated.

WCG sites include all CpG sites (as 'W' can be A or T), covering potential methylation sites more comprehensively. Since M.CviPI enzyme does not act on CpG sites, its use does not affect the state of endogenous genome methylation, allowing researchers to independently assess endogenous DNA methylation (mainly occurring at CpG sites) and artificially introduced methylation (occurring at GCH sites) in the same cell.

Thus, by detecting Cs that are not converted after bisulfite treatment (indicating methylation), one can distinguish artificially methylated GCH sites and endogenously methylated WCG sites:

GCH sites detected with methylation - This is introduced by M.CviPI enzyme, indicating the region has open chromatin.
WCG sites detected with methylation - This represents the natural state of endogenous methylation.
This approach enables simultaneous detection of chromatin openness and endogenous methylation levels.

### **The scCOOL-seq Protocol Key Steps**

- **In vitro Methylation**  
  This step methylates all the cytosines in the genome, allowing differentiation between endogenous methylation (which will be protected) and open chromatin (which will be methylated).

- **Bisulfite Conversion**  
  This converts unmethylated cytosines to uracils, while leaving methylated cytosines unchanged.

- **Amplification with P5-N6-oligo1**  
  This adds the P5 adapter sequence needed for Illumina sequencing to the molecules. The N6 random sequence helps maintain library complexity. This initial amplification is done with linear amplification to maintain the integrity of the molecules.

- **Second Strand Synthesis with P7-N6-oligo2**  
  This converts the molecules to double-stranded DNA and adds the P7 adapter sequence.

- **PCR Amplification**  
  This does a standard PCR amplification of the library using primers complementary to the P5 and P7 adapter sequences to generate sufficient material for sequencing.

### **Why Only One Primer in the First Amplification Round**

- The first round of amplification using only one primer, such as scBS-seq-P5-N6-oligo1, is feasible in some special PCR applications. This strategy is commonly used for specific types of sequence amplification, known as Single Primer Amplification (SPA). In this case, the PCR process is different from traditional double-primer PCR. Here is a brief explanation of its working principle:
  
  - **Template Specificity**: In single primer amplification, the primer is designed to specifically bind to one end of the target DNA. This usually involves specific sequence design, enabling the primer to effectively bind with DNA treated with bisulfite.
  
  - **Amplification Mechanism**: Single primer PCR only synthesizes one strand in the first cycle. In subsequent cycles, the newly synthesized DNA strand itself serves as the template, while the original template strand is no longer involved in amplification. This method produces single-stranded DNA (ssDNA), not double-stranded DNA (dsDNA) as in traditional PCR.

### **Amplification Target Areas in the First Step**

- The template for the first step of linear amplification in scCOOL-seq is DNA that has undergone bisulfite treatment.
  
- Bisulfite treatment converts unmethylated cytosines to uracils, while methylated cytosines remain unchanged.
  
- Therefore, on the template DNA, only the originally methylated-protected cytosines remain as normal C, while most of the other cytosines are converted to U.
  
- During linear amplification, P5-N6-oligo1 as a primer hybridizes with the template. It can only form stable base pairing with C bases on the template.
  
- It cannot form stable pairing with U. Therefore, P5-N6-oligo1 actually can only hybridize with regions on the template that were originally methylated.
  
- These methylation-protected areas include:
  
  - Endogenously methylated regions
  - Regions methylated by the M.CviPI methyltransferase enzyme that has no base preference (i.e., open chromatin regions)
  
- Therefore, the main areas amplified by linear amplification are these two areas, not the entire genome.

### **Differentiating Open Chromatin and Endogenous Methylation**

- **Open Chromatin Areas**:  
  These areas are targeted for methylation by M.CviPI enzyme prior to bisulfite treatment. M.CviPI is a methyltransferase enzyme that specifically recognizes GCH sequences (where H represents A, C, or T) and adds a methyl group to cytosine. Since bisulfite treatment does not change methylated cytosines, these methylated GCH sequences remain unchanged after treatment and can be detected by subsequent sequencing.

- **Endogenous Methylation Areas**:  
  These areas are naturally occurring methylation regions, typically in CpG islands (WCG sequences, where W represents A or T). In bisulfite treatment, unmethylated cytosines are converted to uracils (including areas methylated artificially in the previous step), while methylated cytosines remain unchanged. Therefore, after treatment, methylated CpG islands can be identified through the preserved WCG sequences.

## scNanoCOOL-seq

<div align=center>
<img src="/imgs/scnanocool-seq.png">
</div

### General Process

1. **In vitro GpC Methylation (IVM) of Individual Cells**
   - Each cell is placed in a PCR tube containing a cell lysis and methylase reaction mixture, including RNase inhibitor, GC reaction buffer, S-adenosylmethionine, GpC Methyltransferase (M.CviPI), and IGEPAL CA-630.
   - The lysate is incubated to allow in vitro methylation of the cell nuclei, with subsequent heat inactivation of M.CviPI.
   - This step is critical for differentiating endogenously methylated DNA and artificially methylated open chromatin areas.

2. **Physical Separation of Single-Cell DNA and RNA**
   - A DNA/RNA separation mixture is added to the lysate for efficient separation of nucleic acids.
   - The mixture contains RNase inhibitor, Triton X-100, Tween 20, DTT, superscript II first-strand buffer, and Dynabeads Myone Carboxylic Acid.
   - This results in the separation of RNA (in the supernatant) and relatively intact nuclei (on the beads).

### DNA Part

1. **Cell Lysis and Genomic DNA Release**
   - Protein lysis buffer is used to release genomic DNA from the nuclei contained within the beads.
   - Buffer components: M-digestion buffer, Lambda DNA, and protease K.

2. **Bisulfite Conversion**
   - Genomic DNA undergoes bisulfite treatment to convert unmethylated cytosines to uracils, preserving methylated cytosines.
   - This step is crucial for distinguishing methylated and unmethylated regions in the genome.

3. **Random Priming with Oligo 1-N6**
   - Oligo 1-N6 random primers are used to tag bisulfite-converted DNA, enhancing library complexity.
   - Four rounds of random priming on bisulfite-converted DNA to obtain oligo 1-tagged DNAs and skipped the oligo 2 priming,
   - It is random priming not PCR. The purpose of this method is to capture and amplify sulfite-treated DNA fragments, especially those regions containing methylated cytosines.

4. **Library Amplification and Purification**
   - PCR is used to amplify the tagged DNA, introducing a 24-nucleotide barcode for sample tracking.
   - PCR conditions and AMPure XP beads purification steps are optimized to enrich longer DNA fragments.

5. **Sequencing Preparation**
   - The library is prepared for Nanopore sequencing using a Ligation Sequencing Kit.
   - End-repair, dA-tailing, and adapter ligation steps are included before loading the DNA library onto the PromethION platform.

### RNA Part

1. **Cell Lysis and RNA Isolation**
   - Cells are lysed, and RNA is isolated from the cytoplasmic fraction.
   - Special buffers are used to ensure RNA integrity.

2. **Reverse Transcription**
   - RNA is reverse-transcribed into cDNA using oligo(dT) primers with unique barcodes for each cell.
   - This step includes the use of SuperScript II reverse transcriptase and RNase inhibitors.

3. **cDNA Synthesis and Amplification**
   - Second-strand cDNA synthesis is performed, followed by PCR amplification.
   - Biotinylated pre-indexed primers and IS primers are used for efficient and specific amplification.

4. **cDNA Library Purification and Sequencing**
   - Amplified cDNA is purified using AMPure XP beads and biotin-streptavidin capture system.
   - The cDNA library is then prepared for sequencing on the Illumina platform, using specific adapters and PCR conditions.

## scNanoCOOL-seq and Bisulfite-Treated Genome Sequencing

### Bisulfite Treatment Impact

- **Transformation by Bisulfite**: Bisulfite treatment converts unmethylated cytosines (C) to uracils (U), while methylated cytosines remain unchanged.
- **Selective Preservation**: Post-treatment, only the originally methylated-protected cytosines remain as normal C, most other cytosines are converted to U.

### Amplification Strategy

- **Random Priming with Oligo 1-N6**: scNanoCOOL-seq uses Oligo 1-N6 random primers (not P5-N6-oligo1) for random priming and amplification.
- **Selective Hybridization and Amplification**: This primer can initiate DNA synthesis across a broader region, but due to bisulfite treatment characteristics, it primarily hybridizes and amplifies with the remaining unconverted C (i.e., originally methylated regions).

### Limitation of Sequencing Scope

- **Partial Genomic Coverage**: Thus, scNanoCOOL-seq primarily measures and analyzes those regions in the genome that were methylated prior to bisulfite treatment, not the entire genome.
- **Focused Analysis**: This means the method mainly captures specific, methylation-status-related genomic information, rather than providing a comprehensive overview of the entire genome.

## scNanoCOOL-seq core idea

scnanocool-seq is an innovative sequencing method derived from bisulfite sequencing techniques, specifically adapted to single-cell epigenetic analysis. It leverages unique properties of DNA primers to selectively enrich for longer DNA fragments suitable for high-resolution methylation profiling.

## **Key Insight into Primer Design**

**In the initial Oligo 1-N6 random priming step of scnanocool-seq, the Oligo 1-N6 primer can anneal to both the template and extension products, generating a complex mixture of DNA fragments (different combination of oligo1 + template seq; oligo1 + oligo1). This mixture includes fragments where the target sequence is flanked by Oligo 1-N6 adapters on both ends. When these adapters are reverse-complementary, short sequences can self-anneal to form closed structures, preventing their amplification in subsequent PCR cycles. This selective process enriches the library with longer, linear fragments.**

### Primers:

- **Oligo 1-N6 Random Primer**: `5‘ - CTACACGACGCTCTTCCGATCTNNNNNN - 3’`
  - Used for initial random priming, this primer introduces adapter sequences at both ends of the DNA fragments.

- **24bp Barcode-Oligo1 Primer**: `5‘ - [24 bp barcode]CTACACGACGCTCTTCCGATCT - 3’`
  - Despite having a sequence identical to the Oligo 1-N6 adapter, this primer is used in PCR amplification to target longer fragments that have not self-annealed. The barcode allows for subsequent identification and multiplexing of samples.

### Amplification Process:

1. **Initial Priming Event**: The Oligo 1-N6 primer anneals to bisulfite-converted DNA at random sites, initiating the synthesis of various DNA products.
2. **Selective Enrichment**: Only those longer DNA fragments that are not self-annealed are efficiently amplified in PCR, due to the inability of circularized short fragments to act as templates.
3. **Barcode PCR Amplification**: PCR amplification using the 24bp Barcode-Oligo1 primer ensures that only the desired DNA fragments with the correct adapter configuration are amplified.

This strategic use of identical adapter sequences on both the Oligo 1-N6 primer and the Barcode-Oligo1 primer is critical for the scnanocool-seq method's ability to enrich for longer fragments and exclude short, self-annealed fragments from amplification, thereby optimizing the sequencing library for high-throughput analysis.

## scNMT-seq (single-cell nucleosome, methylation and transcription sequencing)

The main process of scNMT-seq can be summarized as follows:

1. Cell collection and methyltransferase reaction: Cells are collected into a methyltransferase reaction mixture and incubated at 37°C for 15 minutes. The reaction is stopped, and RNA is preserved using RLT plus.

2. scM&T-seq library preparation: mRNA is captured using Smart-seq2 oligo-dT, and the lysate containing gDNA is transferred to a separate PCR plate. Reverse transcription, PCR amplification, and cDNA purification are performed. Libraries are prepared using the Nextera XT Kit. In parallel, genomic DNA is purified and bisulfite-converted using the EZ Methylation Direct MagBead kit.

3. First and second strand synthesis of converted DNA, and library amplification: The converted DNA undergoes first strand synthesis using klenow exo- enzyme and is purified using AMPure XP beads. The second strand synthesis is performed, and products are purified again. The final PCR amplification is carried out, and scBS-seq libraries are purified using AMPure XP beads before pooling and sequencing.

The amplification will not remove the already-fixed methylation information.

<div align=center>
<img src="imgs/scNMT-seq.png">
</div

.
## Slippage during polymerase chain reaction amplification  

Slippage during PCR, also known as replication slippage, is a form of mutation that causes trinucleotide or dinucleotide amplification and sometimes even contraction during DNA replication, resulting in tandem repeat sequences.

## SMARCA1

SMARCA1, also known as SNF2L (Sucrose Non-Fermenting 2-Like), is a crucial gene in the human genome that encodes a protein playing a significant role in chromatin remodeling. It is a part of the SWI/SNF (SWItch/Sucrose Non-Fermentable) family of proteins, known for their involvement in altering chromatin structure to regulate gene expression.

### Function and Role

- **Chromatin Remodeling**: SMARCA1 encodes an ATP-dependent helicase. This enzyme is integral to modifying the structure of chromatin, which consists of DNA wrapped around histone proteins. By altering chromatin architecture, SMARCA1 facilitates or inhibits access to certain DNA regions, thereby regulating gene transcription.

- **Gene Regulation**: Through its action on chromatin, SMARCA1 influences various cellular processes including cell growth, differentiation, and response to signaling. Its role in gene regulation is crucial for normal cellular function and development.

### Expression and Activity

- **Tissue Specificity**: The expression of SMARCA1 varies across different tissues and cell types, reflecting its diverse roles in various cellular contexts.

- **Regulatory Mechanisms**: SMARCA1 activity is finely tuned by various regulatory mechanisms, ensuring that chromatin remodeling occurs in a controlled manner in response to cellular needs.

### Clinical Significance

- **Disease Association**: Alterations or mutations in the SMARCA1 gene have been implicated in certain genetic disorders and cancers. Understanding its function and regulation provides insights into these diseases and potential therapeutic approaches.

- **Research Interest**: SMARCA1 remains an area of active research, particularly in understanding its specific roles in different types of cells and how its dysfunction contributes to disease pathology.

### Conclusion

SMARCA1 is a pivotal component of the chromatin remodeling complex, essential for regulating gene expression and maintaining cellular homeostasis. Its diverse roles and clinical significance make it a key subject in the study of molecular biology and genetics.

## SMARCA4

SMARCA4, also known as BRG1, is a pivotal ATPase component of the SWI/SNF chromatin remodeling complex, instrumental in regulating gene expression by altering chromatin structure.

### Function in Chromatin Remodeling
#### Modulating Gene Accessibility
- SMARCA4 utilizes the energy from ATP hydrolysis to reposition nucleosomes, which controls the accessibility of transcriptional machinery to gene regulatory regions.

### Role in Gene Regulation
#### Partnering with Transcription Factors
- While SMARCA4 itself does not directly recognize specific DNA sequences, it is recruited to promoters and enhancers by transcription factors such as SOX4. The SOX4-SMARCA4 complex then collaborates to enhance chromatin accessibility, particularly at genes like TGFBR2, facilitating their expression.

### ChIP Experiments with SMARCA4 and SOX4
#### Capturing DNA-Protein Interactions
- In Chromatin Immunoprecipitation (ChIP) assays using antibodies specific to SMARCA4, DNA fragments bound by SMARCA4, including those at the promoters and enhancers of target genes, can be isolated. 
- Similarly, ChIP with SOX4-specific antibodies will precipitate DNA regions where SOX4 is bound, which may include the enhancer and promoter regions of genes such as TGFBR2.

<div align=center>
<img src="imgs/smarca.png">
</div

## SMART (Switching Mechanism at the 5′ end of RNA Template) technology
leverages the template-switching capability of certain reverse transcriptases (RTs) to capture full-length sequence information from RNA or DNA templates and incorporate adapter sequences during first-strand cDNA synthesis. The adapter sequences serve as primer-binding sites during subsequent rounds of PCR amplification. It has full-length template coverage, high-quality sequencing libraries.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/64097fd0b647151e7fe5ae47e03b89e3bee775c2/imgs/SMART-1-cDNA.png">
</div
.
        
.
## SMART (Switching Mechanism at the 5′ end of RNA Template) seq 2

Conventional cDNA construction methods usually result in an underrepresentation of the 5’ ends of cDNA.

1 - Single cell sorting: Smart-seq uses flow cytometry for cell sorting, with a maximum throughput of 96 cells at a time.

2 - Cell Lysis: Cell lysis is performed in cell lysis solution.

3 - Reverse transcription (one-strand synthesis): RNA with polyA tails (mainly mRNA) is reverse transcribed using Oligo(dT) primer. Since a mouse-derived reverse transcriptase ( Moloney Murine Leukemia Virus ) is used for reverse transcription, three Cs are added to the 3 end of the cDNA strand.

4 -Template replacement (second-strand synthesis): This step synthesizes the second strand of cDNA using TSO (template-switching oligo) primers, thus replacing the RNA complementary to the first-strand cDNA. TSO also carries the primers required for PCR (green part of the figure).

5 - PCR amplification: This step performs a light cDNA enrichment, and amplification of cDNA to the ng level is sufficient.

<b>For the following 6 and 7, please see Tn5 transposase (Tnp) for NGS sequencing library construction. They are the same.</b>

6 - Labeling: The DNA is interrupted using a modified, highly active Tn5 transposase while adding the junction to both ends of the cDNA. The DNA fragment after labeling is usually 200-600bp.

7 - PCR Enrichment & On-board Sequencing: After the last PCR amplification, the DNA is ready for onboard sequencing. Smart-seq2 uses the usual Illumina sequencer, as its library ends up being similar to the usual bulk RNA-seq.

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/5da88495006b561f1479e44ebf1276f96638c555/imgs/smart2.png">
</div
        
<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/e620aa051aa306e5c66ff75790b5e005988465d5/imgs/SMART-single-cell.png">
</div
.

..       
## SNV (Single nucleotide variant)
A SNV can be rare in one population but common in a different population. Sometimes SNVs are known as single nucleotide polymorphisms (SNPs), although SNV and SNPs are not interchangeable. To qualify as a SNP, the variant must be present in at least 1% of the population.

## Soft clipping and hard clipping

About soft clipping and hard clipping, it means that when query matching, some sequences are not matched completely, but soft clipping will keep the unmatched part afterwards, and hard clipping can remove the matched part completely. For example, when there is only part of the fragment is compared to the reference sequence when comparing. The difference is that a Soft Clip will eventually retain the corresponding sequence in the sequence that follows, while a Hard Clip will delete the fragment directly in the sequence that follows. BWA-MEM uses soft clipping CIGAR operation for supplementary alignments. By default, BWA-MEM uses soft clipping for the primary alignment and hard clipping for supplementary alignments.

Hard Clip exists with the intention of reducing the redundancy of BAM file sequences, for example, there is a read which can be compared to two places A, B. In place A, it is 60M90S, and in place B it is 60H90M, at this time a read actually already has the complete sequence information in position A, and the information in position B is actually redundant. So a marker form like Hard Clip can be introduced at location B, and it will be able to mark the sequence at location B as secondary.

## Splicing junctions

Key to defining the complexity of alternative splicing within a gene is the identification of splice junctions (SJs), which occur at exon-exon boundaries and are typically characterized in pairs representing both the donor site (5’ intron boundary to 3’ upstream exon boundary) and acceptor site (3’ intron boundary to 5’ downstream exon boundary).

## Stages and Treatments of Prostate Cancer

### Localized Prostate Cancer
- **Description**: Cancer confined to the prostate.
- **Treatment**: Surgery, radiation, or active surveillance. Early stage, high cure rate.

### Metastatic Castration-Sensitive Prostate Cancer (mCSPC)
- **Description**: Cancer has spread beyond the prostate but responds to androgen deprivation therapy (ADT).
- **Treatment**: ADT, often combined with chemotherapy (e.g., docetaxel) or novel hormone therapies (e.g., abiraterone). More serious than localized cancer but controllable with ADT.

### Castration-Resistant Prostate Cancer (CRPC)
- **Description**: Cancer progresses despite low testosterone levels from ADT.
- **Treatment**: Second-line hormonal therapies (e.g., enzalutamide, abiraterone), chemotherapy (e.g., docetaxel, cabazitaxel), immunotherapy, and targeted treatments. Most severe stage, requires diverse and potent treatments.

## Sticky ends
One strand is longer than the other (typically by at least a few nucleotides), such that the longer strand has bases which are left unpaired. The sticky ends, a.k.a. cohesive ends, have unpaired DNA nucleotides on either 5’- or 3’- strand, which are known as overhangs. These overhangs are most often generated by a staggered cut of restriction enzymes. Sticky ends are generally more desired in cloning technology where a DNA ligase is used to join two DNA fragments into one, because the yield and specificity of ligation using sticky ends is significantly higher that with blunt ends.

## Strand of Transcripts in IGV Genome Browser

In the Integrated Genome Viewer (IGV), the strand of transcripts can be determined by the directionality indicated by color coding. This feature allows researchers to quickly identify which strand of the DNA the transcript is being read from.

- **Right-clicking** on the feature track in IGV allows users to select the option to color by strand. This option applies color coding to the reads to indicate their orientation relative to the reference genome.
  
- **Red-colored reads** indicate that the sequencing reads are mapped in the same direction as the positive or sense strand of the DNA. This means that the read's 5' end aligns with the 5' end of the reference, and the read's 3' end aligns with the 3' end of the reference.

- **Blue-colored reads** signify that the sequencing reads are mapped in the opposite direction to the negative or antisense strand of the DNA. This means that the read's 5' end aligns with the 3' end of the reference, and the read's 3' end aligns with the 5' end of the reference.

By using this color-coding system, one can quickly ascertain not just the location of a transcript on the genome but also the directionality of its transcription.
  
## Structural variation  
Genomic structural variation is the variation in structure of an organisms chromosome. It consists of many kinds of variation in the genome of one species, and usually includes microscopic and submicroscopic types, such as deletions, duplications, copy-number variants, insertions, inversions and translocations. Originally, a structure variation affects a sequence length about 1kb to 3Mb, which is larger than SNPs and smaller than chromosome abnormality (though the definitions have some overlap). However, the operational range of structural variants has widened to include events > 50bp.The definition of structural variation does not imply anything about frequency or phenotypical effects. Many structural variants are associated with genetic diseases, however many are not.Recent research about SVs indicates that SVs are more difficult to detect than SNPs.

## Structural variation - Inversion

Why are inversions defined as the reverse complement and not just the reverse of the reference? If an inversion were just reversed then there would be 3' -> 3' bonds and 5' -> 5' bonds. That's why inversions are reverse complemented, you then maintain the normal 5'->3' direction.

## Structural variation - Novel sequence insertion

An insertion the sequence of which cannot be mapped to the reference genome.

## Structural variation - SURVIVIOR software merge translocation

Here is a translocation result from SVIM:

chr1	876903	svim.BND.3	N	[chr20:29368734[N	6	PASS	SVTYPE=BND;SUPPORT=5;STD_POS1=.;STD_POS2=.	GT:DP:AD	./.:.:.,.

Here is a translocation result from pbsv:

chr1	876903	pbsv.BND.chr1:876903-chr20:29368734	G	[chr20:29368734[G	.	PASS	SVTYPE=BND;CIPOS=0,0;MATEID=pbsv.BND.chr20:29368734-chr1:876903	GT:AD:DP	0/1:30,8:38

Here is merged result from SURVIVOR:

chr1	876903	svim.BND.3	N	N[chr20:29368734[	6	PASS	SUPP=2;SUPP_VEC=11;SVLEN=0;SVTYPE=TRA;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=29368734;CIPOS=0,0;CIEND=0,0;STRANDS=++	GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO	0/1:NA:28491831:0,0:++:.:TRA:pbsv.BND.chr1_876903-chr20_29368734:NA:NA:chr1_876903-chr20_29368734	./.:NA:28491831:0,0:++:6:TRA:svim.BND.3:NA:NA:chr1_876903-chr20_29368734

Here is a detailed explanation of each field in the merged VCF line:

1. chr1: This is the name of the chromosome on which the variant is located. In this case, the variant is located on chromosome 1.
2. 876903: This is the position of the variant on the chromosome, in base pairs. In this case, the variant is located at position 876903 on chromosome 1.
3. svim.BND.3: This is the ID of the variant as called by the SVIM software. It is a breakend variant with ID 3.
4. N: This is the reference allele at the position of the variant on chromosome 1.
5. N[chr20:29368734[G: This is the alternate allele at the position of the variant on chromosome 1. It indicates that the variant involves a breakend on chromosome 1 and a breakpoint on chromosome 20 at position 29368734.
6. 6: This is the quality score of the variant.
7. PASS: This field indicates whether the variant passed filtering. In this case, the variant passed filtering.
8. SUPP=2;SUPP_VEC=11;SVLEN=0;SVTYPE=TRA;SVMETHOD=SURVIVOR1.0.7;CHR2=chr20;END=29368734;CIPOS=0,0;CIEND=0,0;STRANDS=++: This field contains additional information about the variant. The SUPP field indicates the number of supporting reads for the variant. The SUPP_VEC field indicates the number of supporting read pairs.
9. The SVLEN field indicates the length of the variant, in base pairs. In this case, the value is 0, which indicates that the variant is 0 base pairs in length.
10. The SVTYPE field indicates the type of structural variant. In this case, the value is TRA, which indicates that the variant is a translocation.
11. The SVMETHOD field indicates the software that was used to call the variant. In this case, the value is SURVIVOR1.0.7, which indicates that the variant was called using the SURVIVOR software.
12. The CHR2 field indicates the second chromosome involved in the variant, if applicable. In this case, the value is chr20, which indicates that chromosome 20 is also involved in the variant.
13. The END field indicates the end position of the variant on the chromosome. In this case, the value is 29368734, which indicates that the variant ends at position 29368734 on chromosome 20.
14. The CIPOS and CIEND fields indicate the confidence intervals for the start and end positions of the variant, respectively. In this case, the values are 0,0, which indicates that the confidence intervals for the start and end positions are 0 base pairs.
15. The STRANDS field indicates the strands on which the variant was observed. In this case, the value is ++, which indicates that the variant was observed on both strands.

GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO	0/1:NA:28491831:0,0:++:.:TRA:pbsv.BND.chr1_876903-chr20_29368734:NA:NA:chr1_876903-chr20_29368734	./.:NA:28491831:0,0:++:6:TRA:svim.BND.3:NA:NA:chr1_876903-chr20_29368734

the explanation for the genotype, phase set, and other information for the two samples:

For the first sample:

+ The ID field indicates the ID of the variant. In this case, the value is pbsv.BND.chr1_876903-chr20_29368734, which is the ID of the variant as called by the PBSV software.
+ The RAL field indicates the reference alleles of the variant. In this case, the value is NA, which indicates that the reference alleles are not available.
+ The AAL field indicates the alternate alleles of the variant. In this case, the value is NA, which indicates that the alternate alleles are not available.
+ The CO field indicates the components of the variant. In this case, the value is chr1_876903-chr20_29368734, which indicates that the variant is a breakend involving chromosome 1 and chromosome 20 and spanning from position 876903 on chromosome 1 to position 29368734 on chromosome 20.

For the second sample:

+ The GT field indicates the genotype of the sample. In this case, the value is ./, which indicates that the genotype of the sample is missing.
+ The PSV field indicates the phase set of the sample. In this case, the value is NA, which indicates that the phase set is not available.
+ The LN field indicates the length of the variant. In this case, the value is 28491831, which indicates that the variant is 28491831 base pairs in length.
+ The DR field indicates the number of supporting reads for the variant. In this case, the value is 0,0, which indicates that there are no supporting reads for the variant.
+ The ST field indicates the strands on which the variant was observed. In this case, the value is ++, which indicates that the variant was observed on both strands.
The QV field indicates the quality score of the variant. In this case, the value is 6, which is the quality score of the variant as called by the SVIM software.
+ The TY field indicates the type of variant. In this case, the value is TRA, which indicates that the variant is a translocation.
+ The ID field indicates the ID of the variant. In this case, the value is svim.BND.3, which is the ID of the variant as called by the SVIM software.
+ The RAL field indicates the reference alleles of the variant. In this case, the value is NA, which indicates that the reference alleles are not available.
+ The AAL field indicates the alternate alleles of the variant. In this case, the value is NA, which indicates that the alternate alleles are not available.
+ The CO field indicates the components of the variant. In this case, the value is chr1_876903-chr20_29368734, which indicates that the variant is a breakend involving chromosome 1 and chromosome 20 and spanning from position 876903 on chromosome 1 to position 29368734 on chromosome 20.

## Structural variation detection algorithm

<div align=center>
<img src="imgs/SV_algo.jpeg">
</div

## Super enhancer

Super-enhancers (SEs) are collections of closely spaced genomic regions that exhibit enhancer hallmarks. These enhancers are usually stitched together by constituent enhancers, but there can be gaps of up to 12.5 kb between them.

## Supplementary Alignment

This is the SA tag. A chimeric reads but not a representative reads. 

## SWI/SNF Complex 

### Core Description
The SWI/SNF complex is a type of chromatin modifier. It can open the chromation. When transcription factors require assistance, they recruit the SWI/SNF complex, which then binds with the transcription factors. This interaction leads to the unraveling of chromatin, exposing enhancers or promoters, thereby facilitating the binding of transcription factors to enhancers or promoters, crucial for effective gene regulation.

### Composition
The SWI/SNF complex is composed of multiple subunits, which can vary among different cell types and species. Key components include:

- **SMARCA4/BRG1 or SMARCA2/BRM**: Both subunits have ATPase activity and serve as the main driving forces of the complex. They play a central role in different SWI/SNF complexes but usually do not coexist in the same complex.
- **SMARCB1/SNF5, SMARCC1/BAF155, SMARCC2/BAF170**: These subunits provide structural support and are involved in the stability and assembly of the complex.
- **SMARCE1/BAF57, SMARCD1/BAF60, SMARCD2/BAF60B, SMARCD3/BAF60C**: These components also participate in the assembly and functional regulation of the complex.
- **Other Components**: Including specific subunits like ARID1A/BAF250A and ARID1B/BAF250B, which impart specific functions and tissue specificity to the complex.


### (a) Attraction of SWI/SNF through a transcriptional activator

The transcriptional activator directly binds to the promoter region upstream of the inflammatory gene. This allows the recruitment of the SWI/SNF complex to the same promoter region. Although SWI/SNF does not directly interact with DNA, its chromatin remodeling activity enables an open chromatin structure, facilitating transcription factor binding and gene activation.

## SWI/SNF binding to proinflammatory gene promoters

### (a) Attraction of SWI/SNF through a transcriptional activator

The transcriptional activator binds to the enhancer region of the inflammatory gene. This allows the recruitment of the SWI/SNF complex to the promoter region upstream of the gene. Although SWI/SNF does not directly interact with DNA, its chromatin remodeling activity enables an open chromatin structure, facilitating transcription factor binding and gene activation.

### (b) Attraction of SWI/SNF due to modifications of the N-ends of histones

Modifications (such as acetylation) of histone N-terminal tails can help recruit SWI/SNF. These modifications loosen chromatin structure at target gene promoters. SWI/SNF recognizes and binds to modified histones through distinct domains. Its ATP-dependent remodeling further opens up chromatin, enabling gene transcription.

### (c) Attraction of SWI/SNF through the MEDIATOR complex

The Mediator complex provides a bridge to recruit SWI/SNF to target promoters. Transcriptional activators bind to enhancers and help recruit Mediator, which then interacts with SWI/SNF. SWI/SNF is thus brought to promoter upstream regions where it remodels chromatin structure to allow assembly of the transcription machinery.

### (d) Attraction of SWI/SNF through long non-coding RNA (lncRNA)

Some lncRNAs can interact with transcriptional activators and Mediator. Through this interaction, lncRNAs help bring SWI/SNF to inflammatory gene promoters by facilitating activator and Mediator binding. SWI/SNF then remodels local chromatin structure, removes the nucleosome barrier, and activates inflammatory gene expression.

## T4 DNA Ligase Summary

T4 DNA Ligase is a critical enzyme for DNA and RNA manipulation, known for its ability to ligate double-stranded nucleic acids and repair nicks in complex structures.

### Capabilities

- **Universal Ligation**: T4 DNA Ligase is capable of facilitating ligation between any two of the following molecule types: dsDNA, dsRNA, and DNA/RNA hybrids. 

- **Nick Repair**: Additionally, T4 DNA Ligase excels at repairing nicks within the double-stranded regions of DNA, RNA, and DNA/RNA hybrids, thereby maintaining or restoring molecular integrity.

### Requirements

- **Acts only on double-stranded nucleic acids.**
- **Efficiently ligates blunt-ended molecules. The only one ligates blunt end**
- **Requires a 5' phosphate group for ligation.**

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
This procedure does not include whole genom amplification, could not used ditectly for PCR.

A schematic overview of the SMRT-Tag approach. Hairpin adaptor-loaded triple mutant Tn5 transposase is loaded and used to fragment DNA into 2 – 10 kilobase (kb) fragments. After removing Tn5 transposase, an optimized gap repair is used to fill the resulting 9 bp gaps on either side of the molecule, and an exonuclease treatment is used to purify repaired covalently closed templates, which are sequenced on the PacBio Sequel II instrument. 
  
Crucially, combining tagmentation with optimized gap repair allowed the streamlined creation of PacBio libraries from 80 – 100 ng DNA (a minimum of 15,000 human cell equivalents) compared to current protocols that require > 0.5 – 5 μg DNA (a minimum of ∼200,000 human cell equivalents). 

<div align=center>
<img src="https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/3cc296e4e1b223f9e3d836fbd5f6fb0b54ca2c3e/imgs/pc.jpg">
</div

 .
 
## Tn5 transposase (Tnp) for PacBio sequencing library construction (PCR-based) a. k. a. SNOOTH-SEQ  
Individual cells were collected with a microcapillary connected to a mouth pipette and washed by transferring them into droplets of 1 mg/mL phosphate-buffered saline-bovine serum albumin for three times before lysis. The 2.5-μL lysis reaction consists of 0.25-μL 100mM Tris-EDTA (1M Tris + 0.1M EDTA), 0.125μL Qiagen protease, 0.075μL 10% triton X-100, 0.05μL1MKCL,and2μLH2O. The cell lysis was carried out at 50°C for 3 h to digest the proteins binding on the gDNA and then 70°C for 30 min to inactivate the protease. After that, a 7.5-μL tagmentation mixture including 2 μL 5×TAPS_PEG8K (50 mM TAPS-NaOH (or KOH), pH 8.3 (RT), 25 mM MgCl2, 40% PEG8K), and 1μL0.2ng/μL adaptor conjuncted Tn5 enzyme (Vazyme, Cat. S601-01) was added into each cell lysate. The tagmentation reaction was carried out at 55°C for 10min, followed by adding 2.5-μL0.2%SDSand standing at room temperature for 5min to stop tagmentation, releasing the fragmented gDNA. Then, strand displacement of the Tn5 adaptors and amplification of the fragmented gDNA was carried out using 0.025U/μLTksGflexDNA Polymerase (TAKARA, Cat. R060B), 560nM I5 PCR primer which containing 16 bp cell barcode (5′ AATGATACGGCGACCACCGAGATCTNNNNNNNNNNNNNNN NTCGTCGGCAGCGTC3′). The PCR program was 72°C 3min, 98°C 1 min, and then 20 cycles of 98°C for 15s, 60°C for 30s, and 68°C for 5min. After that, gDNA amplicons using different barcode primers were pooled together and purified with 0.4 volume of Ampure PB beads (Pacific Biosciences Ref. No. 100-265-900) for twice. These purified amplicons were quantifiedusingQubit,andabout1μg amplified products was used to construct libraries for Pacbio sequencing using SMRTbell Template Prep Kit v.1.0-SPv3 (Pacific Biosciences Ref. No. 100-991-900).

Because of the transposon approach, adaptor has been added, each sequence is ended with part of the sequence of adaptor, so there is no need for random primers to amplify the whole genome, but only the I5 primer that specifically identifies the adaptor is needed, except that barcode needs to be added to label multiple cells to facilitate the difference after hybrid sequencing, saving cost and improving throughput

## Transfection

Transfection introduces nucleic acids into cells to study gene function, essential for RNAi and CRISPR experiments. They have chemical or physical methods.

### Forward Transfection
Forward transfection involves plating cells first and then adding a mixture of siRNA (or DNA) and transfection reagent for 10 cm and 15 cm dishes 

### Reverse Transfection
Adds the nucleic acid and transfection reagent mixture to culture wells before the cells, good for small-volume cultures like 96-well plates.

## Transcription Factors: Key Regulators of Gene Expression

### Overview
Transcription factors are proteins that play a critical role in regulating gene transcription. They bind to specific DNA sequences, particularly in the promoter or enhancer regions of genes, to control the recruitment of RNA polymerase and, consequently, the activity of gene expression.

### Binding to Promoters
#### Direct Activation
Transcription factors can bind directly to promoter regions located upstream of genes. This binding helps position RNA polymerase II and initiates the transcription of the gene.

### Interacting with Enhancers
#### Long-Distance Regulation
Transcription factors can also bind to enhancer regions, which are DNA sequences that can boost the transcription of target genes from a distance. Through the mediation of complexes like the Mediator complex, transcription factors can influence promoters and regulate gene transcription even from afar.

### Overcoming Chromatin Inaccessibility
#### Enhancer and Promoter Accessibility
When chromatin is tightly packed and gene loci are inaccessible, transcription factors can interact with (binding) chromatin modifiers, such as the SWI/SNF complex, to modulate the chromatin state. This modification not only makes promoter regions accessible for transcription initiation but also exposes enhancer regions. Transcription factors may first bind to these enhancers and, through spatial reorganization of the chromatin, facilitate the interaction between enhancers and promoters, thereby regulating gene transcription. This interplay between transcription factors, enhancers, and chromatin state is a key step in the control of gene expression.

### Histone Modification Interactions
#### Modulating Chromatin Accessibility
In addition to altering promoter accessibility, transcription factors can interact with histone modification enzymes like histone acetyltransferases (HATs) and histone deacetylases (HDACs). By modifying the acetylation state of histones, they influence the structure and accessibility of chromatin. These modifications can lead to the exposure of enhancer regions, allowing transcription factors to bind and regulate gene expression through the activation of promoters.

These diverse mechanisms ensure that transcription factors can effectively regulate gene expression across different cellular environments and conditions. By employing these complex regulatory modes, cells can precisely control gene activity in response to changing internal and external environments.

## Transcription factor binding site (TFBS)

TFs recruit additional proteins to either activate or repress gene expression. Some TFs bind to a DNA **promoter** sequence near the transcription start site and help form the transcription initiation complex. Other TFs bind to regulatory sequences, such as **enhancer** sequences, and can either stimulate or repress transcription of the related gene.
 
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
  
## VCF format
VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.

<b>This most important thing is that for deletion, dup, or insertion, the base at the ref is not included. </b>

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









