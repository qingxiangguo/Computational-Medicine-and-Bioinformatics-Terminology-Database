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
## SAM (file format)
![image](https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/a4ddccceb15fd1d1ee05ae6bb0183febb48feae4/imgs/1.png)

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

<b>The third column</b>: RNAME-- Reference sequence NAME of the reference genome corresponding to the name of the chromosome or contig, scaffold. If *, it means there is no result on the comparison.

<b>The fourth column</b>: POS-- 1-based leftmost mapping POSition fragment leftmost mapping to the reference genome, counting from 1. If it is 0, it means there is no match.

<b>The fifth column</b>: MAPQ-- MAPping Quality mapping score, -10 log10 Pr{mapping position is wrong}- logarithmic transformation of the probability that the mapping position is wrong with a base of ten. According to this formula, we can also know that the larger the value, the lower the probability that the pairing is wrong.

<b>The sixth column</b>:  CIGAR - CIGAR string stands for Concise Idiosyncratic Gapped Alignment Report. First of all, there are the following letters that represent different meanings. Most of the meanings are easy to understand, such as sequence matching, insertions and deletions, and skipping the fragment. The document mentions that Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.  
![image](https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/491ce3bf5a3dc72277520ada9c45236498586c4b/imgs/2.png)  
For example, 251M, means the alignment of the full match 

30M3D126M3D58M37S, it means that the reads are 30bp match + 3bp missing + 126bp match + 3bp missing + 58bp match + 27bp soft clipping.  

About soft clipping and hard clipping, it means that when query matching, some sequences are not matched completely, but soft clipping will keep the unmatched part afterwards, and hard clipping can remove the matched part completely.

![image](https://github.com/qingxiangguo/Computational-Medicine-and-Bioinformatics-Terminology-Database/blob/491ce3bf5a3dc72277520ada9c45236498586c4b/imgs/3.png)  

For example, in this example above, you can see that only part of the fragment is compared to the reference sequence when comparing. The difference is that a Soft Clip will eventually retain the corresponding sequence in the sequence that follows, while a Hard Clip will delete the fragment directly in the sequence that follows.   

"Hard Clip exists with the intention of reducing the redundancy of BAM file sequences, for example, there is a read which can be compared to two places A, B. In place A, it is 60M90S, and in place B it is 60H90M, at this time a read actually already has the complete sequence information in position A, and the information in position B is actually redundant. So a marker form like Hard Clip can be introduced at location B, and it will be able to mark the sequence at location B as secondary."

If there is a sequence that comes from a mature mRNA, if this sequence skips exactly 200bp of intron in the middle and 75bp mapped to exon before and after, what should the CIGAR value of this sequence be written?CIGAR:75M200N75M

<b>The seventh column</b>: read2's chromosome name on the reference sequence, if not available use "*", same as "="

<b>The eighth column</b>:: PNEXT: position of read2 on the reference sequence

<b>The ninth column</b>:: TLEN: length of inserted fragment  

TLEN - observed Template LENgth If all read segments are mapped to the corresponding reference sequence, the absolute value of TLEN is equal to the distance (end-start+1) between the mapped end of the template sequence and the mapped start of the template sequence (including both ends). It should be noted that the bases on the comparison do not include sof-clipped bases. If the read segment is compared to the start of the leftmost segment of the template, the TLEN field is positive, and if the comparison is to the start of the rightmost segment, it is actually an antisense strand and the TLEN field is negative. If the starting position of the comparison is the same at both ends, then any positive or negative number is assigned. If there is only a single chain, the value is 0. And the positive and negative numbers of any intermediate segments are undefined.







