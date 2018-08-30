
# *Please check aloft.gersteinlab.org for packaged downloads and updates*

This software is freely licensed under the Creative Commons license 
(Attribution-NonCommerical). The main aspects of this license are that: The 
work can be made available for non-commercial use. Derivatives can be made of 
the work. Derivatives do not have to be made available under the same terms 
that they were first used, and we should be cited.

###############################################################################

ALoFT Annotation of Loss-of-Function Transcripts

README DOCUMENTATION

Gerstein Lab
Molecular Biophysics & Biochemistry

Version 1.0                           
Released 2013               
Yale University

Contact: Dr. Suganthi Balasubramanian (suganthi.bala@yale.edu), Mayur Pawashe (mpawashe@gmail.com), Jeremy Liu (jeremy.liu@yale.edu)

###############################################################################

Table of Contents
A. Preface
B. Description
C. System Requirements
D. Usage
E. Input Files
F. Output Files
G. ALoFT Features in VCF Output
H. ALoFT Features in Tabbed Delineated Output
I. Command Line Options

A. Preface
Installation instructions can be found in INSTALL.
Quick Start instructions can be found in QUICK START.

B. Description
ALoFT takes as input an VCF (variant call format) file, runs it through the 
variant annotation tool (VAT) tools snpMapper and indelMapper, and runs 
the sorted VAT output through aloft, which calculates other variant-
specific features that give functional, evolutionary, mismapping, and other
information. ALoFT only calculates these parameters for frameshift indels, 
loss-of-function (LoF) SNPs, and for SNPs located in splice sites. ALoFT 
is split into separate modules based on features.

C. System Requirements
Python 2.7.x or Python 3.x
Linux (64-bit) or OSX 10.6 (64-bit) and up

D. Input Files (Required):
1) VCF file containing unannotated variants, passed in using the --vcf option.
Alternatively, an annotated VCF file output by a separate manual run of VAT 
should be passed in with the --vat option. See options section below.

2) Reference files (defaults already in place after installation). 
For complete list, see section below regarding options and default values.


E. Output Files:
An output directory may be specified with the --output option.
Otherwise the output files will be written in the ./aloft_output directory.

The following files will be created in this output directory:

1) A vcf file that VAT generates if --vcf is supplied called vat_output.vcf.
This intermediate file is then run through aloft to produce the following 
three files. The variants in this file are listed in order, starting with
the lowest position on chromosome 1.

2) Output VCF file named <input_file_name>.aloft.vcf
Contains a culled list of formatted information calculated by VAT and 
features calculated by ALoFT in Variant Call Format for both lof variants
and spice site variants.

3) Tab-delimited file named <input_file_name>.aloft.lof
Contains the subset of VAT output associated with LoF SNPs and frameshift 
indels. All information from the VAT file is included, plus ALL features 
calculated by aloft for LoF variants.
Variants in the input file that intersect coding exons are included.
For variants with multiple transcripts, the variant data is calculated and 
output on different lines for each transcript.

4) Tab-delimited file named <input_file_name>.aloft.splice
Contains the subset of VAT output associated with splice site SNPs.  
All information from the VAT file is included, plus ALL features 
calculated by aloft for splice site SNPs.  Variants with multiple 
transcripts are output on multiple rows, one for each combination of 
alternate allele and affected transcript.
Variants in the input file that intersect splice sites are included.


F. Usage
ALoFT can be invoked as follows:
 $ cd aloft
 $ ./aloft --vcf=/path/to/file --data=path/to/data/dir --output=path/to/dir [--option4=arg1]...

Example:
 $ ./aloft --vcf=vcf_file.vcf	(unannotated vcf file)
OR 
 $ ./aloft --vat=vat_file.vcf	(VAT annotated vcf file)
OR 
 $ ./aloft --vcf=vcf_file.vcf --output=/path/to/directory   
   (VAT annotated vcf file and custom output destination)


G. ALoFT Features in VCF Output
ALoFT retains input file VCF metaheader and variant information and details.
ALoFT appends a subset of the features listed in Tabbed Delineated Output
into vcf form.

1) The following features are listed for all variants (lof and splice):
- AA, Ancestral, AF, AMR_AF, ASN_AF, AFR_AF, EUR_AF, VT, SNPSOURCE, AC, AF,
AN AVGPOST, ERATE, LDAF, RSQ, THETA, VA are annotations from VAT.
These should be listed and described in the first portion of the VCF metadata.
- Ancestral: Determines whether ancestral allele is the same as reference
- GERPscore: Gives the GERP score associated with the variant position
- SegDup: Gives the number of segmental duplications associated with the 
variant position.
- 1000GPhase1(|_AF|_ASN_AF|_AFR_AF|_EUR_AF): 1000Genomes Phase 1 allele freqs
	blank: Yes or No if associated 1000 Genomes allele frequencies
	AF: overall allele frequences
	ASN_AF: asian subset allele frequences
	AFR_AF: african subset allele frequences
	EUR_AF: european subset allele frequencies
- ESP6500(|_AAF): ESP6500 allele frequences
	blank: Yes or No if associated ESP6500 allele frequencies
	AAF: ancestral allele frequencies
- GERPelement: YES if variant/transcript has associated GERPelement and
no otherwise.
- exoncounts: Gives the number of exons in the transcript and the number
of exons truncated by the premature stop, inclusive. For splice site variants
the number of exons truncated is replaced with the string "NA".

2) These features are listed for lof variants, in addition to those in 1):
- nearstart: YES if variant is within the first 5% of the coding sequence (i.e beginning of the coding sequence), NO otherwise.
- nearend: YES if variant in the last 5% of the coding sequence (i.e. towards the end of the coding sequence), NO otherwise.
- canonical: YES if the 5' flanking splice site and the 3' flanking splice 
site are both canonical, NO otherwise.
- XX/XX: <5' flanking splice site>:<3' flanking splice site of the exon that
the variant intersect, as reported in the reference genome.
- lofposition: calculated one indexed coding sequence position in which the
premature stop occurs in the alternate sequence.
- nmd: YES if the premature stop leads to nonsense mediated decay, NO if not.
- lof_anc: YES if ancestral allele alternate sequence leads to premature stop
and NO otherwise.
- heavilyduplicated: YES if the number of duplicated regions (that the variant
is in) is high, NO otherwise
- disorder_prediction: Gives the percentage of residues that are disordered in
the reference sequence and the percetageof residues that are disordered in the
truncated alternate sequence (after the premature stop) as <%>:<%>. If the
transcript is not associated with disordered regions, "." is output.
- PF: PFAM protein domains. YES if variant intersects region, NO if no 
intersection and NA if no regions exist for the particular transcript.
- SSF: SSF protein domains. YES if variant intersects region, NO if no 
intersection and NA if no regions exist for the particular transcript.
- SM: SM protein domains. YES if variant intersects region, NO if no 
intersection and NA if no regions exist for the particular transcript.
- Tmhmm: Transmembrane helix domains. YES if variant intersects region, NO if
no intersection and NA if no regions exist for the particular transcript.
- Sigp: Signal peptide domains. YES if variant intersects region, NO if no 
intersection and NA if no regions exist for the particular transcript.
- PTM: Post translational modifications. YES if variant intersects any 
post translational modification regions, NO if no intersections and NA 
if no regions exist for the particular transcript.

3) These features are listed for splice variants, in addition to those in 2):
- XX/XX: <acceptor_site>:<donor_site> splice sites in reference genome. 
- is_canonical: YES if variant intersecting splice site is canonical, NO if not
- other_canonical: YES if other splice site in the intron, not the splice site
intersected by the variant is canonical, NO otherwise.
- intron_length: Gives the length of the intron that the varint intersects the
splice site.
- small_intron: YES if the length of intron is less than 15bp, NO otherwise.
- heavily_duplicated: YES if the variant region is heavily duplicated, NO
otherwise.
- lof_anc: YES if the ancestral alternate sequence leads to a premature stop
and NO otherwise.
- alternate_acceptor_site: YES if there are potential neighboring splice sites
that could replace a malfunctioning splice site at the variant location and
NO otherwise. This is the NAGNAG case.


H. ALoFT Features in Tabbed Delineated Output
1) ALOFT calculates the following features for all variants:
- VAT Features: includes all features from VAT snpMapper and indelMapper
This includes allele frequencies, variant type, etc. This is the "details"
column in the tabbed delineated output and the first part of the details
section of each transcript in the vcf output.
- partial/full: full if all transcripts of the affected gene are affected, 
partial otherwise
- transcript_length: length of affected transcript in nucleotides
- longest_transcript: YES if transcript is the longest transcript affected by 
variant, NO otherwise
- shortest path to recessive gene: minimum length of a shortest path to a 
recessive gene in protein interaction network
- recessive neighbors: gives the number of recessive genes that are directly
connected to the LOF-containing gene in a protein-protein interaction network
- shortest path to dominant gene: minimum of length of a shortest path to a 
dominant gene in the protein interaction network
- dominant neighbors: gives the number of dominant genes that are directly
connected to the LOF-contianing gene in a protein-protein interaction network
- GERP score: associated GERP score of the variant position
- GERP element: associated GERP element of the variant position
- GERP rejection: associated GERP rejection score of the variant position
- exon counts: number of exons associated with the variant transcript
- Segmental duplications: Gives the position of associated segdups as a 
bracketed list, or a period if none exist.
- 1000GPhase1(|_AF|_ASN_AF|_AFR_AF|_EUR_AF): 1000Genomes Phase 1 allele freqs
	blank: Yes or No if associated 1000 Genomes allele frequencies
	AF: overall allele frequences
	ASN_AF: asian subset allele frequences
	AFR_AF: african subset allele frequences
	EUR_AF: european subset allele frequencies
- ESP6500(|_AAF): ESP6500 allele frequences
	blank: Yes or No if associated ESP6500 allele frequencies
	AAF: ancestral allele frequencies
- Number of pseudogenes associated to transcript: 
number of pseudogenes or a period if none.
- Number of paralogs associated to gene: number of paralogs or a period if needed
- dN/dS (macaque): Evolutionary score in comparison to macaque species
- dN/dS (mouse): Evolutionary score in comparison to mouse species

2) ALoFT calculates the following features for lof snps and indels:
- is single coding exon?: YES if the variant intersects a transcript with 
only one coding exon, NO if the variant does not.
- indel position in CDS: gives the one indexed position of the indel in 
the coding sequence
- stop position in CDS: gives the one indexed position of the premature stop
in the coding sequence
- causes NMD?: YES if the variant causes nonsense mediated decay, calculated 
by default 50 base pair proximity of the premature stop to the last coding 
exon. NO if the variant does not lead to nonsense mediated decay.
- 5' flanking splice site: Gives the upstream 5' splice site of the exon that
the variant intersects
- 3' flanking splice site: Gives the downstream 3' splice site of the exon 
that the variant intersects
- canonical?: YES if the 5' flanking splice site is 'AG' and the 3' flanking 
splice site is 'GT', NO otherwise if neither of the splice sites matches 
'AG' and 'GT', respectively.
- Number of failed filters: number of call filters failed associated with the 
loss of function variant
- filters failed: list of the call filters failed associated with the variant
	heavily_duplicated: if many segmental duplications exist
	lof_anc: the ancestral allele leads to loss of function
	near_start: variant is in the first coding exon of the transcript
- ancestral allele: Gives the nucleotide at the variant position in
the ancestral reference genome
- Disorder prediction: Gives the percentage of disordered residues in the 
translated nucleotide sequence. Also gives the percentage of disordered
residues in the translated nucleotide sequence after the truncation caused 
by a premature stop. Or if a variant does not have disorder regions.

** Protein Families & Post Translational Modifications **
For the following features, 
	the region_id:count is output if variant intersects feature region
	NO_<feature> is output if variant does not intersect feature regions
	NA_<feature> is output if variant's transcript has no feature regions
- PF and PFtrunacated: Determines whether the variant intersects or truncates 
a PFAM protein segment.
- SSF and SSFtruncated: Determines whether the variant intersects or 
truncates a SFF protein segment.
- SM and SMtruncated: Determines whether the variant intersects or truncates 
a SM protein segment.
- Tmhmm and Tmhmmtruncationed: Determines whether the variant intersects or 
truncates a Tmhmm protein segment.
- Sigp and Sigptruncated: Determines whether the variant intersects or
truncates a Sigp protein segment.
- ACETYLATION(truncated), METHYLATION(truncated), PHOSPHORYLATION(truncated):
Post translational modifications determined from Phophosite. 
Determines whether the premature stop variant intersects or truncates 
post translational modification sites. 
Site types include acetylation, (mono/di/tri)methylation, O-GlcNAc, 
phophorylation, sumolyation, and ubiquitination. 

3) For splice SNPs only:
- Donor: Gives the nucleotide sequence of the donor splice site
- Acceptor: Gives the nucleotide sequence of the acceptor splice site
- SNP in canonical site?: Determines whether the splice site that the SNP
intersects is canonical (YES/NO)
- Other splice site canonical?: Determines whether the other splice site,
that the SNP does not intersect, is canonical (YES/NO)
- SNP location: Determines whether the snp intersects the donor or acceptor
- Alt donor: Gives the nucleotide seqeunce of the donor splice site after
the SNP change has been made
- Alt acceptor: Gives the nucleotide sequence of the acceptor splice site 
after the SNP change has been made
- NAGNAG positions: The NAGNAG case. Determines possible nearby canonical 
splice sites to the SNP location. Alternative splice sites.
- Intron length: Gives the length of the intron bracked by the donor and
acceptor splice site.
- Number filters failed: number of call filters failed associated with the splice
site variant
- Filters failed: list of the call filters failed associated with the variant
	ref_noncanonical: reference splice site, that the variant intersects, 
		is noncanonical
	alt_noncanonical: alternate splice site, that the variant intersects, 
		is noncanonical
	other_noncanonical: other reference splice site, that the variant does
		not intersect, is noncanonical
	heavily_duplicated: variant flagged if the variant intersects a 
		heavily duplicated region
	short_intron: variant flagged if the intron is shorter than 15bp


I. Options
aloft recognizes the following options for altering input and reference
files (default values given):
ALOFT will come packaged with most of the necessary reference files.

--version
Will output the ALoFT version number.

--vcf=""
Specifies path to VCF input file.  Set to empty string by default.  If none 
specified, ALoFT will try to skip VAT and run directly on the file 
given to the --vat option. This or --vat option is needed for proper execution.

--vat=aloft_output/vat_output.vcf
Specifies path to VAT output file to run aloft on.

--cache=cache/
Specifies path to directory containing cache of GERP score information and 
protein-protein interaction information.  
Directory will be created if it doesn't already exist.

--nmd_threshold=50
Distance from premature stop to last exon-exon junction; used to predict
NMD. Default distance is 50bp.

--output=aloft_output/
Specifies path to tabbed output files and VCF file from ALoFT.

—data=data/
Specifies path to data directory containing a data.txt file and other data dependencies.
data.txt contains paths to all data files that ALoFT requires.
See data/data.txt bundled with ALoFT for more information on these files.

--verbose
Will run ALoFT in verbose mode.
