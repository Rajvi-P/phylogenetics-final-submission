--------------------------------------------------------------------------------------------------------------------------------------------
## Background:

- The Yersinia pestis bacterium causes the plague. Its strains can be categorized into 3 biovars (strains that are similar genetically but differ physiologically) that correlate to 3 distinct waves of the plague. The oldest biovar is Antiqua, causing the first wave of plague in Africa. The second wave of plague was caused by the biovar Medievalis in Asia and Europe. The biovar Orientalis defines the third wave of plague that was transmitted to the Americas and to Australia. 

- An important genetic element that has contributed to the pathogenicity of this Yersinia pestis is the pMT1 plasmid. This plasmid has been thought to aid in deep tissue invasion (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC108724/).  However, the pMT1 plasmid's evolutionary history has not been explored, making the goal to create a phylogenetic tree of the pMT1 plasmid so that we can learn how pMT1 evolved to cause plague worldwide. 



--------------------------------------------------------------------------------------------------------------------------------------------

## Goal:
**Understand the evolutionary history of Yersinia pestis based on pMT plasmid evolution**

--------------------------------------------------------------------------------------------------------------------------------------------
## Step 0: Create an organized directory and learn how to install software

**Completed - Organizing folders**
- Created a "phylogenetics-final" folder. Inside of the "phylogenetics-final" folder, created (1) "notebook-log.md" file, which explained all steps taken; (2) "data" folder, which
contained all collected data, softwares, and input and output files; (3) "README.md" file, which provided a quick summary of the project.

**Completed - Software downloading training**
- First, install the software from the official website. Move the software to the folder containing your data. Double click on the software file to untar. Make sure a path is created to the software using the following code on terminal:
- _Code:_
```
nano ~/.rajvipatel

#now add the following lines:

export PATH="/opt/homebrew/bin:$PATH"

export PATH="/Users/rajvipatel/phylogenetics-final/data/softwarename:$PATH"

#close out via ^X and write the following:

source ~/.rajvipatel
```


## Step 1: Collect data

**Completed - Collecting data**

- Several whole genome nucleotide sequences of various pMT1 plasmid strains were gathered from GenBank and NCBI's RefSeq. Searches of the strains on Joint Genome Institute produced scientific study citations that mentioned the origin of some strains. Strains found to be having unknown origins or made in a lab or were removed. The remaining 13 sequences of data are stored separately in the "raw-data" subfolder found within the "data" folder.


## Raw Sequence Information

**Yersinia pestis Antiqua plasmid pMT, complete sequence**

- GenBank: CP000309.1

- Info: Isolated from a human infection in Africa (Republic of Congo in 1965)

- Sources that provide strain information: (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1482938/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=637000350)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/CP000309.1?report=fasta



**Yersinia pestis Nepal516 plasmid pMT, complete sequence**

- GenBank: CP000306.1

- Info: Isolated from a human infection in Nepal (possibly from a 1967 outbreak of pneumonic plague)

- Sources that provide strain information: (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1482938/, genbank=refseq: https://www.ncbi.nlm.nih.gov/assembly/GCA_000013805.1, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=645058759)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/CP000306.1?report=fasta



**Yersinia pestis biovar Microtus str. 91001 plasmid pMT1, complete sequence**

- GenBank: AE017045.1

- Info: Isolated from Microtus brandti in Inner Mongolia, China

- Sources that provide strain information: (https://pubmed.ncbi.nlm.nih.gov/15368893/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=637000354) 

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/AE017045.1?report=fasta



**Yersinia pestis biovar Medievalis str. Harbin 35 plasmid pMT, complete sequence**

- NCBI Reference Sequence: NC_017266.1

- Info: Isolated from a human in China

- Sources that provide strain information: (https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA30505, refseq=genbank: https://www.ncbi.nlm.nih.gov/assembly/GCF_000186725.1/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=650377987)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NC_017266.1?report=fasta



**Yersinia pestis strain Cadman plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NZ_CP016275.1

- Info: Isolated in US 1965 from a boy's CSF

- Sources that provide strain information: (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5084870/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2718217770)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP016275.1?report=fasta



**Yersinia pestis strain S19960127 plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NZ_CP045637.1

- Info: Isolated from a pneumonic patient’s organs at necropsy during a plague outbreak that occurred in 1996 in Qayü village, Lhünze County, Shannan Prefecture in Tibet, China

- Sources that provide strain information: (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8096067/)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP045637.1?report=fasta



**Yersinia pestis CA88-4125 plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NC_009596.1

- Info: Isolated from human 1988 California

- Sources that provide strain information: (refseq=genbank: https://www.ncbi.nlm.nih.gov/assembly/GCF_000181455.2/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=640963031)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NC_009596.1?report=fasta



**Yersinia pestis CO92 plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NC_003134.1

- Info: Isolated from a US human who got it from a cat at or before 2001

- Sources that provide strain information: (https://core.ac.uk/reader/13110408?utm_source=linkout, genbank=refseq: https://www.ncbi.nlm.nih.gov/assembly/GCF_000009065.1/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=637000351)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NC_003134.1?report=fasta




**Yersinia pestis D106004 plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NC_017155.1

- Info: Isolated from Apodemus chevrieri in Yulong county in 2006

- Sources that provide strain information: (https://www.ajtmh.org/view/journals/tpmd/81/4/article-p714.xml, refseq=genbank: https://www.ncbi.nlm.nih.gov/assembly/GCF_000022805.1/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=646862350)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NC_017155.1?report=fasta



**Yersinia pestis Z176003 plasmid pMT1, complete sequence**

- NCBI Reference Sequence: NC_014022.1

- Info: Isolated from dead marmot in Tibet Autonomous Region, China

- Sources that provide strain information: (refseq=genbank: https://www.ncbi.nlm.nih.gov/assembly/GCF_000022845.1/, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=646564590)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/NC_014022.1?report=fasta



**Yersinia pestis Pestoides G plasmid pMT sequence**

- GenBank: CP010248.1

- Info: Isolated from Microtus montanus in Republic of Georgia

- Sources that provide strain information: (https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2648501335)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/CP010248.1?report=fasta



**Yersinia pestis subsp. pestis bv. Medievalis strain SCPM-O-B-6530 plasmid pMT, complete sequence**

- GenBank: CP045159.1

- Info: Isolated from Citellophilus tesquorum fleas in Central-Caucasian high-mountain, Russia 2000

- Sources that provide strain information: (https://journals.asm.org/doi/10.1128/mra.01115-21, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2891959310)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/CP045159.1?report=fasta



**Yersinia pestis Pestoides F plasmid pMT, complete sequence**

- GenBank: CP009714.1

- Info: Isolated from host in former Soviet Union

- Sources that provide strain information: (https://link.springer.com/chapter/10.1007/978-0-387-72124-8_2, https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=2630968993)

- Fasta source: https://www.ncbi.nlm.nih.gov/nuccore/CP009714.1?report=fasta



## Step 2: Quality control

**Completed - Quality control the data**
- Data received in FASTA form so already quality check was already completed and confirmed by NCBI.



## Step 3: Alignment

**Completed - using ClustalW**

- Description: a software that uses progressive alignment methods based on pairwise alignment

- Strengths: being commonly used to align sequences, handling large datasets, and being macOS friendly

- Weaknesses and Limitations: not performing well on sequences with high levels of divergence, producing biased results due to its progressive alignment strategy, takes very long to align a large dataset

- Assumptions: assumes all sequences are evolutionarily related and the rate of change in the sequences is constant

_Code:_

```
#Downloaded clustalw version 2.1, identified as "clustalw-2.1-macosx.dmg" from http://www.clustal.org/download/current/, into the software folder. Entered the "data" folder via terminal. Created a fasta file containing all sequences called "pestis-pMT1-all.fasta" to be the input file. Ran alignment with the following commands to get an output file called "pestis-pMT1-all-clustalw.fasta".

clustalw2 -ALIGN -INFILE=pestis-pMT1-all.fasta -OUTFILE=pestis-pMT1-all-clustalw.fasta -OUTPUT=FASTA

```


**Completed - using Mafft**

- Description: a multiple sequence alignment tool that progressively aligns sequences quickly and accurately based on the Fast Fourier transform strategy; uses several different methods for aligning sequences, such as progressive methods, iterative refinement, and consistency-based approaches

- Strengths: handling large datasets, aligning quickly, and being macOS friendly

- Weaknesses and Limitations: not performing well when sequences have a high degree of divergence

- Assumptions: the input sequences are homologous and all lineages have the same evolutionary rate

_Code:_

```
#Downloaded mafft version 7, using the steps from https://anaconda.org/bioconda/mafft, into the software folder. Entered the "data" folder via terminal. Used previously created input fasta file containing all sequences called "pestis-pMT1-all.fasta". Ran alignment with the following commands to get an output file called "pestis-pMT1-all-mafft.fasta".

mafft pestis-pMT1-all.fasta > pestis-pMT1-all-mafft.fasta

```

## Step 4: TrimAl

**Completed - using Trimal**

- Description: a software that trims and filters multiple sequence alignments, removing poorly aligned positions using gap filtering, gap thresholding, and column removal

- Strengths: users can select how sequences are trimmed to remove poorly aligned regions of multiple sequence alignments, which can improve the accuracy of phylogenetic inference

- Weaknesses and Limitations: potential loss of important data if poorly aligned positions are wrongly removed, especially in cases that a user is knowledgeable in setting the filtering parameter

- Assumptions: assumes poorly aligned positions don’t provide much information thus should be removed from the final alignment and that the user has already aligned the sequences before input

_Code:_

```
#mafft-aligned Trimal trimming

#Downloaded Trimal version 1.2, from http://trimal.cgenomics.org/downloads, into the "data" folder. Entered the "data" folder via terminal. Ran Trimal using the previously created input file "pestis-pMT1-all-mafft.fasta" to get an output file called "pestis-pMT1-all-mafft-trimal.fasta".

trimal -in pestis-pMT1-all-mafft.fasta -out pestis-pMT1-all-mafft-trimal.fasta

```

```
#clustalw-aligned Trimal trimming

#Trimal was already downloaded in the "data" folder. Entered the "data" folder via terminal. Ran Trimal using the previously created input file "pestis-pMT1-all-clustalw.fasta" to get an output file called "pestis-pMT1-all-clustalw-trimal.fasta".

trimal -in pestis-pMT1-all-clustalw.fasta -out pestis-pMT1-all-clustalw-trimal.fasta

```


## Step 5: Distance-based tree and parsimony-based tree using the ape and phangorn R packages

- Description: Ape and Phangorn packages in the R software allow phylogenetic tree visualization under the use of various phylogenetic methods 

- Strengths: can use various algorithms for phylogenetic reconstruction, such as maximum likelihood, parsimony, and Bayesian inference; can handle both nucleic acid and protein sequences; handle large data sets and using parallel computing, macOS friendly

- Weaknesses and Limitations: the accuracy of the tree can be affected by the quality of the sequence alignment, choice of evolutionary model, and choice of method of phylogenetic reconstruction

- Assumptions: assume your input is correct; the creation of parsimony-based and distance-based trees assume that sequences evolved independently and that mutations occurred at a constant rate over time

- Note of how to root trees: root all trees using root(tre.ini,outgroup=“name of clade plasmids”) or root(tre.ini, node=number)


**Completed- Parsimony Tree per each Clustalw/Trimal and Mafft/Trimal Alignment with selecting CP000309.1 as root**

_Code:_

```

#mafft-aligned parsimony tree:

#Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-mafft-trimal-parsimony.pdf"

library(ape)

library(adegenet)

library(phangorn)

dna <- fasta2DNAbin(file="/Users/rajvipatel/phylogenetics-final/data/pestis-pMT1-all-mafft-trimal.fasta")

dna2 <- as.phyDat(dna)

tre.ini <- nj(dist.dna(dna,model="raw"))

parsimony(tre.ini, dna2)

tre.pars <- optim.parsimony(tre.ini, dna2)

plot(tre.ini)

nodelabels()

tre2 = tre2 <- root(tre.ini, outgroup = "CP000309.1 263294 bp")

plot(tre2, cex = 0.6)

title("pMT mafft-aligned parsimony-based tree")

```

```

#clustalw-aligned parsimony tree:

#Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-clustalw-trimal-parsimony.pdf"

library(ape)

library(adegenet)

library(phangorn)

dna <- fasta2DNAbin(file="/Users/rajvipatel/phylogenetics-final/data/pestis-pMT1-all-clustalw-trimal.fasta")

dna2 <- as.phyDat(dna)

tre.ini <- nj(dist.dna(dna,model="raw"))

parsimony(tre.ini, dna2)

tre.pars <- optim.parsimony(tre.ini, dna2)

plot(tre.ini)

nodelabels()

tre2 = root(tre.ini,outgroup="CP000309.1 210104 bp")

plot(tre2, cex = 0.6)

title("pMT clustalw-aligned parsimony-based tree")

```

**Completed- Distance Tree per each Clustalw/Trimal and Mafft/Trimal Alignment with selecting CP000309.1 as root**

_Code:_

```

#mafft-aligned distance tree:

#Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-mafft-trimal-distance.pdf"

library(ape)

library(adegenet)

library(phangorn)

dna <- fasta2DNAbin(file="/Users/rajvipatel/phylogenetics-final/data/pestis-pMT1-all-mafft-trimal.fasta")

D <- dist.dna(dna, model="TN93")

tre <- nj(D)

plot(tre)

nodelabels()

tre2 = root(tre,outgroup="CP000309.1 263294 bp")

tre <- ladderize(tre)

plot(tre, cex=.6)
title("pMT mafft-aligned distance-based tree")

```

```

#clustalw-aligned distance tree:

#Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-clustalw-trimal-distance.pdf"

library(ape)

library(adegenet)

library(phangorn)

dna <- fasta2DNAbin(file="/Users/rajvipatel/phylogenetics-final/data/pestis-pMT1-all-clustalw-trimal.fasta")

D <- dist.dna(dna, model="TN93")

tre <- nj(D)

plot(tre)

nodelabels()

tre2 = root(tre,outgroup="CP000309.1 210104 bp")

tre <- ladderize(tre)

plot(tre, cex=.6)
title("pMT clustalw-aligned distance-based tree")

```


## Step 6: Maximum Likelihood Inference Method Per each Clustalw and Mafft-Aligned Tree

**Completed- RAxML-NG Maximum Likelihood Tree per each Clustalw/Trimal and Mafft/Trimal Alignment with selecting CP000309.1 as root**

- Description: a maximum likelihood software that combines rapid bootstrapping and search algorithms to explore the tree space to find the best tree topology

- Strengths: it’s fast and consistent, including with large datasets; can support a wide range of substitution models and can handle multiple data inputs

- Weaknesses and Limitations: it may not always find the optimal tree topology, no information on convergence

- Assumptions: assumes the input sequences are homologous and evolved under a nucleotide substitution model with a constant rate of evolution, assumes the model choices are correct 

_Code:_

```
#Example:

Rajvis-Air:raxml-ng_v1.1.0_macos_x86_64 rajvipatel$ ./raxml-ng --check --msa ng-tutorial/bad.fa --model GTR+G
Rajvis-Air:raxml-ng_v1.1.0_macos_x86_64 rajvipatel$ ./raxml-ng --check --msa ng-tutorial/bad.fa.raxml.reduced.phy --model GTR+G
Rajvis-Air:raxml-ng_v1.1.0_macos_x86_64 rajvipatel$ ./raxml-ng --parse --msa ng-tutorial/prim.phy --model GTR+G
Rajvis-Air:raxml-ng_v1.1.0_macos_x86_64 rajvipatel$ ./raxml-ng --msa ng-tutorial/prim.phy --model GTR+G --prefix T3 --threads 2 --seed 2

```

```

#mafft-aligned raxml-ng tree:

#Downloaded raxml-ng version 1.1, from https://github.com/amkozlov/raxml-ng, into the "data" folder. Using macos "Finder" app, entered the "phylogenetics-final folder" then went to the "data" folder. Duplicated the "pestis-pMT1-all-mafft-trimal.fasta" file and renamed it "pestis-pMT1-all-mafft-trimal-raxml.fasta". Dragged the "pestis-pMT1-all-mafft-trimal-raxml.fasta" file into the "raxml-ng_v1.1.0_macos_x86_64" folder. Entered the "raxml-ng_v1.1.0_macos_x86_64" folder on terminal. Ran the following commands on the "pestis-pMT1-all-mafft-trimal-raxml.fasta" input file.

./raxml-ng --check --msa pestis-pMT1-all-mafft-trimal-raxml.fasta --model GTR+G

./raxml-ng --check --msa pestis-pMT1-all-mafft-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G

./raxml-ng --parse --msa pestis-pMT1-all-mafft-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G

./raxml-ng --msa pestis-pMT1-all-mafft-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G --prefix T5 --threads 2 --search autoMRE=0.01 (without seed so it runs until likelihood score between 2 iterations is less than 1% of current likelihood score)

#Using macos "Finder" app, found the "phylogenetics-folder" with the objective of entering the "raxml-ng_v1.1.0_macos_x86_64" folder. Inside of the "raxml-ng_v1.1.0_macos_x86_64" folder, opened the "T5.raxml.bestTree" file. Copied the Newick formatted data to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-mafft-trimal-raxml.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='((((CP010248.1_263294_bp:0.000043,CP009714.1_263294_bp:0.000030):0.000124,((((((AE017045.1_263294_bp:0.000135,NC_017266.1_263294_bp:0.000168):0.003031,(CP045159.1_263294_bp:0.000142,NC_003134.1_263294_bp:0.000055):0.014568):0.195463,NC_017155.1_263294_bp:0.000067):0.000012,CP000306.1_263294_bp:0.000068):0.000017,NC_014022.1_263294_bp:0.000001):0.000001,CP000309.1_263294_bp:0.000054):0.000107):0.324203,NZ_CP045637.1_263294_bp:0.006637):0.000001,NC_009596.1_263294_bp:0.000001,NZ_CP016275.1_263294_bp:0.000001);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1_263294_bp")

plot(tre2, cex = 0.6)
title("pMT mafft-aligned raxml maximum likelihood-based tree")

```

```

#clustalw-aligned raxml-ng tree:

#Raxml-ng was already downloaded into the "data" folder. Using macos "Finder" app, entered the "phylogenetics-final folder" then went to the "data" folder. Duplicated the "pestis-pMT1-all-clustalw-trimal.fasta" file and renamed it "pestis-pMT1-all-clustalw-trimal-raxml.fasta". Dragged the "pestis-pMT1-all-clustalw-trimal-raxml.fasta" file into the "raxml-ng_v1.1.0_macos_x86_64" folder. Entered the "raxml-ng_v1.1.0_macos_x86_64" folder on terminal. Ran the following commands on the "pestis-pMT1-all-clustalw-trimal-raxml.fasta" input file.

./raxml-ng --check --msa pestis-pMT1-all-clustalw-trimal-raxml.fasta --model GTR+G

./raxml-ng --check --msa pestis-pMT1-all-clustalw-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G

./raxml-ng --parse --msa pestis-pMT1-all-clustalw-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G

./raxml-ng --msa pestis-pMT1-all-clustalw-trimal-raxml.fasta.raxml.reduced.phy --model GTR+G --prefix T6 --threads 2 --search autoMRE=0.01 (without seed so it runs until likelihood score between 2 iterations is less than 1% of current likelihood score)

#Using macos "Finder" app, found the "phylogenetics-folder" with the objective of entering the "raxml-ng_v1.1.0_macos_x86_64" folder. Inside of the "raxml-ng_v1.1.0_macos_x86_64" folder, opened the "T6.raxml.bestTree" file. Copied the Newick formatted data to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-clustalw-trimal-raxml.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(((NC_009596.1_210104_bp:0.000001,NZ_CP016275.1_210104_bp:0.000001):0.000001,NZ_CP045637.1_210104_bp:0.018248):100.000000,(CP000306.1_210104_bp:0.000085,NC_017155.1_210104_bp:0.000095):0.000008,((((CP009714.1_210104_bp:0.000033,CP010248.1_210104_bp:0.000034):0.221188,(((CP045159.1_210104_bp:0.000399,NC_003134.1_210104_bp:0.000343):0.832370,NC_017266.1_210104_bp:0.000153):0.000001,AE017045.1_210104_bp:0.000126):0.198481):0.215744,CP000309.1_210104_bp:0.000067):0.000001,NC_014022.1_210104_bp:0.000001):0.000006);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1_210104_bp")

plot(tre2, cex = 0.6)
title("pMT clustalw-aligned raxml maximum likelihood-based tree")

```

**Completed- IQ-Tree 2 Maximum Likelihood Tree per each Clustalw/Trimal and Mafft/Trimal Alignment with selecting CP000309.1 as root**

- Description: a software that uses fast maximum likelihood calculation algorithms and efficient search strategies to find the best tree topology

- Strengths: fast and efficient, including with large datasets, shown to have the best tree inference accuracy, supports a wide range of substitution models and can handle different types of input data

- Weaknesses and Limitations: it may not always find the optimal tree topology, only information on convergence is amount of time it took

- Assumptions: assumes the model choices are correct, the input sequences are homologous and evolved under a nucleotide substitution model with a constant rate of evolution

_Code:_

```

#Example:

Rajvis-MacBook-Air:iqtree-2.2.0-MacOSX rajvipatel$ bin/iqtree2 -s input.fasta


```

```

#mafft-aligned iqtree2 tree:

#Downloaded iqtree2 version 2.2.0, from http://www.iqtree.org/#download, into the "data" folder. Using macos "Finder" app, entered the "phylogenetics-final folder" then went to the "data" folder. Duplicated the "pestis-pMT1-all-mafft-trimal.fasta" file and renamed it "pestis-pMT1-all-mafft-trimal-iqtree.fasta". Dragged the "pestis-pMT1-all-mafft-trimal-iqtree.fasta" file into the "iqtree-2.2.0-MacOSX" folder. Entered the "iqtree-2.2.0-MacOSX" folder on terminal. Ran the following commands on the "pestis-pMT1-all-mafft-trimal-iqtree.fasta" input file.

bin/iqtree2 -s pestis-pMT1-all-mafft-trimal-iqtree.fasta

#Using macos "Finder" app, found the "phylogenetics-folder" with the objective of entering the "iqtree-2.2.0-MacOSX" folder. Inside of the "iqtree-2.2.0-MacOSX" folder, opened the pestis-pMT1-all-mafft-trimal-iqtree.fasta.treefile file. Copied the Newick formatted data to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-mafft-trimal-iqtree.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(CP000309.1:0.0000541585,((CP000306.1:0.0000676246,(((AE017045.1:0.0001346761,NC_017266.1:0.0001681082):0.0030345631,(NC_003134.1:0.0000547773,CP045159.1:0.0001426097):0.0145850616):0.1952032729,NC_017155.1:0.0000670555):0.0000121843):0.0000172803,NC_014022.1:0.0000003798):0.0000003798,(((NZ_CP016275.1:0.0000007861,NC_009596.1:0.0000003798):0.0000008566,NZ_CP045637.1:0.0066396013):0.3234548284,(CP010248.1:0.0000429104,CP009714.1:0.0000296895):0.0001242404):0.0001075467);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1")

plot(tre2, cex = 0.6)
title("pMT mafft-aligned iqtree maximum likelihood-based tree")

```

```

#clustalw-aligned iqtree2 tree:

#Iqtree2 was already downloaded into the "data" folder. Using macos "Finder" app, entered the "phylogenetics-final folder" then went to the "data" folder. Duplicated the "pestis-pMT1-all-clustalw-trimal.fasta" file and renamed it "pestis-pMT1-all-clustalw-trimal-iqtree.fasta". Dragged the "pestis-pMT1-all-clustalw-trimal-iqtree.fasta" file into the "iqtree-2.2.0-MacOSX" folder. Entered the "iqtree-2.2.0-MacOSX" folder on terminal. Ran the following commands on the "pestis-pMT1-all-clustalw-trimal-iqtree.fasta" input file.

bin/iqtree2 -s pestis-pMT1-all-clustalw-trimal-iqtree.fasta

#Using macos "Finder" app, found the "phylogenetics-folder" with the objective of entering the "iqtree-2.2.0-MacOSX" folder. Inside of the "iqtree-2.2.0-MacOSX" folder, opened the pestis-pMT1-all-clustalw-trimal-iqtree.fasta.treefile file. Copied the Newick formatted data to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-clustalw-trimal-iqtree.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(AE017045.1:0.0001455069,NC_017266.1:0.0001286336,((NC_003134.1:0.0004308973,CP045159.1:0.0005177772):1.2539311755,((((CP000309.1:0.0000669229,NC_014022.1:0.0000004754):0.0000010413,(CP000306.1:0.0000830916,NC_017155.1:0.0000933742):0.0000128895):0.0312880351,((NZ_CP016275.1:0.0000009620,NC_009596.1:0.0000004754):0.0000009591,NZ_CP045637.1:0.0189830752):6.4242885658):0.2307907686,(CP010248.1:0.0000203969,CP009714.1:0.0000452803):0.2998268146):0.2222226247):0.0292935461);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1")

plot(tre2, cex = 0.6)
title("pMT clustalw-aligned iqtree maximum likelihood-based tree")

```

## Step 7: Bayesian Method Per each Parsimony and Distance Tree

**Completed- Mr.Bayes per each Clustalw/Trimal and Mafft/Trimal Alignment**

- Description: a statistical tool to infer the parameters of a probability distribution. It involves creating a prior parameter distribution based on existing knowledge, updating the prior distribution with observed data using Bayes’ theorem, and sampling the parameter’s posterior distribution using the Markov Chain Monte Carlo technique

- Strengths: is macos friendly and it’s strength relies on that it allows the incorporation of prior knowledge into the analysis to make predictions

- Weaknesses and Limitations: can be sensitive to the specification of the prior distribution, resulting in biased outcomes; challenging to interpret results

- Assumptions: assumes the data is independent and identically distributed, assumes the data is in correct Nexus format

_Code:_

```
#Downloaded MrBayes version 3.2.7a, from http://nbisweden.github.io/MrBayes/download.html, into the "data" folder. Using terminal, entered the "data" folder. Created a text file called "mbblock.txt" with the following commands pasted inside:
begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=10000000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
 outgroup CP000309.1;
 mcmc;
 sumt;
end;

```

```

#mafft-aligned Mr.Bayes tree:

#From the "data" folder, created a nexus file called "pestis-pMT1-all-mafft-trimal.nex" using the data from "pestis-pMT1-all-mafft-trimal.fasta". Ran MrBayes using the following commands to get the output file containing tree info called "pestis-pMT1-all-mafft-trimal-mrbayes.nex.con.tre":

cat pestis-pMT1-all-mafft-trimal.nex mbblock.txt > pestis-pMT1-all-mafft-trimal-mrbayes.nex

mb pestis-pMT1-all-mafft-trimal-mrbayes.nex

# Copied the Newick formatted data from the "pestis-pMT1-all-mafft-trimal-mrbayes.con.tre" file to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-mafft-trimal-mrbayes.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(1[&prob=1.00000000e+00,prob(percent)="100"]:6.109002e-05[&length_mean=6.46978020e-05,length_median=6.10900200e-05,length_95%HPD={1.90869600e-05,1.16839500e-04}],10[&prob=1.00000000e+00,prob(percent)="100"]:7.352911e-06[&length_mean=1.06237295e-05,length_median=7.35291100e-06,length_95%HPD={1.00393800e-08,3.18477700e-05}],(((5[&prob=1.00000000e+00,prob(percent)="100"]:3.043328e-05[&length_mean=4.36593589e-05,length_median=3.04332800e-05,length_95%HPD={1.00286800e-08,1.30036900e-04}],7[&prob=1.00000000e+00,prob(percent)="100"]:7.378469e-06[&length_mean=1.06008654e-05,length_median=7.37846900e-06,length_95%HPD={1.01094200e-08,3.15506800e-05}])[&prob=8.48460202e-01,prob(percent)="85"]:1.263909e-04[&length_mean=1.66326746e-04,length_median=1.26390900e-04,length_95%HPD={1.64039700e-08,4.61297500e-04}],6[&prob=1.00000000e+00,prob(percent)="100"]:6.440966e-03[&length_mean=6.43845949e-03,length_median=6.44096600e-03,length_95%HPD={5.71263400e-03,7.17116200e-03}])[&prob=1.00000000e+00,prob(percent)="100"]:3.248014e-01[&length_mean=3.24799386e-01,length_median=3.24801400e-01,length_95%HPD={3.20330100e-01,3.29297100e-01}],(11[&prob=1.00000000e+00,prob(percent)="100"]:4.831269e-05[&length_mean=5.13740157e-05,length_median=4.83126900e-05,length_95%HPD={1.12190400e-05,9.64576100e-05}],13[&prob=1.00000000e+00,prob(percent)="100"]:3.508671e-05[&length_mean=3.77737753e-05,length_median=3.50867100e-05,length_95%HPD={7.55251100e-06,7.37461200e-05}])[&prob=9.97648003e-01,prob(percent)="100"]:1.328197e-04[&length_mean=1.39114450e-04,length_median=1.32819700e-04,length_95%HPD={3.08060800e-05,2.58791600e-04}])[&prob=9.08572122e-01,prob(percent)="91"]:1.220281e-04[&length_mean=1.30119512e-04,length_median=1.22028100e-04,length_95%HPD={1.69033900e-05,2.57650000e-04}],(2[&prob=1.00000000e+00,prob(percent)="100"]:7.145393e-05[&length_mean=7.45666301e-05,length_median=7.14539300e-05,length_95%HPD={2.20339600e-05,1.33405100e-04}],(((3[&prob=1.00000000e+00,prob(percent)="100"]:1.419433e-04[&length_mean=1.45143483e-04,length_median=1.41943300e-04,length_95%HPD={7.18535500e-05,2.20163600e-04}],4[&prob=1.00000000e+00,prob(percent)="100"]:1.779643e-04[&length_mean=1.84492460e-04,length_median=1.77964300e-04,length_95%HPD={7.58723300e-05,3.04809300e-04}])[&prob=1.00000000e+00,prob(percent)="100"]:3.025431e-03[&length_mean=3.02818047e-03,length_median=3.02543100e-03,length_95%HPD={2.55603300e-03,3.49463300e-03}],(8[&prob=1.00000000e+00,prob(percent)="100"]:9.423502e-05[&length_mean=1.12166111e-04,length_median=9.42350200e-05,length_95%HPD={1.89326600e-06,2.66418400e-04}],12[&prob=1.00000000e+00,prob(percent)="100"]:1.623760e-04[&length_mean=1.74515944e-04,length_median=1.62376000e-04,length_95%HPD={3.98599300e-05,3.32807800e-04}])[&prob=1.00000000e+00,prob(percent)="100"]:1.452346e-02[&length_mean=1.45272369e-02,length_median=1.45234600e-02,length_95%HPD={1.32958700e-02,1.57206700e-02}])[&prob=1.00000000e+00,prob(percent)="100"]:1.955124e-01[&length_mean=1.95513894e-01,length_median=1.95512400e-01,length_95%HPD={1.92734700e-01,1.98367000e-01}],9[&prob=1.00000000e+00,prob(percent)="100"]:7.498609e-05[&length_mean=7.85892821e-05,length_median=7.49860900e-05,length_95%HPD={2.56026000e-05,1.38156400e-04}])[&prob=5.90703212e-01,prob(percent)="59"]:2.808193e-05[&length_mean=3.51049500e-05,length_median=2.80819300e-05,length_95%HPD={1.02345300e-08,9.08773700e-05}])[&prob=8.98378802e-01,prob(percent)="90"]:2.889475e-05[&length_mean=3.45114312e-05,length_median=2.88947500e-05,length_95%HPD={6.24967900e-07,8.25825300e-05}]);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="1")

plot(tre2, cex = 0.6)
title("pMT mafft-aligned mrbayes bayesian-based tree")

```

```
#clustalw-aligned Mr.Bayes tree:

#From the "data" folder, created a nexus file called "pestis-pMT1-all-clustalw-trimal.nex" using the data from "pestis-pMT1-all-clustalw-trimal.fasta". Ran MrBayes using the following commands to get the output file containing tree info called "pestis-pMT1-all-clustalw-trimal-mrbayes.nex.con.tre":

cat pestis-pMT1-all-clustalw-trimal.nex mbblock.txt > pestis-pMT1-all-clustalw-trimal-mrbayes.nex

mb pestis-pMT1-all-clustalw-trimal-mrbayes.nex

# Copied the Newick formatted data from the "pestis-pMT1-all-clustalw-trimal-mrbayes.nex.con.tre" file to use in the following commands on software R. Opened software R and inputted the following commands to obtain a tree, which was saved to the desktop under "pMT-clustalw-trimal-mrbayes.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(1[&prob=1.00000000e+00,prob(percent)="100"]:7.240588e-05[&length_mean=7.66487209e-05,length_median=7.24058800e-05,length_95%HPD={2.01034600e-05,1.42040000e-04}],6[&prob=1.00000000e+00,prob(percent)="100"]:9.036560e-06[&length_mean=1.30786045e-05,length_median=9.03656000e-06,length_95%HPD={1.00044100e-08,3.90810900e-05}],(((((2[&prob=1.00000000e+00,prob(percent)="100"]:1.403295e-04[&length_mean=1.44625044e-04,length_median=1.40329500e-04,length_95%HPD={4.25510100e-05,2.49203800e-04}],3[&prob=1.00000000e+00,prob(percent)="100"]:1.464387e-04[&length_mean=1.51815671e-04,length_median=1.46438700e-04,length_95%HPD={4.34266800e-05,2.63860700e-04}])[&prob=9.99520001e-01,prob(percent)="100"]:8.324977e-03[&length_mean=8.33618472e-03,length_median=8.32497700e-03,length_95%HPD={2.24958500e-03,1.42726900e-02}],(4[&prob=1.00000000e+00,prob(percent)="100"]:3.703724e-04[&length_mean=4.28706631e-04,length_median=3.70372400e-04,length_95%HPD={1.49542900e-08,1.02042200e-03}],5[&prob=1.00000000e+00,prob(percent)="100"]:3.206955e-04[&length_mean=3.57370581e-04,length_median=3.20695500e-04,length_95%HPD={1.74027500e-08,7.98419600e-04}])[&prob=1.00000000e+00,prob(percent)="100"]:6.001931e-01[&length_mean=6.00310241e-01,length_median=6.00193100e-01,length_95%HPD={5.80192800e-01,6.20003900e-01}])[&prob=9.98760002e-01,prob(percent)="100"]:1.694790e-01[&length_mean=1.68271182e-01,length_median=1.69479000e-01,length_95%HPD={1.60306100e-01,1.79079200e-01}],(9[&prob=1.00000000e+00,prob(percent)="100"]:3.745253e-05[&length_mean=4.15261999e-05,length_median=3.74525300e-05,length_95%HPD={1.12355300e-08,9.10737400e-05}],10[&prob=1.00000000e+00,prob(percent)="100"]:4.049384e-05[&length_mean=4.40806918e-05,length_median=4.04938400e-05,length_95%HPD={1.06106300e-08,9.38682400e-05}])[&prob=1.00000000e+00,prob(percent)="100"]:1.956228e-01[&length_mean=1.93193702e-01,length_median=1.95622800e-01,length_95%HPD={1.88482800e-01,2.02761400e-01}])[&prob=9.54146728e-01,prob(percent)="95"]:1.640189e-01[&length_mean=1.54649941e-01,length_median=1.64018900e-01,length_95%HPD={9.02291700e-02,1.91691400e-01}],((11[&prob=1.00000000e+00,prob(percent)="100"]:3.382107e-05[&length_mean=4.89017432e-05,length_median=3.38210700e-05,length_95%HPD={1.08282600e-08,1.46069100e-04}],12[&prob=1.00000000e+00,prob(percent)="100"]:1.513242e-05[&length_mean=2.19476710e-05,length_median=1.51324200e-05,length_95%HPD={1.02157800e-08,6.59221700e-05}])[&prob=9.94124008e-01,prob(percent)="99"]:6.680517e-03[&length_mean=7.35942847e-03,length_median=6.68051700e-03,length_95%HPD={2.55138500e-07,1.57124000e-02}],13[&prob=1.00000000e+00,prob(percent)="100"]:9.977361e-03[&length_mean=9.39816515e-03,length_median=9.97736100e-03,length_95%HPD={3.17944800e-04,1.64647000e-02}])[&prob=1.00000000e+00,prob(percent)="100"]:2.501677e+01[&length_mean=2.50406055e+01,length_median=2.50167700e+01,length_95%HPD={2.23066800e+01,2.77436000e+01}])[&prob=9.91709344e-01,prob(percent)="99"]:2.454702e-02[&length_mean=3.95370679e-02,length_median=2.45470200e-02,length_95%HPD={1.96495600e-07,1.60210100e-01}],(7[&prob=1.00000000e+00,prob(percent)="100"]:9.049509e-05[&length_mean=9.49793043e-05,length_median=9.04950900e-05,length_95%HPD={2.97742700e-05,1.72798800e-04}],8[&prob=1.00000000e+00,prob(percent)="100"]:9.989825e-05[&length_mean=1.04327758e-04,length_median=9.98982500e-05,length_95%HPD={3.75695800e-05,1.77694200e-04}])[&prob=6.78033763e-01,prob(percent)="68"]:2.200754e-05[&length_mean=2.65112957e-05,length_median=2.20075400e-05,length_95%HPD={1.35652500e-08,6.51271300e-05}])[&prob=5.38485949e-01,prob(percent)="54"]:1.827197e-05[&length_mean=2.26384687e-05,length_median=1.82719700e-05,length_95%HPD={1.06513600e-08,5.82389300e-05}]);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="1")

plot(tre2, cex = 0.6)
title("pMT clustalw-aligned mrbayes bayesian-based tree")

```

## Step 8: Coalescent Method on each ClustalW and Mafft-Aligned IQ-Tree

**Completed- Astral per each Clustalw/Trimal/IQ-Tree and Mafft/Trimal/IQ-Tree Alignment with selecting Antiqua as root**

- Note: in future studies, multiple loci of species trees are combined not just one species tree input to make one gene tree output

- Description: a software that uses the coalescent method to reconstruct evolutionary history, involves using posterior probability to determine the probability of a species tree using the gene tree

- Strengths: had the most simple and straightforward instructions for a beginner to run it on mac

- Weaknesses and Limitations: it’s sensitive to missing data and general gene tree estimation error, instructions weren't straighforward enough to realize that the input involves information to multiple loci

- Assumptions: assumes that the gene tree is correct and the data in input in Newick format

_Code:_

```

#mafft-aligned iqtree-based astral coalescent tree:

#Installed Astral version 5.7.8, from https://github.com/smirarab/ASTRAL/blob/master/README.md#installation, into "data" folder. Using terminal, went into the "Astral"" folder. Created a tre file containing mafft-aligned IQ-tree data in Newick format named "pestis-pMT1-all-mafft-trimal-iqtree.tre". Ran the following commands to get an output file called "pestis-pMT1-all-mafft-trimal-iqtree-astral.tre"

java -jar astral.5.7.8.jar -i pestis-pMT1-all-mafft-trimal-iqtree.tre -o pestis-pMT1-all-mafft-trimal-iqtree-astral.tre

#Copied the Newick format data from the "pestis-pMT1-all-mafft-trimal-iqtree-astral.tre" output file and ran the following commands in R to get a tree that was saved in desktop under "pMT-mafft-trimal-iqtree-astral.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(CP000309.1,(((CP010248.1,CP009714.1)0.67:0.2876820724517809,(NZ_CP045637.1,(NZ_CP016275.1,NC_009596.1)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809,(NC_014022.1,(CP000306.1,(NC_017155.1,((AE017045.1,NC_017266.1)0.67:0.2876820724517809,(CP045159.1,NC_003134.1)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809):0.0);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1")

plot(tre2, cex = 0.6)
title("pMT mafft-aligned iqtree-based astral coalescent tree")

```

```
#clustalw-aligned iqtree-based astral coalescent tree:

#Astral was already downloaded into "data" folder. Using terminal, went into the "Astral"" folder. Created a tre file containing clustalw-aligned IQ-tree data in Newick format named "pestis-pMT1-all-clustalw-trimal-iqtree.tre". Ran the following commands to get an output file called "pestis-pMT1-all-clustalw-trimal-iqtree-astral.tre"

java -jar astral.5.7.8.jar -i pestis-pMT1-all-clustalw-trimal-iqtree.tre -o pestis-pMT1-all-clustalw-trimal-iqtree-astral.tre

#Copied the Newick format data from the "pestis-pMT1-all-clustalw-trimal-iqtree-astral.tre" output file and ran the following commands in R to get a tree that was saved in desktop under "pMT-clustalw-trimal-iqtree-astral.pdf"

library(ape)

library(adegenet)

library(phangorn)

mytree <- read.tree(text='(AE017045.1,(NC_017266.1,((NC_003134.1,CP045159.1)0.67:0.2876820724517809,((CP010248.1,CP009714.1)0.67:0.2876820724517809,((NZ_CP045637.1,(NZ_CP016275.1,NC_009596.1)0.67:0.2876820724517809)0.67:0.2876820724517809,((CP000306.1,NC_017155.1)0.67:0.2876820724517809,(NC_014022.1,CP000309.1)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809)0.67:0.2876820724517809):0.0);')

plot(mytree)

nodelabels()

tre2 = root(mytree,outgroup="CP000309.1")

plot(tre2, cex = 0.6)
title("pMT clustalw-aligned iqtree-based astral coalescent tree")

```


## Step 9: Co-estimation Method on each ClustalW and Mafft-Aligned IQ-Tree

**Completed- Beast per each Clustalw/Trimal and Mafft/Trimal Alignment with selecting Antiqua as root**

- Description: allows coestimations to be created using prior information from Bayesian framework, uses Markov chain Monte chain sampling to create trees

- Strengths: easy to install and run on mac

- Weaknesses and Limitations: several programs to download, long instructions, needing multiple runs to ensure convergence, difficult to interpret

- Assumptions: assumes sequences are correct and aligned properly, assumes each branch evolves at the same evolutionary rate

_Code:_

```

# Downloaded Beast version 2.7.4, from https://beast.community/installing, into the "data" folder and followed the first and second tutorial from the following website: https://beast.community/second_tutorial. When creating the .xml file from .nex mrbayes file using software Beauti version 2.7.4, the following parameters were chosen: Jc69 as the site substitution model, strict clock with mean clock rate of 1.0, Yule model of speciation, checked the "Selected sample from prior" box. Ran beast with 10,000,000 iterations 5 times for both the clustalw-aligned and mafft-aligned data.

# Was unable to visualize the MCC tree containing multiple replicates. Nevertheless, I still attempted to visualize the tree by performing multiple attempts of the following steps found from https://beast.community/second_tutorial. Combined all output data for .log files into one using LogCombiner version 2.7.4. Combined all output data for .trees files into one using LogCombiner. The output files are found in "BEAST 2.7.4" file named as follows: "pestis-pMT1-all-clustalw-trimal-beast-sum.log", "pestis-pMT1-all-mafft-trimal-beast-sum.log", "pestis-pMT1-all-clustalw-trimal-beast-sum.trees", "pestis-pMT1-all-mafft-trimal-beast-sum.trees". Tracer was downloaded through http://beast.community/tracer to summarize parameter estimates for "pestis-pMT1-all-clustalw-trimal-beast-sum.log" and "pestis-pMT1-all-mafft-trimal-beast-sum.log". For the next step in building the MCC Tree, TreeAnnotator was used with the following parameters: burn in percentage of 10, target tree type of maximum clade credibility tree, and node height of median height. The input data were "pestis-pMT1-all-clustalw-trimal-beast-sum.trees" and "pestis-pMT1-all-mafft-trimal-beast-sum.trees" to get the output, containing the tree, on "pestis-pMT1-all-clustalw-trimal-beast-sum-tree.mcc.tre” and "pestis-pMT1-all-mafft-trimal-beast-sum-tree.mcc.tre”. To visualize the MCC tree, the .dmg FigTree version 1.4.4 file was downloaded from https://github.com/rambaut/figtree/releases and inputted the .mcc.tre files. Was unable to see anything on the screen.

# Due to the previously discussed issue, instead of replicate combinations being visualized, replicate 5 from mafft and clustalw run data was turned into an MCC tree. Both trees were rooted using strain CP000309.1 (Antiqua). The tree tips were annotated and saved to Desktop under "pMT-mafft-trimal-beast.pdf" and "pMT-clustalw-trimal-beast.pdf".


```

## Step 10: Annotation of all phylogenetic trees
- For each tip, noted strain ID, strain name, and origin of samples. Highlighted the strain according to the historically identified waves of plague: the first wave involves the biovar Antiqua, which was highlighted as orange on the tree; the second wave involves the biovar Medievalis, which was highlighted as red on the tree; and the third wave involves the biovar Orientalis, which was highlighted as gray on the tree.

