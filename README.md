# JumpORF
This package is developed for identifying additional open reading frames (aORF) of complete LTR retrotransposons (LTR-RTE), DNA transposons (DNA-TE), LINE, DIRS and Penelope-like elements within plant genomes. 

# Workflow
Complete transposable elements (TE) are identified using domain-based approach. For LTR-RTE detection protein domains identified by DANTE tool searched for structural features (such as long terminal repeats, PBS and PPT) by DANTE_LTR tool. For other TEs protein domains are identified using TEsorter. At the next step nucleotide sequence of the most conserved protein domain within each TE type are clustered to get individual TE family. Next, the 5' and 3' regions of the candidate domains within families are extended for dynamic boundary detection using BLASTn. If conserved part beyond domain edges is found this region is marked as TE element candidate and searched for structural features such as terminal inverted repeats for DNA-TEs (subclass 1) or polyA tails for LINE elements. Finally, for gene prediction within TEs, Augustus tool is used, gene prediction, contained basal protein domains are filtered out to obtain additional ORF set. 

# Installation with Docker
```
git clone https://github.com/soyboy-hub/JumpORF
cd JumpORF
docker compose build
```

# Usage
```
docker run -v $(pwd):/home jump_orf JumpORF genome.fna outpath out_prefix threads gene_model 0.8 5 0.5 no
```
# Testing

Before using JumpORF you should test the installation with a sample genome. If your test finished, installation is complete. If the test fails, feel free to open the new issue. 
```
cd JumpORF/testing
#get TAIR10.1 Arabidopsis thaliana genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz
docker run -v $(pwd):/home jump_orf JumpORF GCF_000001735.4_TAIR10.1_genomic.fna ./out_testing arabidopsis arabidopsis 50 arabidopsis 0.8 5 0.5 no
```
