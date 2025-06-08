# Run dCNS for Panicoids

#### Get K-mer distribution
`jellyfish  count -t 20 -C -m 20 -s 50G -o REF_20mer path_to_reference_genome.fa`

`jellyfish histo -o REF_20mer.histo REF_20mer`

`jellyfish dump REF_20mer > REF_20mer_dumps.fa`

#### Obtain the longest transcript for each gene model for masking
`python3 dCNS/scripts/longestTranscript.py -g path_to_reference_genome.fa -f path_to_reference_annot.gff3 -t false -o REF_GENE.fa`

`minimap2 -ax splice -a -uf -C 1 -k 12 -P -t 12 --cs path_to_reference_genome.fa REF_GENE.fa > REF.sam`

#### Mask Genome
As suggested by Song et al., we chose -f parameter based on the distribution obtained from k-mer histogram (REF_20mer.histo)

`dCNS/dCNS maskGenome -i path_to_reference_genome.fa -o masked_REF_k20_25.fa -s REF.sam -c REF_GENE.fa -k REF_20mer_dumps.fa -f 25` 
