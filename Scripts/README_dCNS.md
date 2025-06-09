# Run dCNS for Panicoids

#### 1. Get K-mer distribution for all reference and query genomes
```
jellyfish  count -t 20 -C -m 20 -s 50G -o REF_20mer path_to_reference_genome.fa

jellyfish histo -o REF_20mer.histo REF_20mer

jellyfish dump REF_20mer > REF_20mer_dumps.fa
```

#### 2. Obtain the longest transcript for each gene model for masking
```
python3 dCNS/scripts/longestTranscript.py -g path_to_reference_genome.fa -f path_to_reference_annot.gff3 -t false -o REF_GENE.fa

minimap2 -ax splice -a -uf -C 1 -k 12 -P -t 12 --cs path_to_reference_genome.fa REF_GENE.fa > REF.sam
```

#### 3. Mask Genomes
As suggested by Song et al., we chose -f parameter based on the distribution obtained from k-mer histogram (REF_20mer.histo)
K=25,21,30,19 for PhFIL, PhHAL, SV, and PV

```
dCNS/dCNS maskGenome -i path_to_reference_genome.fa -o masked_REF_k20_25.fa -s REF.sam -c REF_GENE.fa -k REF_20mer_dumps.fa -f 25
```

#### 4. Extract sequences
```
python3 dCNS/scripts/extractInterGeneticSequence/sequenceUpStreamGeneAndDownStreamV2.py \
-g REF.gff3 -r path_to_reference_genome.fa -c REF_GENE.fa -q path_to_query_genome.fa \
-s PHFIL.sam -o dCNS_QNAME_REF
```
#### 5. Write the commands to run parallelly
```
ls | awk '{print("dCNS/dCNS cut1Gap -ra masked_REF_k20_xx_cds.fa -qa masked_Query_k20_xx.fa -i "$1" -r reference -o "$1".5")}' > command1
```
#### 6. Concat and clean the output
```
perl dCNS/scripts/combineCnsSamFiles.pl dCNS_QNAME_REF > 5.sam

cat 5.sam| sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}'  | samtools view -O BAM --reference masked_REF_k20_xx_cds.fa | samtools sort > 5.bam

seqkit locate -F --only-positive-strand --bed -m 0 -p n masked_REF_k20_xx_cds.fa > ns.bed

bedtools merge -i ns.bed > ns_megered.bed

bamToBed -i 5.bam | bedtools sort -i | bedtools merge > 5.bed

bedtools subtract -a 5.bed -b ns_megered.bed > 5_nons.bed
```
