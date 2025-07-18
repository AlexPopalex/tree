#Please, note that the genomes of B. r. redtenbacheri and B. g. grandii available at ncbi under PRJNA962493 and PRJNA1251886 underwent the ncbi quality check upon upload. This lead to removal of 2 (Brsri) and 3 (Bgigi) short sequences from the scaffolds, respectively.  

#adapter autodetection and removal and quality trimming qual < 20 from 5' end (trimgalore/0.6.6, cutadapt/2.10)
../../../Software/TrimGalore-master/trim_galore -j 20 --paired Bas_R1.fq.gz Bas_R2.fq.gz

#mapping reads Bas and Bgibi to Bgigi ref genome (mod) for pseudoref construction (bwa/0.7.17, samtools/1.15.1)
sed 's/scaffold_/Bgigi_v5_scf/g' Bgigi_v5.fasta > Bgigi_v5_mod.fa

bwa index Bgigi_v5_mod.fa
bwa mem -t 40 -R "@RG\tID:Bas\tSM:Bas\tLB:Bas\tPL:Illumina" \
Bgigi_v5_mod.fa \
../reads/Bas_R1_val_1.fq.gz \
../reads/Bas_R2_val_2.fq.gz > Bas_Bgigi.sam
samtools view -@ 40 -b -h Bas_Bgigi.sam > Bas_Bgigi.bam
rm -f Bas_Bgigi.sam
samtools sort -@ 40 Bas_Bgigi.bam -o Bas_Bgigi_sorted.bam
rm -f Bas_Bgigi.bam

#pseudoref construction with IUPAC chars for het. genotypes, leaving ref base for missing genotypes  (bcftools/1.15.1)
bcftools mpileup --threads 40 -Ou --annotate FORMAT/DP -f Bgigi_v5_mod.fa Bas_Bgigi_sorted.bam | bcftools call --threads 40 -m -Oz -o Bas_Bgigi.vcf.gz
bcftools index Bas_Bgigi.vcf.gz

bcftools consensus -H I -f Bgigi_v5_mod.fa Bas_Bgigi.vcf.gz > Bas_Bgigi_HI.fa

#competitive mapping of long reads (samtools/1.15.1, minimap2/2.14)
sed 's/Bgigi_v5_scf/Bas_cns_scf/g' Bas_Bgigi_HI.fa > Bas_cns_HI.fa
cat Brsri_v3.fasta Bgigi_v5_mod.fa Bas_cns_HI.fa > Brsri_Bgigi_Bas_HI_cmpmp.fa

../../../Software/minimap2-2.26_x64-linux/minimap2 -I 400g -x map-hifi -N 50 -a -t 40 -R "@RG\tID:Blm\tSM:Blm" Brsri_Bgigi_Bas_HI_cmpmp.fa ../reads/Blm.fq.gz > Blm_Brsri_Bgigi_Bas_HI_cmpmp.sam
samtools view -@ 40 -h -S Blm_Brsri_Bgigi_Bas_HI_cmpmp.sam -o Blm_Brsri_Bgigi_Bas_HI_cmpmp.bam
samtools sort -@ 40 Blm_Brsri_Bgigi_Bas_HI_cmpmp.bam -o Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted.bam
samtools index -@ 40 Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted.bam 
rm -f Blm_Brsri_Bgigi_Bas_HI_cmpmp.sam
rm -f Blm_Brsri_Bgigi_Bas_HI_cmpmp.bam


#filter for unambiguously mapping reads i.e. reads that have a primary alignment to one of the genomes and no secondary or supplementary alignments to the other genomes; secondary or supplementary alignments within one genome are allowed though (samtools/1.15.1, seqtk/1.4-r122), needs to be adapted for B. whitei accordingly.
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted.bam -h -F 2308 -@ 20 -o Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_prim.bam
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted.bam -h --rf 4 --rf 256 --rf 2048 -@ 20 -o Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_sec.bam
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_prim.bam | awk '{if ($3 ~/Brsri/) print $1}' > Blm_rsri_prim.txt
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_prim.bam | awk '{if ($3 ~/Bgigi/) print $1}' > Blm_gigi_prim.txt
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_prim.bam | awk '{if ($3 ~/Bas/) print $1}' > Blm_as_prim.txt
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_sec.bam | awk '{if ($3 ~/Brsri/) print $1}' > Blm_rsri_sec.txt
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_sec.bam | awk '{if ($3 ~/Bgigi/) print $1}' > Blm_gigi_sec.txt
samtools view Blm_Brsri_Bgigi_Bas_HI_cmpmp_sorted_sec.bam | awk '{if ($3 ~/Bas/) print $1}' > Blm_as_sec.txt
cat Blm_gigi_sec.txt Blm_as_sec.txt > Blm_gigias_sec.txt
cat Blm_rsri_sec.txt Blm_as_sec.txt > Blm_rsrias_sec.txt
cat Blm_rsri_sec.txt Blm_gigi_sec.txt > Blm_rsrigigi_sec.txt
grep -v -f Blm_gigias_sec.txt Blm_rsri_prim.txt > Blm_rsri_una.txt
grep -v -f Blm_rsrias_sec.txt Blm_gigi_prim.txt > Blm_gigi_una.txt
grep -v -f Blm_rsrigigi_sec.txt Blm_as_prim.txt > Blm_as_una.txt

seqtk subseq ../reads/Blm.fq.gz Blm_as_una.txt | bgzip > Blm_as.fq.gz
seqtk subseq ../reads/Blm.fq.gz Blm_gigi_una.txt | bgzip > Blm_gigi.fq.gz
seqtk subseq ../reads/Blm.fq.gz Blm_rsri_una.txt | bgzip > Blm_rsri.fq.gz


# competitive mapping for paired short reads yielding unambiguously mapping reads (bbmap/38.63)
bbsplit.sh threads=40 ref=Brsri_v3.fasta,Bgigi_v5_mod.fa in=../reads/SEhaplome_R1_val_1.fq.gz in2=../reads/SE_haplome_R2_val_2.fq.gz basename=../reads/SEhaplome-%-bbsplit_R#.fq.gz refstats=SEhaplome-refstats scafstats=SEhaplome-scafstats ambiguous2=toss


#mapping (phased) short reads to maternal parent (Brsri_v3.fasta) for reference (bwa/0.7.17, samtools/1.15.1)
bwa index Brsri_v3.fasta
bwa mem -t 40 -R "@RG\tID:SEhaplome\tSM:SEhaplome\tLB:SEhaplome\tPL:Illumina" \
Brsri_v3.fasta \
SEhaplome-Brsri_v3-bbsplit_R1.fq.gz \
SEhaplome-Brsri_v3-bbsplit_R2.fq.gz > SEhaplome_Brsri.sam
samtools view -@ 40 -b -h SEhaplome_Brsri.sam > SEhaplome_Brsri.bam
rm -f SEhaplome_Brsri.sam
samtools sort -@ 40 SEhaplome_Brsri.bam -o SEhaplome_Brsri_sorted.bam
rm -f SEhaplome_Brsri.bam


#mapping (phased) long reads to maternal parent (Brsri_v3.fasta) for reference (samtools/1.15.1, minimap2/2.14)
../../../Software/minimap2-2.26_x64-linux/minimap2 -I 400g -x map-hifi -N 50 -a -t 40 -R "@RG\tID:Brsri\tSM:Brsri" Brsri_v3.fasta ../reads/Brsri.fq.gz > Brsri_Brsri.sam
samtools view -@ 40 -h -S Brsri_Brsri.sam -o Brsri_Brsri.bam
samtools sort -@ 40 Brsri_Brsri.bam -o Brsri_Brsri_sorted.bam
samtools index -@ 40 Brsri_Brsri_sorted.bam 
rm -f Brsri_Brsri.sam
rm -f Brsri_Brsri.bam


#filter alignments for uniquely mapping reads (paired end including the option -f 2 for retaining properly paired reads only; samtools/1.15.1)
samtools view Bas_Brsri_sorted.bam -h -f 2 -F 2308 -q 20 -@ 40 -o Bas_Brsri_sorted_uniq.bam
samtools index Bas_Brsri_sorted_uniq.bam -@ 40
samtools view Bgigi_Brsri_sorted.bam -h -F 2308 -q 20 -@ 40 -o Bgigi_Brsri_sorted_uniq.bam
samtools index Bgigi_Brsri_sorted_uniq.bam -@ 40


#variant calling for pseudoreference construction from (phased) reads for subsequent tree reconstruction in gvcf format and per scaffold (bcftools/1.15.1)
bcftools mpileup --threads 40 -Ou --annotate FORMAT/DP,FORMAT/AD -r Brsri_v3_scf1 -f Brsri_v3.fasta -b bam_list_uniq | bcftools call --threads 40 -g 0 -m -Oz -o scf1.vcf.gz
bcftools index scf1.vcf.gz

bcftools concat --threads 20 -f vcf_list -O z -o genome.vcf.gz


#filter variants (bcftools/1.15.1)
bcftools norm --threads 40 -f Brsri_v3.fasta -d all -O z -o genome_norm.vcf.gz genome.vcf.gz
bcftools view --exclude-types indels --threads 40 -O z genome_norm.vcf.gz > genome_norm_snp.vcf.gz
bcftools view -m 1 -M 2 --threads 40 -O z genome_norm_snp.vcf.gz > genome_norm_snp_an12.vcf.gz
bcftools filter -e 'QUAL<20' --threads 40 -O z genome_norm_snp_an12.vcf.gz > genome_norm_snp_an12_ql20.vcf.gz
bcftools filter -e 'FMT/DP<4' --set-GTs . --threads 40 -O z genome_norm_snp_an12_ql20.vcf.gz > genome_norm_snp_an12_ql20_dp4.vcf.gz
bcftools index genome_norm_snp_an12_ql20_dp4.vcf.gz

# calculate per-genotype sequencing depth and remove sites with depth above threshold, i.e. (median(depth) + 2 * standard deviation(depth)) * #samples (bcftools/1.15.1)
zcat genome_norm_snp_an12_ql20_dp4.vcf.gz | awk '{if (($8 ~ /AN/) && ($8 ~ /DP/)) print $8}' | grep -o '[^;]*DP=[^;]*' | sed 's/DP=//g' > DP
zcat genome_norm_snp_an12_ql20_dp4.vcf.gz | awk '{if (($8 ~ /AN/) && ($8 ~ /DP/)) print $8}' | grep -o '[^;]*AN=[^;]*' | sed 's/AN=//g' > AN
paste DP AN > DPAN
cat DPAN | awk '{if ($2!="0") print $1,$2,$1/$2}' > DPANdiv
bcftools filter -e 'INFO/DP>187' --threads 40 -O z genome_norm_snp_an12_ql20_dp4.vcf.gz > genome_norm_snp_an12_ql20_dp4_dp187.vcf.gz


# continue filtering (bcftools/1.15.1)
bcftools filter -e 'FMT/GT="het"' --set-GTs . --threads 40 -O z genome_norm_snp_an12_ql20_dp4_dp187.vcf.gz > genome_norm_snp_an12_ql20_dp4_dp187_hetgt.vcf.gz
bcftools index genome_norm_snp_an12_ql20_dp4_dp187_hetgt.vcf.gz
bcftools convert --gvcf2vcf --threads 40 -f Brsri_v3.fasta -O z -o genome_norm_snp_an12_ql20_dp4_dp187_hetgt_conv.vcf.gz genome_norm_snp_an12_ql20_dp4_dp187_hetgt.vcf.gz
bcftools index genome_norm_snp_an12_ql20_dp4_dp187_hetgt_conv.vcf.gz


#construct pseudoreferences per sample with Ns for missing genotypes and sites (samtools/1.15.1, bcftools/1.15.1)
for i in $(cat sample_list); do samtools faidx Brsri_v3.fasta -r scf_list | bcftools consensus -M N -a N -H I -s $i -p $i -o pseudorefs/$i.fa genome_norm_snp_an12_ql20_dp4_dp187_hetgt_conv.vcf.gz; done


#convert psequdoreferences into alignments per scaffold
cat *.fa > genomes.fa
cat genomes.fa | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > genomes_ni.fa

grep "scf1\b" -A 1 --no-group-separator genomes_ni.fa > genomes_scf1.fa
sed -i 's/Brsri_v3_scf1//g' genomes_scf1.fa


#remove sites containing a percentage of N's (trimal/1.4)
sed 's/N/\-/g;s/n/\-/g' genomes_scf1.fa > genomes_scf1_gaps.fa
../../../../Software/trimal/source/trimal -gt 1 -in genomes_scf1_gaps.fa -out genomes_scf1_1.fa
../../../../Software/trimal/source/trimal -gt 0.75 -in genomes_scf1_gaps.fa -out genomes_scf1_075.fa
../../../../Software/trimal/source/trimal -gt 0.5 -in genomes_scf1_gaps.fa -out genomes_scf1_05.fa
sed -i 's/\-/N/g' genomes_scf1_075.fa
sed -i 's/\-/N/g' genomes_scf1_05.fa


#calculate tree (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_scf1_1.fa -m MFP -B 1000 --redo


#calculate constrained tree with best fitting model for topology testing (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_075.fasta -g hybrido_h1_topo.treefile -m TVM+F+I --prefix hybrido_h1 --redo

#test best fitting ML tree vs the tree with the constrained topology (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_075.fasta -te ../genomes_075.fasta.treefile -z hybrido.treefiles -zb 10000 -zw --prefix hybridozw --redo

