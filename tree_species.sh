#Please, note that the genomes of B. r. redtenbacheri and B. g. grandii available at ncbi under PRJNA962493 and PRJNA1251886 underwent the ncbi quality check upon upload. This lead to removal of 2 (Brsri) and 3 (Bgigi) short sequences from the scaffolds, respectively.  

#adapter autodetection and removal and quality trimming qual < 20 from 5' end, for unpaired RAD data remove the --paired option (trimgalore/0.6.6, cutadapt/2.10)
../../../Software/TrimGalore-master/trim_galore -j 20 --paired Bas_R1.fq.gz Bas_R2.fq.gz


#mapping RAD reads to maternal parent (Brsri_v3.fasta) and filter (bwa/0.7.17, samtools/1.15.1)
bwa index Brsri_v3.fasta
for i in $(cat ../tree_rad/full_list.txt); do
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
Brsri_v3.fasta \
../tree_rad/${i}-trimmed.fq > ${i}_Brsri.sam
samtools view -@ 40 -b -h ${i}_Brsri.sam > ${i}_Brsri.bam
rm -f ${i}_Brsri.sam
samtools sort -@ 40 ${i}_Brsri.bam -o ${i}_Brsri_sorted.bam
rm -f ${i}_Brsri.bam
samtools view ${i}_Brsri_sorted.bam -h -F 2308 -q 20 -@ 40 -o ${i}_Brsri_sorted_uniq.bam


#mapping WGS short reads to maternal parent (Brsri_v3.fasta) and filter (bwa/0.7.17, samtools/1.15.1)
bwa index Brsri_v3.fasta
bwa mem -t 40 -R "@RG\tID:Bas\tSM:Bas\tLB:Bas\tPL:Illumina" \
Brsri_v3.fasta \
../tree_genomes/reads/Bas_R1_val_1.fq.gz \
../tree_genomes/reads/Bas_R2_val_2.fq.gz > Bas_Brsri.sam
samtools view -@ 40 -b -h Bas_Brsri.sam > Bas_Brsri.bam
rm -f Bas_Brsri.sam
samtools sort -@ 40 Bas_Brsri.bam -o Bas_Brsri_sorted.bam
rm -f Bas_Brsri.bam
samtools view Bas_Brsri_sorted.bam -h -f 2 -F 2308 -q 20 -@ 40 -o Bas_Brsri_sorted_uniq.bam
rm -f Bas_Brsri_sorted.bam


#mapping PacBio long reads to maternal parent (Brsri_v3.fasta) and filter (bwa/0.7.17, samtools/1.15.1)
./../Software/minimap2-2.26_x64-linux/minimap2 -I 400g -x map-hifi -N 50 -a -t 40 -R "@RG\tID:Bgigi\tSM:Bgigi" Brsri_v3.fasta ../tree_genomes/reads/Bgigi.fq.gz > Bgigi_Brsri.sam
samtools view -@ 40 -h -S Bgigi_Brsri.sam -o Bgigi_Brsri.bam
samtools sort -@ 40 Bgigi_Brsri.bam -o Bgigi_Brsri_sorted.bam
samtools index -@ 40 Bgigi_Brsri_sorted.bam 
rm -f Bgigi_Brsri.sam
rm -f Bgigi_Brsri.bam
samtools view Bgigi_Brsri_sorted.bam -h -F 2308 -q 20 -@ 40 -o Bgigi_Brsri_sorted_uniq.bam
rm -f Bgigi_Brsri_sorted.bam


#variant calling for pseudoreference construction from (phased) reads for subsequent tree reconstruction in gvcf format and per scaffold (bcftools/1.15.1)
bcftools mpileup --threads 40 -Ou --annotate FORMAT/DP,FORMAT/AD -r Brsri_v3_scf1 -f Brsri_v3.fasta -b full_list_genomes.txt | bcftools call --threads 40 -v -m -Oz -o scf1.vcf.gz
bcftools index scf1.vcf.gz

bcftools concat --threads 20 -f vcf_list -O z -o genome.vcf.gz


#filter variants (bcftools/1.15.1)
bcftools norm --threads 40 -f Brsri_v3.fasta -d all -O z -o genome_norm.vcf.gz genome.vcf.gz
bcftools view --exclude-types indels --threads 40 -O z genome_norm.vcf.gz > genome_norm_snp.vcf.gz
bcftools view -m 2 -M 2 --threads 40 -O z genome_norm_snp.vcf.gz > genome_norm_snp_an22.vcf.gz
bcftools filter -e 'QUAL<20' --threads 40 -O z genome_norm_snp_an22.vcf.gz > genome_norm_snp_an22_ql20.vcf.gz
bcftools filter -e 'FMT/DP<4' --set-GTs . --threads 40 -O z genome_norm_snp_an22_ql20.vcf.gz > genome_norm_snp_an22_ql20_dp4.vcf.gz
bcftools view -i 'F_MISSING<0.5' --threads 40 -O z genome_norm_snp_an22_ql20_dp4.vcf.gz > genome_norm_snp_an22_ql20_dp4_mm05.vcf.gz
bcftools index genome_norm_snp_an22_ql20_dp4_mm05.vcf.gz

bcftools stats --threads 40 -s- genome_norm_snp_an22_ql20_dp4_mm05.vcf.gz > stats
grep -A 577 -w nMissing stats | awk '{print $0}' | tail -n +2 | awk '{print $3,$14/1240318}' > stats.imiss
cat stats.imiss | awk '{if ($2 > 0.75) print $1}' > rm
bcftools view --samples-file ^rm --threads 40 -O z genome_norm_snp_an22_ql20_dp4_mm05.vcf.gz > genome_norm_snp_an22_ql20_dp4_mm05_rm.vcf.gz


#make pca (plink/1.90)
~/Software/plink/plink --pca --allow-extra-chr --vcf genome_norm_snp_an22_ql20_dp4_mm05.vcf.gz --out genome_norm_snp_an22_ql20_dp4_mm05


#plot (R core team 2023)
plot(pca$V4~pca$V3,col="grey",pch=19,xlab="PC1 (28.7%)",ylab="PC2 (19%)",cex.lab=1.5,main="",axes=F,col.lab="#595959")
axis(1,col="#595959",col.ticks="#595959",col.axis="#595959",cex.axis=1.5)
axis(2,col="#595959",col.ticks="#595959",col.axis="#595959",cex.axis=1.5)
box(col="#595959")
points(pca_atticus$V4~pca_atticus$V3,col="#71bf44",pch=19)
points(pca_gbenazzii$V4~pca_gbenazzii$V3,col="#f7941d",pch=19)
points(pca_ggrandii$V4~pca_ggrandii$V3,col="#be1e2d",pch=19)
points(pca_gmaretimi$V4~pca_gmaretimi$V3,col="#f9ed32",pch=19)
points(pca_hybridoNW$V4~pca_hybridoNW$V3,col="#f47f72",pch=19)
points(pca_hybridoSE$V4~pca_hybridoSE$V3,col="#ec008c",pch=19)
points(pca_whitei$V4~pca_whitei$V3,col="#92278f",pch=19)
points(pca_rredtenbacheri$V4~pca_rredtenbacheri$V3,col="#0599ce",pch=19)
points(pca_unknwon$V4~pca_unknwon$V3,col="#595959",pch=19)
points(pca_lynceorum$V4~pca_lynceorum$V3,col="#8b5e3c",pch=19)
points(pca$V4[pca$V1=="Bas"]~pca$V3[pca$V1=="Bas"],col="grey",pch=19,cex=2)
points(pca$V4[pca$V1=="Bgibi"]~pca$V3[pca$V1=="Bgibi"],col="#f7941d",pch=19,cex=2)
points(pca$V4[pca$V1=="Bgigi"]~pca$V3[pca$V1=="Bgigi"],col="grey",pch=19,cex=2)
points(pca$V4[pca$V1=="Blm"]~pca$V3[pca$V1=="Blm"],col="grey",pch=19,cex=2)
points(pca$V4[pca$V1=="Brsri"]~pca$V3[pca$V1=="Brsri"],col="#0599ce",pch=19,cex=2)
points(pca$V4[pca$V1=="SEhaplome"]~pca$V3[pca$V1=="SEhaplome"],col="#0599ce",pch=19,cex=2)
points(pca$V4[pca$V1=="NWhaplome"]~pca$V3[pca$V1=="NWhaplome"],col="#0599ce",pch=19,cex=2)
points(pca$V4[pca$V1=="Bwi"]~pca$V3[pca$V1=="Bwi"],col="grey",pch=19,cex=2)
dev.copy2eps(file="pc1v2.eps",width=8,height=8)

plot(pca$V5~pca$V3,col="grey",pch=19,xlab="PC1 (28.7%)",ylab="PC3 (13%)",cex.lab=1.5,main="",axes=F,col.lab="#595959")
axis(1,col="#595959",col.ticks="#595959",col.axis="#595959",cex.axis=1.5)
axis(2,col="#595959",col.ticks="#595959",col.axis="#595959",cex.axis=1.5)
box(col="#595959")
points(pca_atticus$V5~pca_atticus$V3,col="#71bf44",pch=19)
points(pca_gbenazzii$V5~pca_gbenazzii$V3,col="#f7941d",pch=19)
points(pca_ggrandii$V5~pca_ggrandii$V3,col="#be1e2d",pch=19)
points(pca_gmaretimi$V5~pca_gmaretimi$V3,col="#f9ed32",pch=19)
points(pca_hybridoNW$V5~pca_hybridoNW$V3,col="#f47f72",pch=19)
points(pca_hybridoSE$V5~pca_hybridoSE$V3,col="#ec008c",pch=19)
points(pca_whitei$V5~pca_whitei$V3,col="#92278f",pch=19)
points(pca_rredtenbacheri$V5~pca_rredtenbacheri$V3,col="#0599ce",pch=19)
points(pca_unknwon$V5~pca_unknwon$V3,col="#595959",pch=19)
points(pca_lynceorum$V5~pca_lynceorum$V3,col="#8b5e3c",pch=19)
points(pca$V5[pca$V1=="Bas"]~pca$V3[pca$V1=="Bas"],col="#71bf44",pch=19,cex=2)
points(pca$V5[pca$V1=="Bgibi"]~pca$V3[pca$V1=="Bgibi"],col="#f7941d",pch=19,cex=2)
points(pca$V5[pca$V1=="Bgigi"]~pca$V3[pca$V1=="Bgigi"],col="#be1e2d",pch=19,cex=2)
points(pca$V5[pca$V1=="Blm"]~pca$V3[pca$V1=="Blm"],col="#8b5e3c",pch=19,cex=2)
points(pca$V5[pca$V1=="Brsri"]~pca$V3[pca$V1=="Brsri"],col="#0599ce",pch=19,cex=2)
points(pca$V5[pca$V1=="SEhaplome"]~pca$V3[pca$V1=="SEhaplome"],col="#0599ce",pch=19,cex=2)
points(pca$V5[pca$V1=="NWhaplome"]~pca$V3[pca$V1=="NWhaplome"],col="#0599ce",pch=19,cex=2)
points(pca$V5[pca$V1=="Bwi"]~pca$V3[pca$V1=="Bwi"],col="#92278f",pch=19,cex=2)
dev.copy2eps(file="pc1v3.eps",width=8,height=8)


#species determination
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > -0.05) && ($3 < -0.035) && ($4 > -0.05) && ($4 < 0.05)) print $2}' > rossius_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > 0.05) && ($3 < 0.1) && ($4 > 0.073) && ($4 < 0.095)) print $2}' > benazzii_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > 0.07) && ($3 < 0.08) && ($4 > 0.06) && ($4 < 0.075)) print $2}' > maretimi_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > 0.045) && ($3 < 0.06) && ($4 > -0.07) && ($4 < -0.055)) print $2}' > grandii_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > 0) && ($3 < 0.015) && ($4 > 0.04) && ($4 < 0.06)) print $2}' > NWhybrids_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > 0.05) && ($3 < 0.06) && ($5 > -0.3) && ($5 < -0.1)) print $2}' > atticus_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > -0.005) && ($3 < 0.015) && ($5 > -0.06) && ($5 < -0.02)) print $2}' > lynceorum_list.txt
cat genome_norm_snp_an22_ql20_dp4_mm05_rm.eigenvec | awk '{if (($3 > -0.01) && ($3 < 0.01) && ($5 > 0.015) && ($5 < 0.03)) print $2}' > SEhybrids_list.txt


