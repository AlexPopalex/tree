# phase RAD reads of SEhybrids, NWhybrids, Blm, novel hybrids and simulated hybrids (bbmap/38.6.3)
for i in $(cat SEhybrids_list.txt); do bbsplit.sh threads=40 ref=Brsri_v3.fasta,Bgigi_v5_mod.fa in=$i-trimmed.fq.gz basename=bbsplit/$i-%-bbsplit.fq.gz refstats=bbsplit/$i-refstats ambiguous2=toss; done


# extract and plot competitive mapping results
for i in $(cat lynceorum_list.txt); do cat $i-refstats | grep -v name | grep Brsri_v3 | awk '{print $2}'; done > lynceorum_Brsri
for i in $(cat lynceorum_list.txt); do cat $i-refstats | grep -v name | grep Bgigi_v5 | awk '{print $2}'; done > lynceorum_Bgigi
for i in $(cat lynceorum_list.txt); do cat $i-refstats | grep -v name | grep Bas_cns | awk '{print $2}'; done > lynceorum_Bas
paste lynceorum_list.txt lynceorum_Brsri lynceorum_Bgigi lynceorum_Bas

for i in $(cat SEhybrids_list.txt); do cat $i-refstats | grep -v name | grep Brsri_v3 | awk '{print $2}'; done > SEhybrids_Brsri
for i in $(cat SEhybrids_list.txt); do cat $i-refstats | grep -v name | grep Bgigi_v5 | awk '{print $2}'; done > SEhybrids_Bgigi
paste SEhybrids_list.txt SEhybrids_Brsri SEhybrids_Bgigi > SEhybrids_bbsplit.txt

for i in $(cat NWhybrids_list.txt); do cat $i-refstats | grep -v name | grep Brsri_v3 | awk '{print $2}'; done > NWhybrids_Brsri
for i in $(cat NWhybrids_list.txt); do cat $i-refstats | grep -v name | grep Bgibi_cns | awk '{print $2}'; done > NWhybrids_Bgigi
paste NWhybrids_list.txt NWhybrids_Brsri NWhybrids_Bgigi > NWhybrids_bbsplit.txt

boxplot(SE$V2,SE$V3,col="white",cex.lab=1.5,col.lab="#595959",names=c("Brsri","Bgigi"),outline = F,ylim=c(0,100),ylab="% reads mapping",xlab="to genome:",main="",cex.main=1.5,col.ticks="#595959",col.axis="#595959",cex.axis=1.4,frame=F,col.main="#595959")
stripchart(SE$V2, method = "jitter", col = "#0599ce", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(1),cex=1.3)
stripchart(SE$V3, method = "jitter", col = "#be1e2d", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(2),cex=1.3)
boxplot(NW$V2,NW$V3,col="white",cex.lab=1.5,col.lab="#595959",names=c("Brsri","Bgibi"),outline = F,ylim=c(0,100),ylab="% reads mapping",xlab="to genome:",main="",cex.main=1.5,col.ticks="#595959",col.axis="#595959",cex.axis=1.4,frame=F,col.main="#595959")
stripchart(NW$V2, method = "jitter", col = "#0599ce", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(1),cex=1.3)
stripchart(NW$V3, method = "jitter", col = "#f7941d", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(2),cex=1.3)
boxplot(Blm$V2,Blm$V3,Blm$V4,col="white",cex.lab=1.5,col.lab="#595959",names=c("Brsri","Bgigi","Bas"),outline = F,ylim=c(0,100),ylab="% reads mapping",xlab="to genome:",main="",cex.main=1.5,col.ticks="#595959",col.axis="#595959",cex.axis=1.4,frame=F,col.main="#595959")
stripchart(Blm$V2, method = "jitter", col = "#0599ce", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(1),cex=1.3)
stripchart(Blm$V3, method = "jitter", col = "#be1e2d", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(2),cex=1.3)
stripchart(Blm$V4, method = "jitter", col = "#71bf44", vertical = TRUE, pch = 19,ylim=c(0,100),add=T,at=c(3),cex=1.3)
dev.copy2eps(file="compmap.eps",width=12,height=8)


# calculate means of alternative allele depths per sample of hybridogen 1/B. whitei based on vcf used for species determination
for (( i = 10; i <= 111; i++ )); do echo "zcat genome_norm_snp_an22_ql20_dp4_mm05_rm_SE.vcf.gz | awk '{if (\$$i ~ /0\/1/) print \$$i}' > ${i}_aad"; done 
for (( i = 10; i <= 111; i++ )); do cat ${i}_aad | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print ($1+$2),($2/($1+$2)); else print "NA","NA"}' | awk '{ total += $2 } END { if (NR>0) print total/NR; else print "NA"}'; done > aad_SE_means

# calculate means of minor allele depths per sample of Blm based on vcf used for species determination
for (( i = 10; i <= 38; i++ )); do echo "zcat genome_norm_snp_an22_ql20_dp4_mm05_rm_Blm.vcf.gz | awk '{if (\$$i ~ /0\/1/) print \$$i}' > ${i}_aad_Blm"; done 
for (( i = 10; i <= 38; i++ )); do cat ${i}_aad_Blm | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print $1,$2; else print "NA","NA"}' | sed '/NA/d' | awk '{if ($1<=$2) print ($1+$2),$1/($1+$2); else print ($1+$2),$2/($1+$2)}' | awk '{ total += $2 } END { if (NR>0) print total/NR; else print "NA"}'; done > mad_Blm_means



for (( i = 10; i <= 111; i++ )); do cat ${i}_aad | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print ($2/($1+$2)); else print "NA","NA"}' > ${i}_SE_aad_dist; done
for (( i = 10; i <= 111; i++ )); do cat ${i}_mad_SE | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print $1,$2; else print "NA","NA"}' | sed '/NA/d' | awk '{if ($1<=$2) print $1/($1+$2); else print $2/($1+$2)}' > ${i}_SE_mad_dist; done
for (( i = 10; i <= 38; i++ )); do cat ${i}_aad_Blm | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print ($2/($1+$2)); else print "NA","NA"}' > ${i}_Blm_aad_dist; done
for (( i = 10; i <= 38; i++ )); do cat ${i}_aad_Blm | awk 'BEGIN {FS=OFS=":"} {print $3,$4}' | awk 'BEGIN {FS=OFS=":"} {print $2}' | awk -F ',' '{if (($1+$2) > 0) print $1,$2; else print "NA","NA"}' | sed '/NA/d' | awk '{if ($1<=$2) print $1/($1+$2); else print $2/($1+$2)}' > ${i}_Blm_mad_dist; done

paste *_Blm_mad_dist > Blm_mad_dist_pasted
paste *_Blm_aad_dist > Blm_aad_dist_pasted
paste *_SE_mad_dist > SE_mad_dist_pasted
paste *_SE_aad_dist > SE_aad_dist_pasted

read.table("Blm_aad_dist_pasted",fill=T)-> Blm_aad_dist
read.table("Blm_mad_dist_pasted",fill=T)-> Blm_mad_dist
read.table("SE_aad_dist_pasted",fill=T)-> SE_aad_dist
read.table("SE_mad_dist_pasted",fill=T)-> SE_mad_dist

plot(density(na.omit(Blm_aad_dist$V1),adjust=3),type="n")
for (( i = 1; i <= 29; i++ )); do echo "lines(density(na.omit(Blm_aad_dist\$V${i}),adjust=3))"; done
plot(density(na.omit(Blm_mad_dist$V1),adjust=3),type="n")
for (( i = 1; i <= 29; i++ )); do echo "lines(density(na.omit(Blm_mad_dist\$V${i}),adjust=3))"; done
plot(density(na.omit(SE_aad_dist$V1),adjust=3),type="n")
for (( i = 1; i <= 102; i++ )); do echo "lines(density(na.omit(SE_aad_dist\$V${i}),adjust=3))"; done
plot(density(na.omit(SE_mad_dist$V1),adjust=3),type="n")
for (( i = 1; i <= 102; i++ )); do echo "lines(density(na.omit(SE_mad_dist\$V${i}),adjust=3))"; done






# map phased haplosets of hybrids and unphased reads of parental species to the B. r. redtenbacheri and B. g. grandii genomes, respectively (bwa/0.7.17, samtools/1.15.1) 
bwa index Brsri_v3.fasta
for i in $(cat ../SEhybrids_list.txt); do
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
Brsri_v3.fasta \
${i}-Brsri_v3-bbsplit.fq.gz > rossius/${i}_Brsri.sam
samtools view -@ 40 -b -h rossius/${i}_Brsri.sam > rossius/${i}_Brsri.bam
rm -f rossius/${i}_Brsri.sam
samtools sort -@ 40 rossius/${i}_Brsri.bam -o rossius/${i}_Brsri_sorted.bam
rm -f rossius/${i}_Brsri.bam
samtools view -h rossius/${i}_Brsri_sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -h -F 2308 -q 20 -@ 40 -o rossius/${i}_Brsri_sorted_uniq.bam
rm -f rossius/${i}_Brsri_sorted.bam
samtools index -@ 40 rossius/${i}_Brsri_sorted_uniq.bam
done

# call variants of all rossius haplosets (bcftools/1.15.1)
bcftools mpileup --threads 40 -Ou --annotate FORMAT/DP,FORMAT/AD -r Brsri_v3_scf1 -f Brsri_v3.fasta -b full_list_rad.txt | bcftools call --threads 40 -g 0 -m -Oz -o scf1.vcf.gz
bcftools index scf1.vcf.gz

bcftools concat --threads 20 -f vcf_list -O z -o genome.vcf.gz


# filter variants (bcftools/1.15.1)
bcftools norm --threads 40 -f Brsri_v3.fasta -d all -O z -o genome_norm.vcf.gz genome.vcf.gz
bcftools view --exclude-types indels --threads 40 -O z genome_norm.vcf.gz > genome_norm_snp.vcf.gz
bcftools view -m 1 -M 2 --threads 40 -O z genome_norm_snp.vcf.gz > genome_norm_snp_an12.vcf.gz
bcftools filter -e 'QUAL<20' --threads 40 -O z genome_norm_snp_an12.vcf.gz > genome_norm_snp_an12_ql20.vcf.gz
bcftools filter -e 'FMT/DP<4' --set-GTs . --threads 40 -O z genome_norm_snp_an12_ql20.vcf.gz > genome_norm_snp_an12_ql20_dp4.vcf.gz
bcftools index genome_norm_snp_an12_ql20_dp4.vcf.gz


# calculate per-genotype sequencing depth and remove sites with depth above threshold, i.e. (median(depth) + 2 * standard deviation(depth)) * #samples (bcftools/1.15.1)
zcat genome_norm_snp_an12_ql20_dp4.vcf.gz | awk '{if (($8 ~ /AN/) && ($8 ~ /DP/)  && ($8 !~ /MinDP/)) print $8}' | grep -o '[^;]*DP=[^;]*' | sed 's/DP=//g' > DP
zcat genome_norm_snp_an12_ql20_dp4.vcf.gz | awk '{if (($8 ~ /AN/) && ($8 ~ /DP/)  && ($8 !~ /MinDP/)) print $8}' | grep -o '[^;]*AN=[^;]*' | sed 's/AN=//g' > AN
paste DP AN > DPAN
cat DPAN | awk '{if ($2!="0") print $1,$2,$1/$2}' > DPANdiv
bcftools filter -e 'INFO/DP>12546' --threads 40 -O z genome_norm_snp_an12_ql20_dp4.vcf.gz > genome_norm_snp_an12_ql20_dp4_dp12546.vcf.gz


# continue filtering (bcftools/1.15.1)
bcftools filter -e 'FMT/GT="het"' --set-GTs . --threads 40 -O z genome_norm_snp_an12_ql20_dp4_dp12546.vcf.gz > genome_norm_snp_an12_ql20_dp4_dp12546_hetgt.vcf.gz
bcftools index genome_norm_snp_an12_ql20_dp4_dp12546_hetgt.vcf.gz
bcftools convert --gvcf2vcf --threads 40 -f Brsri_v3.fasta -O z -o genome_norm_snp_an12_ql20_dp4_dp12546_hetgt_conv.vcf.gz genome_norm_snp_an12_ql20_dp4_dp12546_hetg$
bcftools index genome_norm_snp_an12_ql20_dp4_dp12546_hetgt_conv.vcf.gz


#construct pseudoreferences per sample with Ns for missing genotypes and sites (samtools/1.15.1, bcftools/1.15.1)
for i in $(cat xaa_list); do samtools faidx Brsri_v3.fasta -r scf_list | bcftools consensus -M N -a N -H I -s $i -p $i -o /scratch/abrandt1/pseudorefs_rossius/$i.fa genome_norm_snp_an12_ql20_dp4_dp12546_hetgt_conv.vcf.gz; done


#convert psequdoreferences into alignments per scaffold
cat *.fa > genomes.fa
cat genomes.fa | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > genomes_ni.fa

grep "scf1\b" -A 1 --no-group-separator genomes_ni.fa > genomes_scf1.fa
sed -i 's/Brsri_v3_scf1//g' genomes_scf1.fa


#remove sites containing a percentage of N's (trimal/1.4)
sed 's/N/\-/g;s/n/\-/g' genomes_scf1.fa > genomes_scf1_gaps.fa
/work/FAC/FBM/DEE/tschwand/default/abrandt/Software/trimal/source/trimal -gt 1 -in genomes_scf1_gaps.fa -out genomes_scf1_1.fa
/work/FAC/FBM/DEE/tschwand/default/abrandt/Software/trimal/source/trimal -gt 0.75 -in genomes_scf1_gaps.fa -out genomes_scf1_075.fa
/work/FAC/FBM/DEE/tschwand/default/abrandt/Software/trimal/source/trimal -gt 0.5 -in genomes_scf1_gaps.fa -out genomes_scf1_05.fa
sed -i 's/\-/N/g' genomes_scf1_075.fa
sed -i 's/\-/N/g' genomes_scf1_05.fa


#calculate tree (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_scf1_1.fa -m MFP -B 1000 --redo


#calculate constrained tree with best fitting model for topology testing (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_075.fasta -g Blm_gigi_topo_h1.treefile -m TVM+F+I+I+R7 --prefix Blm_gigi_h1 --redo

#test best fitting ML tree vs the tree with the constrained topology (iqtree/2.2.0.5)
iqtree2 -nt 40 -s genomes_075.fasta -te genomes_075.fasta.treefile -z Blm_gigi_1.treefiles -zb 10000 -zw --prefix Blm_1_zw --redo

