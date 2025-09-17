### transition 2

# concatenate and rename 18 largest scaffolds of Brsri (scfs 4 and 9 split in reference genome) according to their synteny to the 17 longest Bgigi scaffolds using geneious v2022.0.1


# competitively map RAD reads of unfertilized eggs of F2 females and only retain primary alignments and reads in proper pair using bwa v0.7.17 and samtools v1.19.2
bwa index Bgigi_Brsri_psscfs.fasta
for i in $(cat list); do
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
Bgigi_Brsri_psscfs.fasta \
${i}_R1.fq.gz \
${i}_R2.fq.gz > psscfs/${i}.sam
samtools view -@ 40 -b -h psscfs/${i}.sam > psscfs/${i}.bam
rm -f psscfs/${i}.sam
samtools sort -@ 40 psscfs/${i}.bam -o psscfs/${i}_sorted.bam
rm -f psscfs/${i}.bam
samtools view -h psscfs/${i}_sorted.bam | samtools view -h -f 2 -F 2308 -@ 40 -o psscfs/${i}_sorted_prim.bam
rm -f psscfs/${i}_sorted.bam
samtools index -@ 40 psscfs/${i}_sorted.bam
samtools index -@ 40 psscfs/${i}_sorted_prim.bam
done


# calculate sliding window coverage and adjust window size for Brsri to keep synteny to Bgigi when plotted using tinycov v0.4.0
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim.txt --res 10000000 --skip 100000; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf1.txt -w Brsri_v3_psscf1 --res 8611095 --skip 86111; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf2.txt -w Brsri_v3_psscf2 --res 7241043 --skip 72410; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf3.txt -w Brsri_v3_psscf3 --res 8885625 --skip 88856; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf4.txt -w Brsri_v3_psscf4 --res 9786708 --skip 97867; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf5.txt -w Brsri_v3_psscf5 --res 11242100 --skip 112421; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf6.txt -w Brsri_v3_psscf6 --res 9948955 --skip 99490; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf7.txt -w Brsri_v3_psscf7 --res 10352514 --skip 103525; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf8.txt -w Brsri_v3_psscf8 --res 10073029 --skip 100730; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf9.txt -w Brsri_v3_psscf9 --res 10183801 --skip 101838; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf10.txt -w Brsri_v3_psscf10 --res 11045952 --skip 110460; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf11.txt -w Brsri_v3_psscf11 --res 9934040 --skip 99340; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf12.txt -w Brsri_v3_psscf12 --res 10399645 --skip 103996; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf13.txt -w Brsri_v3_psscf13 --res 10853505 --skip 108535; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf14.txt -w Brsri_v3_psscf14 --res 10962567 --skip 109626; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf15.txt -w Brsri_v3_psscf15 --res 13517296 --skip 135173; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf16.txt -w Brsri_v3_psscf16 --res 10924266 --skip 109243; done
for i in $(cat list); do tinycov covplot ${i}_sorted_prim.bam --text ${i}_covplot_prim_corr_psscf17.txt -w Brsri_v3_psscf17 --res 11092337 --skip 110923; done
for i in $(cat list); do cat ${i}_covplot_prim_corr_psscf1.txt ${i}_covplot_prim_corr_psscf2.txt ${i}_covplot_prim_corr_psscf3.txt ${i}_covplot_prim_corr_psscf4.txt ${i}_covplot_prim_corr_psscf5.txt ${i}_covplot_prim_corr_psscf6.txt ${i}_covplot_prim_corr_psscf7.txt ${i}_covplot_prim_corr_psscf8.txt ${i}_covplot_prim_corr_psscf9.txt ${i}_covplot_prim_corr_psscf10.txt ${i}_covplot_prim_corr_psscf11.txt ${i}_covplot_prim_corr_psscf12.txt ${i}_covplot_prim_corr_psscf13.txt ${i}_covplot_prim_corr_psscf14.txt ${i}_covplot_prim_corr_psscf15.txt ${i}_covplot_prim_corr_psscf16.txt ${i}_covplot_prim_corr_psscf17.txt > ${i}_covplot_prim_Brsri_corr.txt; done
for i in $(cat list); do grep Bgigi ${i}_covplot_prim.txt | sed 's/Bgigi_v5_scf//g' > ${i}_covplot_prim_Bgigi.txt; done
for i in $(cat list); do sed -i 's/Brsri_v3_psscf//g' ${i}_covplot_prim_Brsri_corr.txt; done
for i in $(cat list); do nl ${i}_covplot_prim_Bgigi.txt > ${i}_covplot_prim_Bgigi_nl.txt; done
for i in $(cat list); do nl ${i}_covplot_prim_Brsri_corr.txt > ${i}_covplot_prim_Brsri_corr_nl.txt; done


# extract reads mapping to Brsri using samtools
for i in $(cat ../list); do samtools view ${i}_sorted_prim.bam | awk '{if ($3 ~ /Brsri/) print $1}' | sort -n -k 1,1 | uniq > ${i}_ros_reads; done
for i in $(cat ../list); do zgrep -f ${i}_ros_reads -A 3 ../../eggs_09_23/${i}_R1.fq.gz | bgzip - > ${i}_R1_ros.fq.gz; done
for i in $(cat ../list); do zgrep -f ${i}_ros_reads -A 3 ../../eggs_09_23/${i}_R2.fq.gz | bgzip - > ${i}_R2_ros.fq.gz; done


# generate pseudoreferences of males (grandfathers of the F2 females) and the hybridogenetic haplome for competitive mapping (shown only for the males here) using bwa, samtools and bcftools v1.21
bwa index Brsri_v3.fasta
bwa mem -t 40 -R "@RG\tID:F0males\tSM:F0males\tLB:F0males\tPL:Illumina" \
Brsri_v3.fasta \
../reads/F0males_R1_val_1.fq.gz \
../reads/F0males_R2_val_2.fq.gz > F0males_Brsri.sam
samtools view -@ 40 -b -h F0males_Brsri.sam > F0males_Brsri.bam
rm -f F0males_Brsri.sam
samtools sort -@ 40 F0males_Brsri.bam -o F0males_Brsri_sorted.bam
rm -f F0males_Brsri.bam

bcftools mpileup --threads 40 -Ou --annotate FORMAT/DP -f Brsri_v3.fasta F0males_Brsri_sorted.bam | bcftools call --threads 40 -m -Oz -o F0males_Brsri.vcf.gz
bcftools index F0males_Brsri.vcf.gz
bcftools consensus -H I -f Brsri_v3.fasta F0males_Brsri.vcf.gz > F0males_Brsri_HI.fa


# competitively map RAD reads of unfertilized eggs of F2 females to pseudoreference of males and hybridogenetic haplome and only retain uniquely mapping reads in proper pair using bwa and samtools
bwa index SEhaplome_F0males_psscfs.fa
for i in $(cat ../list); do
bwa mem -t 40 -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:Illumina" \
SEhaplome_F0males_psscfs.fa \
${i}_R1_ros.fq.gz \
${i}_R2_ros.fq.gz > psscfsros/${i}.sam
samtools view -@ 40 -b -h psscfsros/${i}.sam > psscfsros/${i}.bam
rm -f psscfsros/${i}.sam
samtools sort -@ 40 psscfsros/${i}.bam -o psscfsros/${i}_sorted.bam
rm -f psscfsros/${i}.bam
samtools view -h psscfsros/${i}_sorted.bam | samtools view -h -f 2 -F 2308 -@ 40 -o psscfsros/${i}_sorted_prim.bam
rm -f psscfsros/${i}_sorted.bam
samtools index -@ 40 psscfsros/${i}_sorted.bam
samtools index -@ 40 psscfsros/${i}_sorted_prim.bam
done

for i in $(cat ../list); do
samtools view -h psscfsros/${i}_sorted_prim.bam -@ 40 | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -h -q 20 -@ 40 -o psscfsros/${i}_sorted_prim_uniq.bam
samtools index psscfsros/${i}_sorted_prim_uniq.bam -@ 40
done

for i in $(cat ../../list); do samtools view ${i}_sorted_prim_uniq.bam | awk '{if ($3 ~ /SEhaplome/) print $1}' | sort -n -k 1,1 | uniq > ${i}_SEhaplome_reads; done
for i in $(cat ../../list); do samtools view ${i}_sorted_prim_uniq.bam | awk '{if ($3 ~ /F0males/) print $1}' | sort -n -k 1,1 | uniq > ${i}_F0males_reads; done
for i in $(cat ../../list); do cat ${i}_SEhaplome_reads | wc -l; done > SEhaplome_readnums
for i in $(cat ../../list); do cat ${i}_F0males_reads | wc -l; done > F0males_readnums
paste ../../list SEhaplome_readnums F0males_readnums | awk '{if ($2+$3>0) print$0,$3/($2+$3); else print $0,"NA"}' > readnums


# calculate sliding window coverage and sort scaffolds to generate synteny with Bgigi using tinycov
for i in $(cat ../../list); do tinycov covplot ${i}_sorted_prim_uniq.bam --text ${i}_covplot_prim_uniq.txt --res 10000000 --skip 100000; done

for i in $(cat ../../list); do grep SEhaplome ${i}_covplot_prim_uniq.txt > ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do grep F0males ${i}_covplot_prim_uniq.txt > ${i}_covplot_prim_uniq_F0males.txt; done

for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf1\t/SEhaplome_psscf1\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf1\t/F0males_psscf1\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf5\t/SEhaplome_psscf2\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf5\t/F0males_psscf2\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf3\t/SEhaplome_psscf3\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf3\t/F0males_psscf3\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf2\t/SEhaplome_psscf4\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf2\t/F0males_psscf4\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf4_2\t/SEhaplome_psscf5\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf4_2\t/F0males_psscf5\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf4_1\t/SEhaplome_psscf6\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf4_1\t/F0males_psscf6\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf6\t/SEhaplome_psscf7\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf6\t/F0males_psscf7\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf8\t/SEhaplome_psscf8\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf8\t/F0males_psscf8\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf7\t/SEhaplome_psscf9\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf7\t/F0males_psscf9\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf9_2\t/SEhaplome_psscf10\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf9_2\t/F0males_psscf10\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf9_1\t/SEhaplome_psscf11\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf9_1\t/F0males_psscf11\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf10\t/SEhaplome_psscf12\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf10\t/F0males_psscf12\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf12\t/SEhaplome_psscf13\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf12\t/F0males_psscf13\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf16\t/SEhaplome_psscf14\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf16\t/F0males_psscf14\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf13\t/SEhaplome_psscf15\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf13\t/F0males_psscf15\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf15\t/SEhaplome_psscf16\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf15\t/F0males_psscf16\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf14\t/SEhaplome_psscf17\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf14\t/F0males_psscf17\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf11\t/SEhaplome_psscf18\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf11\t/F0males_psscf18\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf17\t/SEhaplome_psscf19\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf17\t/F0males_psscf19\t/g' ${i}_covplot_prim_uniq_F0males.txt; done
for i in $(cat ../../list); do sed -i 's/SEhaplome_cns_scf18\t/SEhaplome_psscf20\t/g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_cns_scf18\t/F0males_psscf20\t/g' ${i}_covplot_prim_uniq_F0males.txt; done

for i in $(cat ../../list); do sed -i 's/SEhaplome_psscf//g' ${i}_covplot_prim_uniq_SEhaplome.txt; done
for i in $(cat ../../list); do sed -i 's/F0males_psscf//g' ${i}_covplot_prim_uniq_F0males.txt; done

for i in $(cat ../../list); do cat ${i}_covplot_prim_uniq_SEhaplome.txt | sort -n -k 1,1 -k 2,2 | nl > ${i}_covplot_prim_uniq_SEhaplome_nl.txt; done
for i in $(cat ../../list); do cat ${i}_covplot_prim_uniq_F0males.txt | sort -n -k 1,1 -k 2,2 | nl > ${i}_covplot_prim_uniq_F0males_nl.txt; done


