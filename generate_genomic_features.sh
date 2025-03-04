# This file contains commands for generating genomic feature values used in the paper, as annotated in gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13_genome_1kb.txt

# Input files:
#     gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13.tar.gz
#     (large raw files are not included but can be downloaded using the wget command below)
# Tools:
#     UCSC Genome Browser uses a variety of executables
#     http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
#     bedtools
#     https://bedtools.readthedocs.io/en/latest/
#     CrossMap
#     https://crossmap.readthedocs.io/en/latest/

# Define window size
window=1kb
window=1kb_flnk_10k # 10k around every 1kb
window=1kb_flnk_100k
window=1kb_flnk_1M


# ===== GC Content =====
# source data: 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
hgGcPercent -doGaps -win=1000 -file=genome_v31_${window}_GC.txt -verbose=0 hg38 hg38.2bit


# ===== Low-complexity Regions =====
# source data: https://academic.oup.com/bioinformatics/article/30/20/2843/2422145?login=false#supplementary-data
bedtools coverage -a hg38.chrom.${window}.bed -b btu356_LCR-hs38.bed > genome_v31_${window}_LCR.txt


# ===== Repeats =====
# source data: UCSC table browser — Group: Repeats — Track: RepeatMasker
repeat=SINE
repeat=LINE
bedtools coverage -a hg38.chrom.${window}.bed -b repeatMasker_hg38.${repeat}.bed > genome_v31_${window}_${repeat}.txt


# ===== Genomic Location =====
# source data: UCSC table browser — Group: Mapping and Sequencing — Track: Centromeres

bedtools closest -a hg38.chrom.${window}.mid.bed -b centromeres_hg38.bed -d > genome_v31_${window}_dist2cent.tmp.txt

bedtools closest -a hg38.chrom.${window}.mid.bed -b telomeres_hg38.bed -d > genome_v31_${window}_dist2telo.tmp.txt

python genomic_features13/normalize_dist2gaps_1kb.py


# ===== Recombination Rate =====
# source data: https://science.sciencemag.org/content/363/6425/eaau1043
parent=pat sex=male
parent=mat sex=female

sed "/#/d" genetic.map.final.${parent}.gor.txt | sed "/cMperMb/d" | awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2,$3,$4}'> genetic.map.final.${parent}.gor.bed

bedGraphToBigWig genetic.map.final.${parent}.gor.bed hg38.chrom.sizes genetic.map.final.${parent}.gor.bw

bigWigAverageOverBed genetic.map.final.${parent}.gor.bw hg38.chrom.${window}.bed genome_v31_${window}_recomb_${sex}.txt


# ===== Methylation =====
# source data: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81233/matrix/
mkdir methyl_sperm
wget -i Sperm_files.txt

for FILE_PREFIX in GSM2481659_scBS-ZHB-Sp1.Cmet GSM2481660_scBS-ZHB-Sp2.Cmet GSM2481661_scBS-ZHB-Sp3.Cmet GSM2481662_scBS-ZHB-Sp4.Cmet GSM2481663_scBS-ZHB-Sp5.Cmet GSM2481664_scBS-ZHB-Sp6.Cmet GSM2481665_scBS-ZHB-Sp7.Cmet GSM2481666_scBS-ZHB-Sp8.Cmet GSM2481667_scBS-ZHB-Sp9.Cmet GSM2481652_scBS-ZHB-Sp10.Cmet GSM2481653_scBS-ZHB-Sp11.Cmet GSM2481654_scBS-ZHB-Sp12.Cmet GSM2481655_scBS-ZHB-Sp13.Cmet GSM2481656_scBS-ZHB-Sp14.Cmet GSM2481657_scBS-ZHB-Sp15.Cmet GSM2481658_scBS-ZHB-Sp16.Cmet GSM2986301_scBS-hSp10.Cmet GSM2481579_scBS-hSp8.Cmet GSM2481580_scBS-hSp9.Cmet GSM2986302_scBS-hSP3.Cmet GSM2986303_scBS-hSP5.Cmet; do

    zcat $FILE_PREFIX.bed.gz | awk 'BEGIN {FS= "\t"; OFS= "\t"}; ($10 == "CpG") {print $1,$2-1,$2,$4, $5,$6,$8}' > $FILE_PREFIX.CpG.1.bed

    rm $FILE_PREFIX.bed.gz

    bedSort $FILE_PREFIX.CpG.1.bed $FILE_PREFIX.CpG.2.bed

    echo "chrom	start	end	strand	ttl	met	perc" > $FILE_PREFIX.CpG.bed

    cat $FILE_PREFIX.CpG.2.bed >> $FILE_PREFIX.CpG.bed

    rm $FILE_PREFIX.CpG.1.bed
    rm $FILE_PREFIX.CpG.2.bed

done


cat methyl_sperm/* | sed /chrom/d | awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2-1,$2,$3, $4,$5}' > $GSE81233_Sperm.1.bed

bedSort $GSE81233_Sperm.1.bed $GSE81233_Sperm.2.bed

bedtools merge -i $GSE81233_Sperm.2.bed -c 4,5,6 -o sum,sum,mean |  awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2,$3,$4,$5,$6,$5/$4}' > $GSE81233_Sperm.merged.mean2.bed

rm $GSE81233_Sperm.1.bed
rm $GSE81233_Sperm.2.bed

awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $1,$2,$3,$7}' GSE81233_Sperm.merged.mean2.bed > GSE81233_Sperm.merged.mean2.bedGraph
    
bedGraphToBigWig GSE81233_Sperm.merged.mean2.bedGraph hg19.chrom.sizes GSE81233_Sperm.merged.mean2.bw

.local/bin/CrossMap.py bigwig hg19ToHg38.over.chain.gz GSE81233_Sperm.merged.mean2.bw GSE81233_Sperm.merged.mean2.hg38

bigWigAverageOverBed misc/GSE81233_Sperm.merged.mean2.hg38.bw hg38.chrom.${window}.bed genome_v31_${window}_met_sperm.txt


# ===== CpG Island =====
# source data: UCSC table browser — Group: Regulation — Track: CpG Islands

awk 'BEGIN {FS= "\t"; OFS= "\t"}; {print $2,$3,$4}' cpgIslandExt.txt | sed /chrom/d > cpgIslandExt.bed

bedtools coverage -a hg38.chrom.${window}.bed -b cpgIslandExt.bed > genome_v31_${window}_CpGI.txt


# ===== Nucleosome Occupancy =====
# source data: 
wget https://www.encodeproject.org/files/ENCFF000VME/@@download/ENCFF000VME.bigWig
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz

.local/bin/CrossMap.py bigwig hg19ToHg38.over.chain.gz ENCFF000VME.bigWig ENCFF000VME.hg38

bigWigAverageOverBed ENCFF000VME.hg38.bw hg38.chrom.${window}.bed genome_v31_${window}_nucleosome.txt


# ===== DNM clusters =====
# source data: https://www.nature.com/articles/s41588-018-0071-6
# Supplementary Table 5

infile="Goldmann_18_S5_cDNMs_F.lft38.bed"
outfile="Goldmann_18_S5_cDNMs_paternal_flnk05M.lft38.txt"

infile="Goldmann_18_S5_cDNMs_M.lft38.bed"
outfile="Goldmann_18_S5_cDNMs_maternal_flnk05M.lft38.txt"

awk -F'\t' '{
    start = ($2 - 500000) < 0 ? 0 : ($2 - 500000);
    end = $3 + 500000;
    print $1, start, end
}' OFS='\t' "$infile" > "$outfile"

parent=paternal
parent=maternal
bedtools coverage -a hg38.chrom.${window}.bed -b Goldmann_18_S5_cDNMs_${parent}_flnk05M.lft38.txt > genome_v31_${window}_cDNM_${parent}_flnk05M.txt

