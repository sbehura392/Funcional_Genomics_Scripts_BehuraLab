# =============================================================================
# Pipeline: Splice Variant and Differential Exon Usage Analysis
#
# Reference:
# Islam M, Behura SK. Molecular Regulation of Fetal Brain Development in
# Inbred and Congenic Mouse Strains Differing in Longevity. Genes (Basel).
# 2024;15(5):604.
#
# Author: Susanta Behura
#
# Description:
# This pipeline processes paired-end RNA-seq data to identify splice junction
# usage and generate exon–junction count tables for differential exon usage
# analysis.
#
# Workflow Summary:
#   1. Perform quality control of raw RNA-seq reads
#   2. Build a STAR genome index and run first-pass alignment
#   3. Collect splice junction files from all samples
#   4. Build a second-pass STAR index using detected splice junctions
#   5. Run second-pass alignment for improved junction discovery
#   6. Merge and filter splice junction data across samples
#   7. Extract exon coordinates and annotation from the reference GTF
#   8. Intersect exon regions with splice junctions
#   9. Format exon–junction counts for downstream edgeR analysis
#
# Requirements:
#   - fastp
#   - STAR
#   - bedtools
#   - awk
#   - sed
#   - sort
#
# Input:
#   - Paired-end FASTQ files
#   - Reference genome FASTA
#   - Gene annotation GTF
#
# Output:
#   - Filtered splice junction BED file
#   - Exon annotation BED file
#   - Exon–junction count table for R / edgeR analysis
# =============================================================================


# =============================================================================
# Step 1: Perform Read Quality Control with fastp
# =============================================================================

module load fastp

for f in $(ls *_1.fastq | sed 's/_1.fastq//g'); do
    fastp \
        -i ${f}_1.fastq \
        -I ${f}_2.fastq \
        -o ${f}_1.qc.fq.gz \
        -O ${f}_2.qc.fq.gz
done


# =============================================================================
# Step 2: Build STAR Genome Index and Run First-Pass Alignment
# =============================================================================

module load star

# Generate STAR genome index
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir ./Index \
    --genomeFastaFiles Mus_musculus.GRCm39.dna.toplevel.fa \
    --sjdbGTFfile Mus_musculus.GRCm39.110.gtf

# Perform first-pass alignment
for f in $(ls *_1.qc.fq.gz | sed 's/_1.qc.fq.gz//g'); do
    STAR \
        --runThreadN 3 \
        --genomeDir ./Index \
        --readFilesIn ${f}_1.qc.fq.gz ${f}_2.qc.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outFileNamePrefix ${f}.
done


# =============================================================================
# Step 3: Collect First-Pass Splice Junction Files
# =============================================================================

mkdir -p allSpliceSitesReadCounts
mv *.SJ.out.tab allSpliceSitesReadCounts/


# =============================================================================
# Step 4: Build Second-Pass STAR Index Using Detected Junctions
# =============================================================================

STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir ./Index2 \
    --genomeFastaFiles Mus_musculus.GRCm39.dna.toplevel.fa \
    --sjdbGTFfile Mus_musculus.GRCm39.110.gtf \
    --sjdbFileChrStartEnd allSpliceSitesReadCounts/*.SJ.out.tab


# =============================================================================
# Step 5: Run Second-Pass Alignment
# =============================================================================

for f in $(ls *_1.qc.fq.gz | sed 's/_1.qc.fq.gz//g'); do
    STAR \
        --runThreadN 8 \
        --genomeDir ./Index2 \
        --readFilesIn ${f}_1.qc.fq.gz ${f}_2.qc.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outFileNamePrefix 2pass_${f}.
done


# =============================================================================
# Step 6: Collect Second-Pass Splice Junction Files
# =============================================================================

mkdir -p allSpliceSitesReadCounts_2pass
mv *.SJ.out.tab allSpliceSitesReadCounts_2pass/

cd allSpliceSitesReadCounts_2pass/


# =============================================================================
# Step 7: Merge and Filter Splice Junction Data Across Samples
# =============================================================================
#
# Output fields:
#   chromosome, start, end, strand, read count, sample ID, splice type
#
# STAR SJ.out.tab columns used here:
#   $1 = chromosome
#   $2 = intron start
#   $3 = intron end
#   $5 = strand
#   $6 = splice junction annotation (0 = novel, >0 = annotated)
#   $7 = uniquely mapped read count
#   $8 = multi-mapped read count
#   $9 = maximum spliced alignment overhang
#

awk -v OFS='\t' '{print $0, FILENAME}' *.SJ.out.tab | \
    sed 's/2pass_//g' | \
    sed 's/.SJ.out.tab//g' | \
    awk -v OFS='\t' '{
        if ($6=="0") print $0, "Novel_splice";
        else print $0, "Known_splice";
    }' | \
    awk -v OFS='\t' '{
        if ($8=="0" && $9 > 20) {
            print $1, $2, $3, $5, $7, $11, $10
        }
    }' | \
    sort -k1,1 -k2,2n > splice_junctions_mapped_reads.bed


# =============================================================================
# Step 8: Extract Exon Annotation from the GTF File
# =============================================================================

awk -F'[\t;]' -v OFS='\t' '{
    if ($3=="exon")
        print $1, $4, $5, $9, $11, $14, $17, $20
}' ../Mus_musculus.GRCm39.110.gtf | \
    sed 's/gene_id//g' | \
    sed 's/transcript_id//g' | \
    sed 's/gene_name//g' | \
    sed 's/transcript_name//g' | \
    sed 's/exon_id //g' | \
    sed 's/"//g' | \
    awk -v OFS='\t' '{
        if ($8 ~ /^ENSMUSE/)
            print $1, $2, $3, $8"_"$5"_"$4"_"$6"_"$7
    }' | \
    sort -k1,1 -k2,2n > splice_annotation_info.bed


# =============================================================================
# Step 9: Intersect Exons with Splice Junctions
# =============================================================================

module load bedtools2/bedtools2-2.27.1

bedtools intersect \
    -a splice_annotation_info.bed \
    -b splice_junctions_mapped_reads.bed \
    -wa -wb > exon_and_junctions


# =============================================================================
# Step 10: Format Output for edgeR-Based Exon Usage Analysis
# =============================================================================

awk -v OFS='\t' '{
    print $10"_type"$8":"$5":"$6":"$7":"$4, $9, $11
}' exon_and_junctions > exon_and_junctions_for_R


# =============================================================================
# Notes:
# =============================================================================
# - STAR two-pass alignment improves sensitivity for known and novel splice
#   junction detection.
# - Junction filtering criteria (for example, overhang >20 and no multi-mapped
#   reads) can be adjusted depending on sequencing depth and data quality.
# - The final output is formatted for downstream exon-level differential usage
#   analysis in R, including edgeR-based workflows.
# - Additional filtering by read support, transcript confidence, or annotation
#   class may improve robustness in downstream analyses.
# =============================================================================