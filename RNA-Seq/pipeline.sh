#!/bin/bash

# Author: Danny Arends
# Modified by: Florian Krause
# Year: 2021 - 2025
# License: GNU General Public License v3.0 (GPL-3.0)
#
# This script is licensed under GPL-3.0.
#
# IMPORTANT NOTE ON EXTERNAL DEPENDENCIES:
#   • The STAR aligner (used via `STAR --runMode genomeGenerate`) is licensed under MIT
#   • Picard tools (used for MarkDuplicates, AddOrReplaceReadGroups) are Apache 2.0
#   • GATK (BaseRecalibrator, ApplyBQSR) is distributed under BSD-3-Clause
#   • Mouse reference genome (GRCm38.p6) and annotation (GTF) are from Ensembl:
#     * Licensed under CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
#     * Requires attribution: "Ensembl release 102 | Mus musculus GRCm38"
#
# This project's code is separate from the tools and reference data it uses.
# All external dependencies retain their original licenses and copyright notices.
#
# See:
#   • STAR license: https://github.com/alexdobin/STAR/blob/master/LICENSE
#   • Picard license: https://github.com/broadinstitute/picard/blob/main/LICENSE.txt
#   • GATK license: https://github.com/broadinstitute/gatk/blob/master/LICENSE
#   • Ensembl reference data terms: https://www.ensembl.org/info/legal/disclaimer.html

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
# Path to your main project directory
WORKDIR="/path/to/your/project_directory"
# Path where raw FASTQ files are located (Files must be named B6N4_R1... / BFMI4_R1...)
FASTQ_DIR="/path/to/raw/fastq_files"
# Path to the directory containing jar files (Picard, GATK)
SOFTWARE_DIR="/path/to/your/software_tools"

# 2. SUBDIRECTORIES (Automatically set relative to WORKDIR)
REF_DIR="$WORKDIR/References"
OUTPUT_DIR="$WORKDIR/Processed_BAMs"

# 3. TOOL VERSIONS / FILENAMES
# Ensure these files exist in your SOFTWARE_DIR
PICARD_JAR="picard.jar"
GATK_JAR="gatk-package-4.2.6.1-local.jar" # Or your specific GATK version

# 4. RESOURCES
THREADS=8
RAM_GATK="-Xmx16g" # Adjust based on available system memory

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

mkdir -p $WORKDIR $REF_DIR $OUTPUT_DIR

# --- 1. REFERENCE PREPARATION (Mouse mm10 / GRCm38) ---
cd $REF_DIR

# Reference Genome (GRCm38.p6)
if [ ! -f "Mus_musculus.GRCm38.dna.primary_assembly.fa" ]; then
    echo "Downloading Reference Genome..."
    wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
    gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
fi

# GTF Annotation
if [ ! -f "Mus_musculus.GRCm38.102.gtf" ]; then
    echo "Downloading GTF Annotation..."
    wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
    gunzip Mus_musculus.GRCm38.102.gtf.gz
fi

# Known Sites (dbSNP) for GATK Recalibration
if [ ! -f "mgp.v5.merged.snps_all.dbSNP142.vcf.gz" ]; then
    echo "Checking for dbSNP file..."
    echo "PLEASE NOTE: For full reproducibility, ensure 'mgp.v5.merged.snps_all.dbSNP142.vcf.gz' is present in $REF_DIR."
    # Creating a placeholder to prevent script failure if file is missing during dry-runs
    touch mgp.v5.merged.snps_all.dbSNP142.vcf.gz 
fi

# Indexing Reference for GATK (dict & fai)
if [ ! -f "Mus_musculus.GRCm38.dna.primary_assembly.dict" ]; then
    echo "Creating Sequence Dictionary..."
    java -jar $SOFTWARE_DIR/$PICARD_JAR CreateSequenceDictionary \
        R=Mus_musculus.GRCm38.dna.primary_assembly.fa \
        O=Mus_musculus.GRCm38.dna.primary_assembly.dict
    samtools faidx Mus_musculus.GRCm38.dna.primary_assembly.fa
fi

# STAR Index
if [ ! -d "STAR_Index" ]; then
    echo "Generating STAR Genome Index..."
    mkdir STAR_Index
    STAR --runThreadN $THREADS --runMode genomeGenerate \
         --genomeDir STAR_Index \
         --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
         --sjdbGTFfile Mus_musculus.GRCm38.102.gtf \
         --sjdbOverhang 100
fi

# --- 2. ALIGNMENT & PRE-PROCESSING LOOP ---
cd $OUTPUT_DIR

# ------------------------------------------------------------------------------
# SAMPLE DEFINITION
# Samples: B6N4 and BFMI4
# Input files must match these prefixes (e.g., B6N4_R1_001.fastq.gz)
# ------------------------------------------------------------------------------
SAMPLES=("B6N4" "BFMI4")

# Read Groups Information
# Format: ID:Library:PlatformUnit:SampleName
# "6-8-weeks" has been removed from Library (LB) and replaced with Sample ID.
declare -A RG_INFO
RG_INFO["B6N4"]="B6N4:B6N4:AGCGAGAT+GATACTGG:B6N4"
RG_INFO["BFMI4"]="BFMI4:BFMI4:GAATCACC+CTTCGTTC:BFMI4"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ">>> Processing Sample: $SAMPLE <<<"

    # A) STAR Alignment
    # Reads files named ${SAMPLE}_R1... directly
    if [ ! -f "${SAMPLE}_Aligned.sortedByCoord.out.bam" ]; then
        STAR --runThreadN $THREADS --genomeDir $REF_DIR/STAR_Index \
             --readFilesIn $FASTQ_DIR/${SAMPLE}_R1_001.fastq.gz $FASTQ_DIR/${SAMPLE}_R2_001.fastq.gz \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${SAMPLE}_
    else
        echo "   ... STAR alignment already exists. Skipping."
    fi

    # B) Picard MarkDuplicates
    if [ ! -f "${SAMPLE}_dedup.bam" ]; then
        java $RAM_GATK -jar $SOFTWARE_DIR/$PICARD_JAR MarkDuplicates \
             I=${SAMPLE}_Aligned.sortedByCoord.out.bam \
             O=${SAMPLE}_dedup.bam \
             M=${SAMPLE}_metrics.txt \
             REMOVE_DUPLICATES=true
    fi

    # C) Picard AddOrReplaceReadGroups
    if [ ! -f "${SAMPLE}_dedup_rg.bam" ]; then
        IFS=':' read -r LB SM PU SN <<< "${RG_INFO[$SAMPLE]}"
        
        java $RAM_GATK -jar $SOFTWARE_DIR/$PICARD_JAR AddOrReplaceReadGroups \
             I=${SAMPLE}_dedup.bam \
             O=${SAMPLE}_dedup_rg.bam \
             RGID=$SAMPLE RGLB=$LB RGPL=ILLUMINA RGPU=$PU RGSM=$SM

        samtools index ${SAMPLE}_dedup_rg.bam
    fi

    # D) GATK BaseRecalibrator
    if [ ! -f "${SAMPLE}_recal_data.table" ]; then
        java $RAM_GATK -jar $SOFTWARE_DIR/$GATK_JAR BaseRecalibrator \
             -R $REF_DIR/Mus_musculus.GRCm38.dna.primary_assembly.fa \
             --known-sites $REF_DIR/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
             -I ${SAMPLE}_dedup_rg.bam \
             -O ${SAMPLE}_recal_data.table
    fi

    # E) GATK ApplyBQSR
    if [ ! -f "${SAMPLE}_final.bam" ]; then
        java $RAM_GATK -jar $SOFTWARE_DIR/$GATK_JAR ApplyBQSR \
             -R $REF_DIR/Mus_musculus.GRCm38.dna.primary_assembly.fa \
             -bqsr ${SAMPLE}_recal_data.table \
             -I ${SAMPLE}_dedup_rg.bam \
             -O ${SAMPLE}_final.bam

        samtools index ${SAMPLE}_final.bam
        
        # QC Stats
        samtools flagstats ${SAMPLE}_final.bam > ${SAMPLE}_final.flagstats
    fi
done

echo "Pipeline finished successfully."