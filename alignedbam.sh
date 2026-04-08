#!/bin/bash

# Set the reference genome and other parameters
REFERENCE=GCA_000018845.1.fasta
num_threads=64  # Adjust the number of threads as needed

module load minimap2
module load samtools

# Loop over all BAM files in the current directory
for USAM in *.bam; do
    # Define output SAM file name based on the BAM file name
    SAMFILE="${USAM%.bam}.sam"
    BAMFILE="aligned_bam/${USAM%.bam}_aligned.bam"
    SORTED_BAM="aligned_bam/${USAM%.bam}_sorted.bam"
	PILEUP_BED="pileups/${USAM%.bam}_pileup.bed"
    
    # Convert BAM to FASTQ, then align with minimap2 and generate SAM file
    samtools fastq -TMM,ML ${USAM} | minimap2 -x map-ont -a -t ${num_threads} -y --secondary=no ${REFERENCE} - > ${SAMFILE}
    
    # Convert SAM to BAM
    samtools view -Sb ${SAMFILE} > ${BAMFILE}
    
    # Sort the BAM file
    samtools sort -o ${SORTED_BAM} ${BAMFILE}
    
    # Index the sorted BAM file
    samtools index ${SORTED_BAM}
    
    # Clean up intermediate files (optional)
    rm ${BAMFILE}
    
	modkit pileup --filter-threshold A:0.95 --mod-thresholds a:0.95 $SORTED_BAM $PILEUP_BED
	#modkit pileup $SORTED_BAM $PILEUP_BED
    echo "Processed ${USAM} -> ${SORTED_BAM} and indexed."
done
