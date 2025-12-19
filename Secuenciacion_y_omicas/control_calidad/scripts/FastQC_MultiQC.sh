# Ejecutar FastQC sobre los archivos FastQ
fastqc fastq/*.fastq.gz -o qc/

# Generar informe con MultiQC
multiqc qc/ -o qc/
