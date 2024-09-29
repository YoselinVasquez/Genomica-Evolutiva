# Genómica Evolutiva - Bioinformática y Ciencias Ómicas
una colección de códigos para la clase de genómica evolutiva

# codigo 1: Descargar sratoolkit.
```r
mkdir genomas;
grep ">" genomas.fasta ;
ls ;

# Para descargar sratoolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz

# Dar permiso de lectura, escritura, ejecución
chmod 777 stk.tar.gz

# Descomprimir
tar -vxzf stk.tar.gz

export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin
```

# Código 2: Descargar archivos sra., generar archivos de salida en formato fastq. y obtención de archivos en formato fastqc.
```r
prefetch -h 
prefetch --max-size 50G --option-file denv.txt

# Saca los archivos sra al directorio de trabajo
mv */*.sra .

# Eliminar directorios
rm -r ERR12389866/ ERR12543675/

# Generar archivos de salida en formato fastq
fasterq-dump --split-files *.sra 

# Comprimir archivos fastq
gzip *fastq

# Descargar archivos fastqc
fastqc *
```
# Código 3: Mapeado de secuencias
```r
# Descargar de NCBI
prefetch --max-size 50G --option-file accessions_mpox.txt
mv */*.sra .
fasterq-dump --split-files *.sra
gzip *fastq
fastqc *

#1 indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 4 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 4 ${prefix}.bam ;
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls ;

#4# extraer genoams consenso #
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 
ls ;
```
