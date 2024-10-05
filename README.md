# Genómica Evolutiva - Bioinformática y Ciencias Ómicas
Una colección de códigos para la clase de genómica evolutiva.

# codigo 1: Descargar sratoolkit.
```r
# Para crear un directorio, usar el comado mkdir seguido del nombre del directorio
mkdir genomas;

# para ver los encabezados de los archivos .fa
grep ">" *.fa ;

#Para ver todos los directorios y archivos que se encuentran en una ambiente de trabajo
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
# Crear directorios con archivos ".sra"
prefetch --max-size 50G --option-file denv.txt

# Saca los archivos sra. al directorio de trabajo
mv */*.sra .

# Eliminar directorios
rm -r ERR12389866/ ERR12543675/

# Generar archivos de salida en formato fastq (1 forward y 1 reverse).
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

#1 Indexar el genoma de referencia#
bwa index reference.fasta ;

# Ejecutar #2# y #3# juntos.
#2# Preparar las instrucciones generales#
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

#4# Extraer genoma consenso. Se generan archivos "fa." y "qual.txt"
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 
ls ;
```
# Código 4: Instalación y análisis Prokka e intalación de Artemis
```r
# Instalar prokka
conda create -n prokka ;
conda activate prokka ;
conda install bioconda/label/cf201901::prokka

# Otorgar permisos al directorio de trabajo
chmod 777 (nombre del directorio)
```
# Analisis con prokka
mkdir annotation/ ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

#Para salir del environment 
conda deactivate ;

# Instalar artemis
conda create -n art
conda activate art
conda install bioconda::artemis
```
# Código 5: Observación de los resultados en Artemis 
```r
# Es necesario contar con archivos ".fa" y ".gff"
# Primero: crear un ambiente para artemis
conda activate art

# Segundo: entrar a la interfaz colocando "art" en la terminal.

# Tercero: seleccionar el archivo ".gff"
File --> open file manager --> cargar achivo con extensión ".gff".
