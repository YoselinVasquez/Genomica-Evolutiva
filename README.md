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
# Desacargar bwa, samtools, ivar, igv.
# Descargar de NCBI genoma de referencia en formato fasta.
prefetch --max-size 50G --option-file accessions_mpox.txt
mv */*.sra .
fasterq-dump --split-files *.sra
gzip *fastq 
fastqc *

# Conocer el número de secuencias totales.
zcat SRR21103588_1.fast.gz | echo $((`wc -l`/4))

#1 Indexar el genoma de referencia#
bwa index reference.fasta ;

# Ejecutar #2# y #3# juntos. Generar 1 archivo "bam.bai" y 1 archivo "bam".
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

#4# Extraer genoma consenso. Se generan archivos ".fa" y "qual.txt"
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

# Otorgar permisos a los archivo de trabajo
chmod 777 (nombre de archivo)

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
# Código 5: Observar los resultados en Artemis 
```r
# Es necesario contar con archivos ".fa" y ".gff"
# Primero: crear un ambiente para artemis
conda activate art

# Segundo: entrar a la interfaz colocando "art" en la terminal.

# Tercero: seleccionar el archivo ".gff"
File --> open file manager --> cargar achivo con extensión ".gff".
```
# Código 6: Programas para ensamblaje (Nanopore)
```r
# Se recomienda crear un ambiente nuevo dentro de conda, para cada programa instalado.

# Ver calidad de datos con nanoplot
## Instalar Nanoplot
conda install -c conda-forge -c bioconda nanoplot

# Filtrado de datos con nanofilt
## Insatalar Nanofilt
conda install -c bioconda nanofilt

# Ensamblaje de novo con flye
## Insatalar Flye
conda install -c bioconda flye

# Polishing de datos con minimap2 y racon.
## Instalar Minimap2
conda install -c bioconda minimap2

## Instalar Racon
conda install -c bioconda racon

# Requerimientos de MEDAKA (Pyabpoa, bcftools, samtools (v1.11), minimap2)
## Instalación
pip install pyabpoa
sudo apt install bcftools
conda install -c bioconda samtools==1.11

# Anotación de secuencias consenso con medaka
## Instalar medaka
conda install -c conda-forge –c bioconda medaka
```
# Código 7: Ensamblaje Nanopore (pipeline)
```r
# Descargar los códigos de acceso (códigos: SRR17110067 y SRR17110070)
# Crear directorio "sra_files", generar archivos SRA, mover archivos SRA al ambiente de trabajo, generar archivos fastq, comprimir archivos fastqc y mover los archivos SRA al directorio sra_files.
mkdir sra_files ;
prefetch --max-size 50G --option-file accessions.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra 
gzip *.fastq ;
mkdir sra_files ;
mv *.sra sra_files/ ;

# Inspeccionar las longitudes de los reads
zcat SRR17110067.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20
zcat SRR17110070.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20

# Ver calidad de datos: nanoplot
NanoPlot -t 2 -o SRR17110067_QC --fastq SRR17110067.fastq.gz
NanoPlot -t 2 -o SRR17110070_QC --fastq SRR17110070.fastq.gz

# Filtrar datos: nanofilt
gunzip -c SRR17110067.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110067.trim.fastq.gz ;
gunzip -c SRR17110070.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110070.trim.fastq.gz ;
ls -lh ;

# Ensamblado de lecturas: flye
flye -o SRR17110067.genoma --nano-raw SRR17110067.trim.fastq.gz --threads 4 ;
flye -o SRR17110070.genoma --nano-raw SRR17110070.trim.fastq.gz --threads 4 ;
ls -lh ;

# Polishing: minimap2 + racon
minimap2 -x ava-ont -t 4 SRR17110067.genoma/assembly.fasta SRR17110067.trim.fastq.gz > overlaps1.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps1.paf SRR17110067.genoma/assembly.fasta > SRR17110067.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.genoma/assembly.fasta SRR17110070.trim.fastq.gz > overlaps2.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps2.paf SRR17110070.genoma/assembly.fasta > SRR17110070.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110067.racon1.fasta SRR17110067.trim.fastq.gz > overlaps3.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps3.paf SRR17110067.racon1.fasta > SRR17110067.racon2.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.racon1.fasta SRR17110070.trim.fastq.gz > overlaps4.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps4.paf SRR17110070.racon1.fasta > SRR17110070.racon2.fasta ;

# Anotación de seciencias consenso: medaka
medaka_consensus -i SRR17110070.trim.fastq.gz -d SRR17110070.racon2.fasta -o medaka_SRR17110070 -t 4 ;
medaka_consensus -i SRR17110067.trim.fastq.gz -d SRR17110067.racon2.fasta -o medaka_SRR17110067 -t 4 ;

# Evaluar calidad del emsamblado: QUAST
quast.py -o quast_results -m 0 consensus.fasta
```

# Código 8: BLAST
```r
# Instalacion a traves de CONDA
conda install bioconda::blast
or
conda install -c conda-forge -c bioconda -c defaults blast

# Pasos para descargar Datasets de secuecnias de ADN y proteinas, vinculados a Factores de virulencia.
## http://www.mgc.ac.cn/VFs/
## Default webpage accessible to all users worldwide
## Download
## DNA sequences of full dataset
## Protein sequences of full dataset

# Crearemos formatos FASTA de los archivos descargados
gzip -d VFDB_setB_nt.fas.gz 
gzip -d VFDB_setB_pro.fas.gz

# Corrida en BLAST
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl ;
blastn -db VFDB_setB_nt.fas -query GCA_001183825.1.fasta -perc_identity 90 -outfmt 6 -num_threads 4 > blast.csv ;
head blast.csv ;
cat blast.csv ;

# Colocar los Headers en el archivo ".csv" generado
sed '1i query.acc.ver subject.acc.ver perc.identity alignment.length mismatches gap.opens q.start q.end s.start s.end evalue bit.score' blast.csv | tr " " "\t" > blast.2.csv

# Instalar R
conda install -c conda-forge -c bioconda -c defaults r-base

# Analizar los datos en R
## Leer la data obtedida de blast
data <- read.csv("blast.2.csv", sep="\t", header=TRUE)

## Para conocer el número de filas y columnas de la tabla resultante
dim(data)

## Conocer las filas asignadas a una columna determinada
length(data$subject.acc.ver)

# Conocer el número de elementos únicos de esa columna
length(unique(data$subject.acc.ver))
length(unique(data$query.acc.ver))

## Conocer estadísticos básicos en un solo paso
summary(data$query.acc.ver)
summary(data$alignment.length)
summary(data$perc.identity)

## Obtener un boxplot de los porcentajes de identidad
boxplot(data$perc.identity)
boxplot(data$perc.identity, xlab="genoma", ylab="% identidad")

## Obtener solo los encabezados de cada columna
data.frame(names(data))

## Obtenet un plot de longitud de alineamiento vs %identidad
plot(data$alignment.length, data$perc.identity, xlab="length", ylab="% identity", main="BLASTn VFDB vs Chlamydia", pch=16, col="blue", cex=2)

## Generar un archivo ".bed" en Rstudio
write.table(seq, "extract.txt", sep="\t", row.names = F, col.names =F, quot=F)

# Instalar bedtools desde bioconda para extraer las regiones "blasteadas"
 conda install bioconda::bedtools

bedtools getfasta -fi  GCA_001183825.1.fasta -bed extract.txt -fo virulence.fasta

# Traducir las secuencias con virtual ribosome
https://services.healthtech.dtu.dk/services/VirtualRibosome-2.0/
```
# Código 9: ORTHO-ANI: Para analizar distancias
```
# Instalar NCBI-DATASETS
conda create -n ncbi_datasets
conda activate ncbi_datasets

# Descargar los genomas con la lista sugerida. Emplear el codigo "command_ncbidatasets.sh" disponible en https://github.com/Vjimenez-vasquez/NCBI-DATASETS
./command_ncbidatasets.sh accessions.txt

# guardar los genomas en una nueva carpeta y comprimirla

# Ingresar a la pagina de OAT para correr el algoritmo ORTHO-ANI
https://www.ezbiocloud.net/tools/orthoani

# Descargar BLAST+
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
descargar "ncbi-blast-2.16.0+-win64.exe" para instarlar en el sistema
o
descargar "ncbi-blast-2.16.0+-x64-win64.tar.gz" para descargar los ejecutables directamente

# 10.7: Analizar
https://help.ezbiocloud.net/orthoani-genomic-similarity/
https://pypi.org/project/orthoani/
```
# Código 10: Pangenome analysis: Encontrar Homologías
```
conda install -c conda-forge -c defaults -c bioconda roary
conda install -c conda-forge -c defaults -c bioconda snp-sites
conda install -c conda-forge -c defaults -c bioconda raxml
conda install -c conda-forge -c defaults -c bioconda figtree

# Annotation (PROKKA)
conda activate prokka_env
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fasta
do
prefix=$(basename $r1 .fasta)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Bacteria ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;
conda deactivate ;
cp */*.ffn ffn/ ; 
ls ;

# Inferring clusters, core genes and accesory genes (ROARY): Alineaciones
# https://github.com/sanger-pathogens/Roary #
# roary -p 4 -f roary_output -g 200000 -z -r -e -n -v -cd 80 -i 90 annotation/*.gff ; #
roary -p 4 -f roary_output -g 200000 -r -e -n -v -cd 80 -i 90 annotation/*.gff ;
cp roary_output/core_gene_alignment.aln . ;
ls -lh ; 

# SNPs alignment (SNP-SITES)
## Con gaps
snp-sites -m -o snp1.phy core_gene_alignment.aln ;
## Sin gaps
snp-sites -m -c -o snp2.phy core_gene_alignment.aln ; 
ls -lh ;

# Phylogeny (RAXML): Creo un árbol filogenético
raxmlHPC-PTHREADS -p 1254512 -m GTRCAT -s snp2.phy -n nwk -# 20 -T 4 ;
mv RAxML_bestTree.nwk raw_tree.nwk ;
rm RAxML_* ;
mkdir phylogeny ;
mv snp1.phy snp2.phy snp2.phy.reduced raw_tree.nwk core_gene_alignment.aln phylogeny/ ;

# Cargar el programa "pangenome_command_2.R" en R o R-Studio
# Ingresar la ruta correcta en cada caso (donde se encuentran los archivos "gene_presence_absence.csv" y "metadata_1.tsv")
setwd("")
dir()

# Ingresar la ruta correcta hasta donde se encuentra el archivo pangenome_command_2.R
source("../../pangenome_command_2.R")

pangenome <- pres_abs(metadata = "metadata_1.tsv", roary_output = "gene_presence_absence.csv", last_column = "3", output = "out_5.tsv")
head(pangenome)
class(pangenome)
pangenome[,1:10]

# visualizacion en microreact: visualizar el árbol y la metada.
https://microreact.org/
cargar el arbol enraizado (formato .nwk) y la metadata final (out_5.tsv)
```
