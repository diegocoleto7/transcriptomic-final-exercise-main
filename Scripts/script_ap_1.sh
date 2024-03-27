#Moving resources to a new dir
mkdir -p Apartado1/res

mv Apartado1/input/Homo_sapiens.GRCh38.109.chr21.gtf Apartado1/input/Homo_sapiens.GRCh38.dna.chromosome.21.fa Apartado1/res/

#FastQC for files

echo -e "\nGenerating FastQC.........."

mkdir -p Apartado1/out/fastQC

for sample in Apartado1/input/*
do
	fastqc "$sample" -o Apartado1/out/fastQC
done

#Indexing genome

echo -e "\nIndexing the genome.........."

hisat2-build --seed 123 -p 10 Apartado1/res/Homo_sapiens.GRCh38.dna.chromosome.21.fa Apartado1/res/Homo_sapiens.GRCh38.dna.chromosome.21

#Alignment of reads

echo -e "\nAligning the reads.........."

mkdir -p Apartado1/out/hisat2/{SRR479052,SRR479054}

hisat2 --new-summary --summary-file Apartado1/out/hisat2/SRR479052/SRR479052.hisat2.summary  --seed 123 --phred33 -p 10 -k 1 -x Apartado1/res/Homo_sapiens.GRCh38.dna.chromosome.21 -1 Apartado1/input/SRR479052.chr21_1.fastq -2 Apartado1/input/SRR479052.chr21_2.fastq -S Apartado1/out/hisat2/SRR479052/SRR479052.sam

hisat2 --new-summary --summary-file Apartado1/out/hisat2/SRR479054/SRR479054.hisat2.summary  --seed 123 --phred33 -p 10 -k 1 -x Apartado1/res/Homo_sapiens.GRCh38.dna.chromosome.21 -1 Apartado1/input/SRR479054.chr21_1.fastq -2 Apartado1/input/SRR479054.chr21_2.fastq -S Apartado1/out/hisat2/SRR479054/SRR479054.sam

#SAMtools:view, sort and index

echo -e "\n SAMtools: view, sort and index.........."

samtools view -bS Apartado1/out/hisat2/SRR479052/SRR479052.sam > Apartado1/out/hisat2/SRR479052/SRR479052.bam
samtools view -bS Apartado1/out/hisat2/SRR479054/SRR479054.sam > Apartado1/out/hisat2/SRR479054/SRR479054.bam

samtools sort Apartado1/out/hisat2/SRR479052/SRR479052.bam -o Apartado1/out/hisat2/SRR479052/SRR479052.sorted.bam
samtools sort Apartado1/out/hisat2/SRR479054/SRR479054.bam -o Apartado1/out/hisat2/SRR479054/SRR479054.sorted.bam

samtools index Apartado1/out/hisat2/SRR479052/SRR479052.sorted.bam
samtools index Apartado1/out/hisat2/SRR479054/SRR479054.sorted.bam

#HTseq-count to count reads

echo -e "\nCounting reads.........."

mkdir -p Apartado1/out/htseq

htseq-count --format=bam --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name Apartado1/out/hisat2/SRR479054/SRR479054.sorted.bam Apartado1/res/Homo_sapiens.GRCh38.109.chr21.gtf > Apartado1/out/htseq/SRR479054.htseq

htseq-count --format=bam --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id --additional-attr=gene_name Apartado1/out/hisat2/SRR479052/SRR479052.sorted.bam Apartado1/res/Homo_sapiens.GRCh38.109.chr21.gtf > Apartado1/out/htseq/SRR479052.htseq

