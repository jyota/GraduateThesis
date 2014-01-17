if [ -e "$1.sra" ]
then
  echo "Processing SRA file."
  ../sra*/bin/fastq-dump ./$1.sra --split-3
  echo "Converting strand one to phred 64."
  java -Xmx512M -jar ~/picard/picard-code/ubu-master/target/ubu*ies.jar fastq-format --phred33to64 --strip --suffix /1 --in $1_1.fastq --out prep_1.fastq
  echo "Converting strand two to phred 64."
java -Xmx512M -jar ~/picard/picard-code/ubu-master/target/ubu*ies.jar fastq-format --phred33to64 --strip --suffix /2 --in $1_2.fastq --out prep_2.fastq
  echo "Removing intermediate files."
  rm $1.sra
  rm $1_1.fastq
  rm $1_2.fastq
  echo "Running MapSplice."
  python mapsplice.py -p 2 --fusion -c ~/192.168.1.133/data/chromFa -x ~/192.168.1.133/data/chromFa/hg19 -1 prep_1.fastq -2 prep_2.fastq -o $1
  echo "Cleanup unnecessary data."
  cd $1
  if [${PWD##*/} == $1 ]
  then
	  if [ -e "alignments.sam" ]
	  then
	    rm *.txt
	    rm -rf logs
	    rm ../prep_1.fastq
	    rm ../prep_2.fastq
	  fi
	  echo 'Add read groups'
	  java -Xmx2G -jar ~/picard/picard-code/dist/AddOrReplaceReadGroups.jar INPUT=~/MapSplice-v2.1.3/$1/alignments.sam OUTPUT=~/MapSplice-v2.1.3/$1/rg_alignments.bam RGSM=$1 RGID=$1 RGPL=Illumina RGLB=TruSeq RGPU=barcode VALIDATION_STRINGENCY=SILENT TMP_DIR=.
	  echo 'Convert back to Phred33'
	  java -Xmx512M -jar ~/picard/picard-code/ubu-master/target/ubu*ies.jar sam-convert --phred64to33 --in rg_alignments.bam --out phred33_alignments.bam
	  echo 'Sort by coordinate'
	  samtools sort phred33_alignments.bam sorted_genome_alignments
	  echo 'index'
	  samtools index sorted_genome_alignments.bam
	  echo 'sort by chromosome then read ID'
	  perl ~/picard/picard-code/ubu-master/src/perl/sort_bam_by_reference_and_name.pl --input sorted_genome_alignments.bam --output sorted_by_chr_read.bam --temp-dir . --samtools=/usr/bin/samtools
	  echo 'Translate to transcriptome coordinates for RSEM'
	  java -Xms3G -Xmx12G -jar ~/picard/picard-code/ubu-master/target/ubu*ies.jar sam-xlate --bed ~/192.168.1.133/data/hg19transcripts.bed --in sorted_by_chr_read.bam --out transcriptome_alignments.bam --order ~/rsem-1.2.4/hg19.idx.fa --xgtags --reverse
	  echo 'Filter indels, large inserts, zero mapping quality'
	  java -Xms3G -Xmx12G -jar ~/picard/picard-code/ubu-master/target/ubu*ies.jar sam-filter --in transcriptome_alignments.bam --out transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1

	  echo 'Cleanup intermediate files.'
	  if [-e "transcriptome_alignments_filtered.bam"]
	  then
	    rm alignments.sam
	    rm transcriptome_alignments.bam
	    rm sorted_genome_alignments.bam
	    rm sorted_genome_alignments.bam.bai
	    rm rg_alignments.bam
	    rm phred33_alignments.bam
	    find . -type d -delete
            echo "Run complete for $1."
	else
	    echo "MapSplice output directory does not seem to exist for $1! Exiting."
  fi
fi

