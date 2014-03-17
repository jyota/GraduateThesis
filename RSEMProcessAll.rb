require 'fileutils'

if ARGV.length==0 then 
	puts "Arguments needed. Example:\nruby RSEMProcessAll.rb /media/Thesis/MapSplice-v2.1.3\n(directory is that containing folders w/BAM files\n"
else
	puts "Processing..."
	listDirs = Dir.entries(ARGV[0])
	listDirs = listDirs.grep(/ERR/)
	fullFile = Array.new
    listDirs.each{ |fn| fullFile.push(ARGV[0] + "/" + fn + "/transcriptome_alignments_filtered.bam")}
    fullFile.each_with_index{ |fn, idx|
    	FileUtils.copy(fn,"./" + listDirs[idx] + ".bam")
    	system "./rsem-calculate-expression --paired-end --bam --no-bam-output --estimate-rspd " + listDirs[idx] + ".bam hg19 " + listDirs[idx]
    }
end
