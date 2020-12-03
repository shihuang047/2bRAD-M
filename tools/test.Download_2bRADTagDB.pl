#!/usr/bin/perl
#Author:zhangrongchao, zhangrongchaoxx@163.com
use strict;
use warnings;

if($#ARGV!=1){
	print STDERR "perl $0 input.path.list.gz output.mk\n";
	exit 1;
}

my $input_list=$ARGV[0];
my $output=$ARGV[1];

unless($input_list=~/\.gz$/){
	print STDERR "$input_list must be gz.\n";
	exit 1;
}

my %hash_database_path=(#16files
	'PsrI'  =>  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25621073/PsrI.specie.gz',
	);
my %hash_abfh_path=(
	'abfh'  =>  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25623785/abfh_classify.txt.gz'
	);
my %hash_example_path=(#2files
	'simulate_50.fa.gz' => 'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25621832/simulate_50.fa.gz',
	);



chomp(my $total=`gzip -dc  $input_list|wc -l`);#总行数
print STDOUT "genome total line: $total\n";

chomp(my $last_genome=`gzip -dc $input_list|tail -l`);#获取最后一个

open OU,">$output" or die "cannot open $output\n";
#print OU "#Download directory is \$(outdir).\" >&2\n";
print OU "#Author: zhangrongchao; Email: zhangrongchaoxx@163.com.\n";
print OU "#Usage: make -f $output Database_path=yourpath(default:./2B-RAD-M-ref_db/).\n\n\n";
print OU "Database_path?=./2B-RAD-M-ref_db/\n";
print OU "outdir=\$(abspath \$(Database_path))\n\n";

my $num=0;
my %hash;
open IN,"gzip -dc $input_list|" or die "cannot open $input_list\n";
while(<IN>){
	chomp;
	$num++;
	my $now_genome=$_;
	#上一个基因组处理
	$last_genome=(split /\//,$last_genome)[-1]; #GCF_000183285.1_ASM18328v1_genomic.fna.gz
	$last_genome=~s/\.fna\.gz$//;#GCF_000183285.1_ASM18328v1_genomic
	$last_genome="\$(outdir)/genome_ref/$last_genome\_dustmasked.fna.gz";
	#当前基因组处理
	my @tmp=split /\//,$now_genome;
	$tmp[-1]=~s/\.fna\.gz$//;#GCF_000007365.1_ASM736v1_genomic
	if($num==1){
#		print OU "ALL:$last_genome\n";
		print OU "ALL:\$(outdir)/example_data/simulate_50.fa.gz\n";
		print OU "\t\@echo \"Congratulations: The all files have been downloaded\" >&2\n\n";

		print OU "\$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna.gz:\n";

		print OU "\t\@echo \"Create Directory\" >&2\n";
		print OU "\t\@mkdir -p \$(outdir)/genome_ref/\n";
		print OU "\t\@mkdir -p \$(outdir)/example_data/\n";
		print OU "\t\@echo \"MSA1002\t\$(outdir)/example_data/MSA1002_R1.fq.gz\" > \$(outdir)/example_data/list_mock\n";

		print OU "\t\@mkdir -p \$(outdir)/database/\n";
		print OU "\t\@mkdir -p \$(outdir)/database/AlfI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/AloI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/BaeI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/BcgI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/BplI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/BsaXI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/BslFI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/Bsp24I\n";
		print OU "\t\@mkdir -p \$(outdir)/database/CjeI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/CjePI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/CspCI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/FalI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/HaeIV\n";
		print OU "\t\@mkdir -p \$(outdir)/database/Hin4I\n";
		print OU "\t\@mkdir -p \$(outdir)/database/PpiI\n";
		print OU "\t\@mkdir -p \$(outdir)/database/PsrI\n";

		print OU "\t\@echo \"Start downloading the genome\" >&2\n";
	}else{
		print OU "\$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna.gz:$last_genome\n";
	}
	print OU "\t\@wget -q -O \$(outdir)/genome_ref/$tmp[-1].fna.gz $now_genome\n";
	print OU "\t\@gunzip \$(outdir)/genome_ref/$tmp[-1].fna.gz\n";
	print OU "\t\@dustmasker -infmt fasta -in \$(outdir)/genome_ref/$tmp[-1].fna -level 20 -outfmt fasta ";
	print OU ">\$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna\n";
	print OU "\t\@sed -i '/^>/! s/[^AGCT]/N/g' \$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna\n";
	print OU "\t\@gzip -f \$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna\n";
	print OU "\t\@rm -f \$(outdir)/genome_ref/$tmp[-1].fna\n";
	for (my $i=100;$i>0;$i=$i-5){
		if($num/$total*100>=$i){
			print OU "\t\@echo \"Genome download completed $i%\" >&2\n" unless(exists $hash{$i});#仅没输出过的才会输出日志
			$hash{$i}++;
			last;
		}
	}
	print OU "\n";
	$last_genome=$now_genome;
}
close IN;

#上一个基因组处理
$last_genome=(split /\//,$last_genome)[-1]; #GCF_000183285.1_ASM18328v1_genomic.fna.gz
$last_genome=~s/\.fna\.gz$//;#GCF_000183285.1_ASM18328v1_genomic
$last_genome="\$(outdir)/genome_ref/$last_genome\_dustmasked.fna.gz";

#database 下载
$num=0;
for my $enzyme(sort {$a cmp $b} keys %hash_database_path){
	$num++;
	if($num==1){
		print OU "\$(outdir)/database/$enzyme/specie.gz:$last_genome\n";
		print OU "\t\@echo \"Start downloading the 2bRAD-M species database\" >&2\n";
	}else{
		print OU "\$(outdir)/database/$enzyme/specie.gz:$last_genome\n";
	}
	print OU "\t\@wget -q -O \$(outdir)/database/$enzyme/specie.gz $hash_database_path{$enzyme}\n\n";
	$last_genome="\$(outdir)/database/$enzyme/specie.gz";
}

#abfh_classify
print OU "\$(outdir)/abfh_classify.txt:$last_genome\n";
print OU "\t\@wget -q -O \$(outdir)/abfh_classify.txt.gz $hash_abfh_path{abfh}\n";
print OU "\t\@gunzip -f \$(outdir)/abfh_classify.txt.gz\n\n";

#example data
$num=0;
for my $example(sort {$a cmp $b} keys %hash_example_path){
	$num++;
	if($num==1){
		print OU "\$(outdir)/example_data/$example:\$(outdir)/abfh_classify.txt\n";
		print OU "\t\@echo \"Start downloading the sample data\" >&2\n";
	}else{
		print OU "\$(outdir)/example_data/$example:$last_genome\n";
	}
	print OU "\t\@wget -q -O \$(outdir)/example_data/$example $hash_example_path{$example}\n\n";
	$last_genome="\$(outdir)/example_data/$example";
}
close OU;





