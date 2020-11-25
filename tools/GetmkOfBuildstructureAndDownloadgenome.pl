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

chomp(my $total=`gzip -dc  $input_list|wc -l`);#总行数
print STDOUT "total line: $total\n";

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
		print OU "ALL:$last_genome\n";
		print OU "\t\@echo \"Congratulations: The entire genome has been downloaded\" >&2\n\n";
		print OU "\$(outdir)/genome_ref/$tmp[-1]\_dustmasked.fna.gz:\n";
		print OU "\t\@mkdir -p \$(outdir)/database/ \$(outdir)/genome_ref/\n";
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
			print OU "\t\@echo \"complete $i%\" >&2\n" unless(exists $hash{$i});#仅没输出过的才会输出日志
			$hash{$i}++;
			last;
		}
	}
	print OU "\n";
	$last_genome=$now_genome;
}
close IN;

