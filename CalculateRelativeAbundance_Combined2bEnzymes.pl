#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="2020.06.03";

#默认值
my $mark ||="combine";
my $g_score_threshold ||=0;

#Standard output for clearing cache
select STDOUT;$|=1;

my ($list,$site,$outdir);
GetOptions(
		"l:s" => \$list,
		"s:s" => \$site,
		"io:s" => \$outdir,

		"m:s" => \$mark,
		"g:i" => \$g_score_threshold,
		);


sub usage{# helper information
print STDERR "\e[;33;1m
DESCRIPTION
    It computes the relative abundance of taxa identified from each of 2b-RAD samples using a precalcuated taxa-specific 2b-RAD reference database by one or multiple type 2b restriction enzymes.
USAGE
	perl $0
PARAMETERS
	-l  <str> The path of the input filepath list (the line which begins with # will be ignored)
	       eg: sample_name<tab>...
	-s  <int> One or multiple type 2b restriction enzymes (sites). The selected sites should be separated by comma.
	       [1]CspCI  [9]BplI
	       [2]AloI   [10]FalI
	       [3]BsaXI  [11]Bsp24I
	       [4]BaeI   [12]HaeIV
	       [5]BcgI   [13]CjePI
	       [6]CjeI   [14]Hin4I
	       [7]PpiI   [15]AlfI
	       [8]PsrI   [16]BslFI
	       [17]All_Detected_Enzyme
	-io <str> The input and output directory
OPTIONS
	-m  <str> Whether the taxa idenfication or abundance estimation should take into account for the 2b-RAD taxa-specific markers from more than one restriction sites [combine]
	-g  <int> The G-score threshold [$g_score_threshold, it means >=$g_score_threshold] To control the false-positive in the species identification, G score was derived for each species identified within a sample, which is a harmonious mean of read coverage of 2b-RAD markers belongs to a species and number of all possible 2b-RAD markers of this species. Therecommended/default threshold is 5.
AUTHOR:  $author $time\e[0m\n";
}


my %hs_site2enzyme=(#the codes for all restriction enzymes
	'1'  =>  'CspCI',   '2'  =>  'AloI',
	'3'  =>  'BsaXI',   '4'  =>  'BaeI',
	'5'  =>  'BcgI',    '6'  =>  'CjeI',
	'7'  =>  'PpiI',    '8'  =>  'PsrI',
	'9'  =>  'BplI',    '10' =>  'FalI',
	'11' =>  'Bsp24I',  '12' =>  'HaeIV',
	'13' =>  'CjePI',   '14' =>  'Hin4I',
	'15' =>  'AlfI',    '16' =>  'BslFI',
	);

unless($list && $site && $outdir){
	print "1\n";
	&usage;
	exit;
}

# check the availability of the taxa-specific 2b-RAD reference genome database
if($site=~/17/){
	$site="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16";
}
my @site=split /,/,$site;
for $site(@site){
	unless(exists $hs_site2enzyme{$site}){#检测酶切位点是否存在
		&usage;
		print STDERR "Parameter -s is wrong\n";
		exit;
	}
}
=pod
#注释文件检测并读取
my %hs_anno;
unless(-e "$database/species_ID_annotation.txt"){
	&usage;
	print STDERR "cannot find $database/species_ID_annotation.txt\n";
	exit;
}else{
	open AN,"$database/species_ID_annotation.txt" or die "cannot open $database/species_ID_annotation.txt\n";
	while(<AN>){
		next if(/^#/ || /^$/);#去除注释行和空行
		chomp;
		my @tmp=split /\t/;
		$hs_anno{$tmp[2]}="$tmp[1]\t$tmp[3]";
	}
	close AN;
}
=cut

#合并处理
open LI,"$list" or die "cannot open $list\n";
while(<LI>){#循环样品
	next if(/^#/ || /^$/);#去除注释行和空行
	chomp;
	my (%hs_sample,$head);
	my @use_site;#用到的酶 结果
	my $sample_name=(split /\t/)[0];
	my $cnt=0;
	for $site(@site){# iterate all enzymes
		if(-e "$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls"){
			push @use_site,$hs_site2enzyme{$site};
		}else{
			print STDERR "warning: cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls\n";
			next;#跳过没有鉴定的酶 文件
		}
		$cnt++;
		open IN,"$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls";
		while(<IN>){
			chomp;
			my $line=$_;
			next if($line=~/^#/);
			if($line=~/^Kingdom/){
				$head=$line;
				next;
			}
			my @tmp=split /\t/;
			my $class=join("\t",@tmp[0..$#tmp-8]);
			$hs_sample{$class}{-8}+=$tmp[-8];#理论标签数
			$hs_sample{$class}{-7}+=$tmp[-7];#测序得到的标签数
			$hs_sample{$class}{-5}+=$tmp[-5];#测序得到的reads数
			$hs_sample{$class}{-2}+=$tmp[-2];#测到的深度大于1的标签数
		}
		close IN;
	}
	next if($cnt==0);#如果所有的酶都没有结果，那么不继续输出
	open OU,">$outdir/$sample_name/$sample_name.$mark.xls" or die "cannot open $outdir/$sample_name/$sample_name.$mark.xls\n";
	print OU "#@use_site combine\n";
#	print OU "#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecie\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\tSequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\tSequenced_Reads_Num/Sequenced_Tag_Num\tG_Score\ttaxid\tunique_name\n";
	print OU "#$head\n";
	for my $class(keys %hs_sample){
		my @tmp=split /\t/,$class;
		my $Theoretical_Tag_Num=$hs_sample{$class}{-8};
		my $Sequenced_Tag_Num=$hs_sample{$class}{-7};
		my $Sequenced_Tag_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Tag_Num/$Theoretical_Tag_Num*100;
		my $Sequenced_Reads_Num=$hs_sample{$class}{-5};
		my $Sequenced_Reads_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Theoretical_Tag_Num;
		my $Sequenced_Reads_Num2Sequenced_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Sequenced_Tag_Num;
		my $Sequenced_Tag_Num_2=$hs_sample{$class}{-2};
		my $G_Score=sprintf "%.8f",sqrt($Sequenced_Tag_Num*$Sequenced_Reads_Num);
		next if ($G_Score<$g_score_threshold);#过滤gscore阈值
		print OU "$class\t$Theoretical_Tag_Num\t$Sequenced_Tag_Num\t$Sequenced_Tag_Num2Theoretical_Tag_Num%\t";
		print OU "$Sequenced_Reads_Num\t$Sequenced_Reads_Num2Theoretical_Tag_Num\t$Sequenced_Reads_Num2Sequenced_Tag_Num\t";
#		print OU "$G_Score\t$hs_anno{$tmp[6]}\n";
		print OU "$Sequenced_Tag_Num_2\t$G_Score\n";
	}
	close OU;
	undef %hs_sample;
}
close IN;







