#!/usr/bin/env perl
# Authors: Zheng Sun, Rongchao Zhang, Shi Huang
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="2020.12.16";


select STDOUT;$|=1;# cache cleaning

my $remove_redundant ||="no";# 基因组内部是否去冗余 yes or no, default value is "no"

my($list,$site,$type,$outdir,$enzyme_file,$help);
GetOptions(
		"l:s"  => \$list,
		"s:i"  => \$site,
		"t:s"  => \$type,
		"o:s"  => \$outdir,

		"e:s"  => \$enzyme_file,#酶切结果文件，或库文件
		"r:s"  => \$remove_redundant, #基因组内部是否去冗余
		"h|help:s" => \$help,
		);

sub usage{# help
print STDERR "\e[;33;1m
DESCRIPTION
  It constructs the taxa-specific 2b-RAD reference genome database from a whole-genome reference database.
USAGE
  perl $0
PARAMETERS
  -l <file> genome classification list (the line which begins with # will be ignored)
            eg:GCFid<tab>kingdom<tab>phylum<tab>class<tab>order<tab>family<tab>genus<tab>species<tab>strain(<tab>genome_path)
  -e <file> enzyme file or database file
  -s <int>  2b restriction enzymes (sites).
            [1]CspCI  [9]BplI
            [2]AloI   [10]FalI
            [3]BsaXI  [11]Bsp24I
            [4]BaeI   [12]HaeIV
            [5]BcgI   [13]CjePI
            [6]CjeI   [14]Hin4I
            [7]PpiI   [15]AlfI
            [8]PsrI   [16]BslFI
  -t <str>  The database level. One or more taxonomy level of the 2b-RAD reference database can be specified: kingdom,phylum,class,order,family,genus,species,strain. Use 'all' for any levels. (comma separated).
  -o <dir>  outdir (if not exists,it will be created)
OPTION
  -r <str>  whether to delete redundant tags within the genome (yes or no) [default: $remove_redundant]
  -h|help  print this help
Author:  $author
Last update:  $time\e[0m\n";
}

if(defined($help)){
	&usage;
	exit 0;
}

unless($list && $enzyme_file && $site && $type && $outdir){
	&usage;
	print STDERR "para -l -e -s -t or -o error.\n";
	exit 1;
}

#转化为绝对路径
$list=abs_path($list);
$outdir=abs_path($outdir);
$enzyme_file=abs_path($enzyme_file);

#check the parameter -r: using default value "no"
unless($remove_redundant eq "yes" || $remove_redundant eq "no"){
	&usage;
	print STDERR "-r parameter error: $remove_redundant\n";
	exit 1;
}


#所有分类水平
my %hs_type_database=(
			'kingdom' => '1',
			'phylum'  => '2',
			'class'   => '3',
			'order'   => '4',
			'family'  => '5',
			'genus'   => '6',
			'species' => '7',
			'strain'  => '8',
			);
# check the parameter -t: specify the taxonomic level of 2b-RAD genome database
my %hs_type;
if($type eq "all"){
	%hs_type=%hs_type_database;
}else{
	my @tmp=split /,/,$type;
	for my $i(@tmp){
		if(exists $hs_type_database{$i}){
			$hs_type{$i}=$hs_type_database{$i};
		}else{
			&usage;
			print STDERR "-t parameter error: cannot find '$i'\n";
			exit 1;
		}
	}
}

# check the parameter -s:
my (@site,$enzyme);
if( 1 == $site ){#CspCI
	@site = (
		'[AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}',
		'[AGCT]{10}CCAC[AGCT]{5}TTG[AGCT]{11}',
		);
	$enzyme="CspCI";
}elsif( 2 == $site ){#AloI
	@site = (
		'[AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}',
		'[AGCT]{7}GGA[AGCT]{6}GTTC[AGCT]{7}',
		);
	$enzyme="AloI";
}elsif( 3 == $site ){#BsaXI
    @site = (
            '[AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}',
            '[AGCT]{7}GGAG[AGCT]{5}GT[AGCT]{9}',
            );
	$enzyme="BsaXI";
}elsif( 4 == $site ){#BaeI
    @site = (
            '[AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}',
            '[AGCT]{7}G[AG]TAC[AGCT]{4}GT[AGCT]{10}',
            );
	$enzyme="BaeI";
}elsif( 5 == $site ){#BcgI
    @site = (
            '[AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}',
            '[AGCT]{10}GCA[AGCT]{6}TCG[AGCT]{10}',
            );
	$enzyme="BcgI";
}elsif( 6 == $site ){#CjeI
    @site = (
            '[AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}',
            '[AGCT]{9}AC[AGCT]{6}TGG[AGCT]{8}',
            );
	$enzyme="CjeI";
}elsif( 7 == $site ){#PpiI
    @site = (
            '[AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}',
            '[AGCT]{8}GAG[AGCT]{5}GTTC[AGCT]{7}',
            );
	$enzyme="PpiI";
}elsif( 8 == $site ){#PsrI
    @site = (
            '[AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}',
            '[AGCT]{7}GTA[AGCT]{6}GTTC[AGCT]{7}',
            );
	$enzyme="PsrI";
}elsif( 9 == $site ){#BplI
    @site = (
            '[AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}', #palindromes
            );
	$enzyme="BplI";
}elsif( 10 == $site ){#FalI
    @site = (
            '[AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}', #palindromes
            );
	$enzyme="FalI";
}elsif( 11 == $site ){#Bsp24I
    @site = (
            '[AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}',
            '[AGCT]{7}CCA[AGCT]{6}GTC[AGCT]{8}',
            );
	$enzyme="Bsp24I";
}elsif( 12 == $site ){#HaeIV
    @site = (
            '[AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}',
            '[AGCT]{9}GA[CT][AGCT]{5}[AG]TC[AGCT]{7}',
            );
	$enzyme="HaeIV";
}elsif( 13 == $site ){#CjePI
    @site = (
            '[AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}',
            '[AGCT]{8}GA[AGCT]{7}TGG[AGCT]{7}',
            );
	$enzyme="CjePI";
}elsif( 14 == $site ){#Hin4I
    @site = (
            '[AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}',
            '[AGCT]{8}GA[CTG][AGCT]{5}[AG]TC[AGCT]{8}',
            );
	$enzyme="Hin4I";
}elsif( 15 == $site ){#AlfI
    @site = (
            '[AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}', #palindromes
            );
	$enzyme="AlfI";
}elsif( 16 == $site ){#BslFI ??some question?? single enzyme
    @site = (
            '[AGCT]{6}GGGAC[AGCT]{14}',
            '[AGCT]{14}GTCCC[AGCT]{6}',
            );
	$enzyme="BslFI";
}else{
	&usage;
	print STDERR "The parameter -s is wrong\n";
	exit 1;
}

#提供酶切文件，检查文件是否存在
if(defined($enzyme_file)){
	unless(-e $enzyme_file){
		print STDERR "[ERROR] $enzyme_file does not exist,please check.\n";
		exit 1;
	}
}


#统计总的基因组个数
my $genome_total_num=0;
#my (%hash_gcf2class,%hash_gcf_rank);
my %hash_gcf2class;
if($list=~/\.gz/){
	open LI,"gzip -dc $list|" or die "cannot open $list\n";
}else{
	open LI,"$list" or die "cannot open $list\n";
}
while(<LI>){
	next if(/^#/ || /^$/);# remove blank lines or lines starting with #
	chomp;
	my @tmp=split /\t/;
	$genome_total_num++;#总基因组个数
	$hash_gcf2class{$tmp[0]}=$_;#gcf对应的分类信息
#	$hash_gcf_rank{$tmp[0]}=$genome_total_num;#记录基因组在列表中的排序，便于打印日志
	if(defined($enzyme_file)){#提供酶切文件
		;
	}else{#不提供酶切文件
		$tmp[-1]=abs_path($tmp[-1]);
		unless(-e $tmp[-1]){#check the availability of a genome fasta file
			print STDERR "[ERROR] $tmp[-1] does not exist,please check your genome file\n";
			exit 1;
		}
	}
}
close LI;

if($genome_total_num==0){
	print STDERR "[warning] There is no genome in the List file.\n";
	exit 0;
}


&CheckDir("$outdir");# create the output directory
#&CheckDir("$outdir/database");

print STDOUT "###($enzyme) Record the taxonomies of each 2b-RAD tag and identification of taxa-specifc 2b-RAD tags -- start, ",`date`;#STDOUT
for my $level(sort {$hs_type{$a}<=>$hs_type{$b}} keys %hs_type){ #iterate all taxonomic levels of 2b-RAD database
	print STDOUT "###($level) Record the taxonomies of each 2b-RAD tag -- start, ",`date`;#STDOUT
	my (%hash_ingenome,%hash,%complete,%hash_gcf_rank);
	my %hash_seq;#记录基因组酶切所有标签
	$/=">";
	if($enzyme_file=~/\.gz$/){#打开酶切文件
		open IN,"gzip -dc $enzyme_file|" or die "cannot open $enzyme_file\n";
	}else{
		open IN,"$enzyme_file" or die "cannot open $enzyme_file\n";
	}
	<IN>;
	while(<IN>){
		chomp;
		my @tmp=split /\n/;
		#GCF号|基因组内部标签排序|scaffoldid|startpos|正反向酶切|是否为指定水平下unique&&noredundancy标签
		my ($gcfid,$ingenome_tag_num,$scaid,$start,$chain,$unique)=split /\|/,$tmp[0];
		my $tag=$tmp[1];
		next unless(exists $hash_gcf2class{$gcfid});#gcfid不在列表中则跳过
		$hash_gcf_rank{$gcfid}++;
		$hash_seq{$gcfid}{$ingenome_tag_num}=join("\n",@tmp[0..1]);#记录列表中，基因组酶切的所有标签

		my @a=split /\t/,$hash_gcf2class{$gcfid};#分类
		my $class=join("\t",@a[1..$hs_type{$level}]);#concatenate the full taxonomic annotation
		if(exists $hash{$tag}){#判断在哈希中是否存在
			$hash{$tag}{$class}++;#记录标签分类信息
			$hash_ingenome{$gcfid}{$tag}++ if($remove_redundant eq "yes");#如果需要去除基因组内部冗余，则记录标签在基因组内部是否冗余
		}else{#反向互补处理
			$tag=~tr/ATCG/TAGC/;
			$tag=reverse($tag);
			$hash{$tag}{$class}++;#记录标签分类信息
			$hash_ingenome{$gcfid}{$tag}++ if($remove_redundant eq "yes");#如果需要去除基因组内部冗余，则记录标签在基因组内部是否冗余
		}
		for (my $i=100;$i>0;$i=$i-1){ #每完成1%则输出进度
			if((keys %hash_gcf_rank)/$genome_total_num*100>=$i){
				print STDOUT "$i% " unless(exists $complete{$i});#仅没输出过的才会输出日志
				$complete{$i}++;
				last;
			}
		}
	}
	close IN;
	$/="\n";
	undef %hash_gcf_rank;
	print STDOUT "\n###($level) Record the taxonomies of each 2b-RAD tag -- complete, ",`date`;# STDOUT
	
	print STDOUT "###($level) Identification of taxa-specifc 2b-RAD tags -- start, ",`date`;# STDOUT
	undef %complete;#完成进度清空
	my (%hash_genome_tag_num,%hash_genome_unique_tag_num);
	my $complete=0;
	if($list=~/\.gz/){
		open LI,"gzip -dc $list|" or die "cannot open $list\n";
	}else{
		open LI,"$list" or die "cannot open $list\n";
	}
	open OU,"|gzip > $outdir/$enzyme.$level.fa.gz" or die "cannot open $outdir/$enzyme.$level.fa.gz\n";
	while(<LI>){
		next if(/^#/ || /^#/);# remove blank lines or lines starting with #
		my $line=$_;
		chomp($line);
		my $gcfid=(split /\t/,$line)[0];
		next unless(exists $hash_seq{$gcfid});
		for my $i(sort {$a<=>$b} keys %{$hash_seq{$gcfid}}){#循环标签
			my @tmp=split /\n/,$hash_seq{$gcfid}{$i};#id\nseq
			#GCF号|基因组内部标签排序|scaffoldid|startpos|正反向酶切|是否为指定水平下unique&&noredundancy标签
			my ($gcfid,$ingenome_tag_num,$scaid,$start,$chain,$unique)=split /\|/,$tmp[0];
			my $tag=$tmp[1];
			$gcfid=~s/^>//;
			unless(exists $hash{$tag}){#如果不存在，则进行反向互补
				$tag=~tr/ATCG/TAGC/;
				$tag=reverse($tag);
				#反向互补后，改变链的方向
				if($chain==0){
					$chain=1;
				}elsif($chain==1){
					$chain=0;
				}
			}
			if(keys %{$hash{$tag}}==1){#指定水平下为unique
				if($remove_redundant eq "yes"){#基因组内部需要去冗余
					if($hash_ingenome{$gcfid}{$tag}==1){#在基因组内部只出现过一次（noredundancy）
						$unique=1;
						$hash_genome_unique_tag_num{$gcfid}++;#基因组电子酶切unique标签数
					}else{
						$unique=0;
					}
				}else{#基因组内部不需要去冗余
					$unique=1;
					$hash_genome_unique_tag_num{$gcfid}++;#基因组电子酶切unique标签数
				}
			}else{
				$unique=0;
			}
			print OU ">$gcfid|$ingenome_tag_num|$scaid|$start|$chain|$unique\n$tag\n" if($unique==1);
		}
		$complete++;
		for (my $i=100;$i>0;$i=$i-1){ #每完成1%则输出进度
			if($complete/$genome_total_num*100>=$i){
				print STDOUT "$i% " unless(exists $complete{$i});#仅没输出过的才会输出日志
				$complete{$i}++;
				last;
			}
		}

	}
	close LI;
	close OU;
	#统计输出
	open STAT,"> $outdir/$enzyme.$level.stat.xls" or die "cannot open $outdir/$enzyme.$level.stat.xls\n";
	print STAT "#Unique_Name\tAll_Tag_Num\tUnique_Tag_Num\n";
	if($list=~/\.gz/){
		open LI,"gzip -dc $list|" or die "cannot open $list\n";
	}else{
		open LI,"$list" or die "cannot open $list\n";
	}
	while(<LI>){
		next if(/^#/ || /^$/);# remove blank lines or lines starting with #
		my $line=$_;
		chomp($line);
		my @tmp=split /\t/,$line;
		print STAT "$tmp[0]";
		if(exists $hash_seq{$tmp[0]}){#酶切标签数
			my $genome_tag_num=keys %{$hash_seq{$tmp[0]}};
			print STAT "\t$genome_tag_num";
		}else{
			print STAT "\t0";
		}
		if(exists $hash_genome_unique_tag_num{$tmp[0]}){#unique标签数
			print STAT "\t$hash_genome_unique_tag_num{$tmp[0]}\n";
		}else{
			print STAT "\t0\n";
		}
#		print STAT "$tmp[0]\t$hash_genome_tag_num{$tmp[0]}\t$hash_genome_unique_tag_num{$tmp[0]}\n";
	}
	close LI;
	close STAT;
	print STDOUT "\n###($level) Identification of taxa-specifc 2b-RAD tags -- complete, ",`date`; #STDOUT
}

print STDOUT "###($enzyme) Record the taxonomies of each 2b-RAD tag and identification of taxa-specifc 2b-RAD tags -- complete, ",`date`;#STDOUT




sub CheckDir{# create the directory
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
		}
		return 1;
}
