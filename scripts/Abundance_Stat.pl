#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my ($list,$outdir,$prefix,$mock,$control,$help);
GetOptions(
		"l:s" => \$list,
		"o:s" => \$outdir,
		"p:s" => \$prefix,

		"m:s" => \$mock,
		"c:s" => \$control,

		"h|help:s" => \$help,
		);
sub usage{
	print STDERR "\e[;33;1m
	DESCRIPTION
		2b微生物丰度计算程序
	USAGE
		perl $0
	PARAMETERS
		-l <file> list
		          eg:sample<tab>qualitative/quantitative_file
		-o <dir>  outdir file
		-p <str>  output prefix
	OPTIONS
		-m <str>  mock sample (separated by commas)
		-c <str>  negative control sample (separated by commas)
		-h|help   help
	AUTHOR:  ZRC 2020.09.14
	\e[0m\n";
}

if(defined($help)){
	&usage;
	exit 0;
}

unless($list && $outdir && $prefix){
	&usage;
	print STDERR "para -l -o  or -p error.\n";
	exit 1;
}

&CheckDir("$outdir");
#记录mock样品
my (%hash_mock,%hash_control);
if(defined($mock)){
	my @mock=split /,/,$mock;
	for(@mock){
		$hash_mock{$_}++;
	}
}
#记录control样品
if(defined($control)){
	my @control=split /,/,$control;
	for(@control){
		$hash_control{$_}++;
	}
}

#读取定性/定量计算结果文件
my (%hash_specie,%hash_all,@sample_sort,$classify_col,$head);
#循环样品
open LI,"$list" or die "cannot open $list\n";
while(<LI>){
	next if (/^#/ || /^$/);#去除注释行和空行
	chomp;
	my ($sample,$path)=split /\t/;
	$path=abs_path($path);
	unless(-e $path){
		print STDERR "warning: $sample $path not exist, cannot be calculate Abundance.\n";
		next;
	}
	push @sample_sort,$sample;#记录样品顺序
	open IN,"$path" or die "cannot open $path\n";
	while(<IN>){
		chomp;
		my @tmp=split /\t/;
		if(/^#Kingdom/){
			for my $i(0..$#tmp){#确定分类列
				if($tmp[$i] eq "Theoretical_Tag_Num"){
					$classify_col=$i-1;
					$head=join("\t",@tmp[0..$classify_col]);
					last;
				}
			}
		}
		next if(/^#/);#跳过注释行
		my $id=join("\t",@tmp[0..$classify_col]);
		$hash_specie{$id}{$sample}=$tmp[-4];#记录Sequenced_Reads_Num/Theoretical_Tag_Num值
		$hash_all{$sample}+=$tmp[-4];#记录总数
	}
	close IN;
}
close LI;

#输出所有样品丰度计算结果
open OU,">$outdir/$prefix.all.xls" or die "cannot open $outdir/$prefix.all.xls\n";
print OU "$head\t",join("\t",@sample_sort),"\n";#表头
for my $id(sort {$a cmp $b} keys %hash_specie){#循环检测到的物种
	my $judge=0;
	my $print=$id;
	for my $sample(@sample_sort){#循环样品
		if(exists $hash_specie{$id}{$sample}){
			my $percent=$hash_specie{$id}{$sample}/$hash_all{$sample};
			if($percent==0){
				$print .="\t0";
			}else{
				$print .="\t$percent";
				$judge++;
			}
		}else{
			$print .="\t0";
		}
	}
	print OU "$print\n" if ($judge!=0);
}
close OU;

#输出 删除mock和阴性对照样品，以及阴性对照检测出来的菌 的结果
open OU,">$outdir/$prefix.filter.xls" or die "cannot open $outdir/$prefix.filter.xls\n";
#表头处理
print OU "$head";
for(@sample_sort){
	next if(exists $hash_mock{$_} || exists $hash_control{$_});#过滤掉mock和阴性对照样品
	print OU "\t$_";
}
print OU "\n";
for my $id(sort {$a cmp $b} keys %hash_specie){#循环检测到的物种
	my $judge=0;#整行判断
	my $judge_control=0;#阴性对照判断
	my $print=$id;
	for my $sample(@sample_sort){
		next if(exists $hash_mock{$sample});#过滤掉mock样品
		next if(exists $hash_control{$sample});#过滤掉阴性对照样品
#		if(exists $hash_control{$sample}){#阴性对照样品处理
#			$judge_control++ if(exists $hash_specie{$id}{$sample} && $hash_specie{$id}{$sample}!=0);
#			next;
#		}
		if(exists $hash_specie{$id}{$sample}){
			my $percent=0;
			$percent=$hash_specie{$id}{$sample}/$hash_all{$sample} if($hash_all{$sample}!=0);
			if($percent==0){
				$print .="\t0";
			}else{
				$print .="\t$percent";
				$judge++;
			}
		}else{
			$print .="\t0";
		}
	}
	print OU "$print\n" if ($judge!=0 && $judge_control==0);
}

sub CheckDir{#创建目录
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}


























