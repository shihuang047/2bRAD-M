#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use MyModule;

my $author="Zhangrongchao";
my $time="20200924";

my $g_score_threshold ||=5;#对定性的合并结果，进行分类筛选 gscore阈值
my $GCF_threshold ||=1;#鉴定到某个基因组几个标签以上，该基因组才会被纳入定量建库


my($list,$database,$outdir,$qual_dir,$help);
GetOptions(
		"l:s"   => \$list, #待处理样品列表
		"d:s"   => \$database, #数据库目录
		"o:s"   => \$outdir, #输出目录
		"qualdir:s" => \$qual_dir, #所有样品定性结果总目录

		"gscore:i"   => \$g_score_threshold,#筛选分类
		"gcf:i" => \$GCF_threshold,#筛选分类中的基因组

		"h|help:s" => \$help,
		);

sub usage{#帮助
print STDERR "\e[;33;1m
	DESCRIPTION
		2b微生物根据定性结果筛选定量基因组生成建库所需list
	USAGE
		perl $0
	PARAMETERS
		-l  <file> sample list (the line which begins with # will be ignored)
		           eg: sample(<tab>...<tab>...)
		-d  <dir> database path
		-o  <dir> outdir (if not exists,it will be created)
		-qualdir <dir> dir of qualitative
	OPTIONS
		-gscore <i> G score threshold of classify in qualitative analysis, it decides quantitative database. [$g_score_threshold, it means >$g_score_threshold]
		-gcf    <i> detected tag threshold of GCF in qualitative analysis, it decides quantitative database. [$GCF_threshold, it means >$GCF_threshold]
		-h|help     print help
	AUTHOR:  $author $time\e[0m\n";
}


#参数检测
unless($list && $database && $outdir && $qual_dir){
	&usage;
	print STDERR "para -l -d -o or -qualdir error.\n";
	exit 1;
}

#数据库文件检测
unless(-d "$database/database" && -e "$database/abfh_classify.txt" && -d "$database/genome_ref"){
	print STDERR "incomplete database, please check.\n";
	exit 1;
}

#记录数据库中gcf转化为全部信息
my %gcf2classify_path;
open DB,"$database/abfh_classify.txt" or die "cannot open $database/abfh_classify.txt\n";
while(<DB>){
	chomp;
	next if(/^#unique_name/); #unique_name    kingdom phylum  class   order   family  genus   specie  strain  genome_path
	my @tmp=split /\t/;
	$gcf2classify_path{$tmp[0]}=$_;
}
close DB;


&CheckDir("$outdir");

open IN,"$list" or die "cannot open $list\n";
while(<IN>){
	next if(/^#/ || /^$/);#跳过注释行和空行
	chomp;
	my $sample_name=(split /\t/)[0];
	#样品定量列表
	##单样品重建库list准备
	##记录通过G_score阈值的分类
	my (%hs_pass_Gscore_class,@enzyme_use);
	if(-e "$qual_dir/$sample_name/$sample_name.combine.xls"){
		open XI,"$qual_dir/$sample_name/$sample_name.combine.xls" or die "cannot open $qual_dir/$sample_name/$sample_name.combine.xls\n";
	}else{
		print STDERR "!!!$sample_name does not have $qual_dir/$sample_name/$sample_name.combine.xls, can't do quantitative analysis\n";
		next;
	}
	while(<XI>){
		chomp;
		my @tmp=split /\t/;
		next if(/^#Kingdom/);#跳过表头
		if(/^#/){#记录 合并使用的酶 组合
			my @a=split /\s+/,$tmp[0];
			for my $enzyme(@a){
				$enzyme=~s/^#//;
				next if($enzyme eq "combine");#跳过 combine字段
				push @enzyme_use,$enzyme;
			}
			next;
		}
		my $class=join("\t",@tmp[0..$#tmp-8]);#获取分类信息
		$hs_pass_Gscore_class{$class}++ if($tmp[-1]>$g_score_threshold);#通过gscore阈值的分类
	}
	close XI;
	&CheckDir("$outdir/$sample_name");#建立每个样品的文件夹
	open OU,"|sort|uniq > $outdir/$sample_name/sdb.list" or die "cannot open $outdir/$sample_name/sdb.list\n"; #输出选出的基因组列表，并排序去重
	for my $enzyme(@enzyme_use){
		if(-e "$qual_dir/$sample_name/$sample_name.$enzyme.GCF_detected.xls"){
			open QU,"$qual_dir/$sample_name/$sample_name.$enzyme.GCF_detected.xls" or die "cannot open $qual_dir/$sample_name/$sample_name.$enzyme.GCF_detected.xls\n";
		}else{
			print STDERR "warning: $sample_name does not have $qual_dir/$sample_name/$sample_name.$enzyme.GCF_detected.xls\n";
			next;
		}
		while(<QU>){
			chomp;
			my @tmp=split /\t/;
			my $class=join("\t",@tmp[0..$#tmp-4]);
			if(exists $hs_pass_Gscore_class{$class} && $tmp[-2]>$GCF_threshold){
				my @all_class=split /\t/,$gcf2classify_path{$tmp[-4]};
				print OU join("\t",@all_class[0..8]),"\t$database/genome_ref/$all_class[-1]\n";
			}
		}
		close QU;
	}
	close OU;

}
close IN;









