#!/usr/bin/env perl
# Authors: Zheng Sun, Rongchao Zhang, Shi Huang
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use Parallel::ForkManager;

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="2020.06.03";

my $cpu ||=10;
my $remove_redundant ||="no";

select STDOUT;$|=1;#标准输出清楚缓存

my($list,$site,$type,$outdir,$verbose);
GetOptions(
		"l:s"  => \$list,
		"s:i"  => \$site,
		"t:s"  => \$type,
		"o:s"  => \$outdir,

		"c:i"  => \$cpu,
		"v:s"  => \$verbose,
		"r:s"  => \$remove_redundant,#移除冗余标签 yes or no
		);

sub usage{#帮助
print STDERR "\e[;33;1m
DESCRIPTION
  It constructs the taxa-specific 2b-RAD reference genome database from a whole-genome reference database.
USAGE
  perl $0
PARAMETERS
  -l <str> genome classification list (the line which begins with # will be ignored)
         eg:unique_name<tab>kingdom<tab>phylum<tab>class<tab>order<tab>family<tab>genus<tab>specie<tab>strain<tab>genome_path
  -s <int>  One or multiple type 2b restriction enzymes (sites).
         [1]CspCI  [9]BplI
         [2]AloI   [10]FalI
         [3]BsaXI  [11]Bsp24I
         [4]BaeI   [12]HaeIV
         [5]BcgI   [13]CjePI
         [6]CjeI   [14]Hin4I
         [7]PpiI   [15]AlfI
         [8]PsrI   [16]BslFI
  -t <str> The database level. One or more taxonomy level of the 2b-RAD reference database can be specified: kingdom,phylum,class,order,family,genus,specie,strain. Use 'all' for any levels. (comma separated).
  -o <str> outdir (if not exists,it will be created)
OPTION
  -c <int> cpu [$cpu]
  -v <str> verbose
AUTHOR:  $author $time\e[0m\n";
}

unless($list && $site && $type && $outdir){
	&usage;
	exit;
}

my %hs_site2enzyme=(#酶切位点对应表
	'1'  =>  'CspCI',	'2'  =>  'AloI',
	'3'  =>  'BsaXI',	'4'  =>  'BaeI',
	'5'  =>  'BcgI',	'6'  =>  'CjeI',
	'7'  =>  'PpiI',	'8'  =>  'PsrI',
	'9'  =>  'BplI',	'10' =>  'FalI',
	'11' =>  'Bsp24I',	'12' =>  'HaeIV',
	'13' =>  'CjePI',	'14' =>  'Hin4I',
	'15' =>  'AlfI',	'16' =>  'BslFI',
	);

my %hs_type_database=(
			'kingdom' => '1',
			'phylum'  => '2',
			'class'   => '3',
			'order'   => '4',
			'family'  => '5',
			'genus'   => '6',
			'specie'  => '7',
			'strain'  => '8',
			);
my @head=('#Name','Kingdom','Phylum','Class','Order','Family','Genus','Specie','Strain');

#判断需要分析的类别
my %hs_type;
if($type=~/all/){
	%hs_type=%hs_type_database;
}else{
	my @tmp=split /,/,$type;
	for my $i(@tmp){
		if(exists $hs_type_database{$i}){
			$hs_type{$i}=$hs_type_database{$i};
		}else{
			&usage;
			print STDERR "-t parameter error: cannot find '$i'\n";
			exit;
		}
	}
}
#判断酶切位点是否存在
unless(exists $hs_site2enzyme{$site}){
	&usage;
	print STDERR "-s parameter error: cannot find $site\n";
	exit;
}

#判断-r $remove_redundant参数
unless($remove_redundant eq "yes" || $remove_redundant eq "no"){
	&usage;
	print STDERR "-r parameter error: cannot find $remove_redundant\n";
	exit;
}

&CheckDir($outdir);#创建文件夹
&CheckDir("$outdir/$hs_site2enzyme{$site}");#创建某个酶的文件夹

#电子酶切
&CheckDir("$outdir/$hs_site2enzyme{$site}/enzyme_result");#创建电子酶切结果文件夹

print STDOUT "###Electron digestion start, ",`date`;#输出日志
open LIST,"$list" or die "cannot open $list\n";
my $pm=new Parallel::ForkManager($cpu); #多线程
while(<LIST>){
	my $pid=$pm->start and next;
	my $line=$_;
	if(/^#/ || /^$/){;}else{#去除注释行和空行
		chomp($line);
		my @tmp=split /\t/,$line;
		if(-e $tmp[-1]){#检测文件是否存在
			if(defined($verbose)){#冗长模式
				print STDOUT "[Electron digestion] Analyze $tmp[0]\n";
			}
			system("$Bin/EeTt.pl -i $tmp[-1] -t 1 -s $site -od $outdir/$hs_site2enzyme{$site}/enzyme_result/ -op $tmp[0] -gz no 1> /dev/null; mv $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].$hs_site2enzyme{$site}.fa $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa");#电子酶切，生成$tmp[0].fa
		}else{
			print STDERR "[ERROR] $tmp[-1] does not exist,please check your genome file\n";
			exit;
		}
	}
	$pm->finish;
}
$pm->wait_all_children;
close LIST;
print STDOUT "###Electron digestion complete, ",`date`;#输出日志

#标签记录和鉴定
print STDOUT "###Record tag and Identification label start, ",`date`;#输出日志
for my $level(sort {$hs_type{$a}<=>$hs_type{$b}} keys %hs_type){
	#记录每个标签的分类
	print STDOUT "###($level) Record tag start, ",`date`;#输出日志
	my %hash;
	open LIST,"$list" or die "cannot open $list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);#去除注释行和空行
		my $line=$_;
		chomp($line);
		my @tmp=split /\t/,$line;
		my $class=join("\t",@tmp[1..$hs_type{$level}]);
		if(defined($verbose)){#冗长模式
			print STDOUT "[($level) Record tag] Analyze $tmp[0]\n";
		}
		open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n";
		while(<EN>){
			chomp;
			next if(/^>/);
			my $tag=$_;
			if(exists $hash{$tag}){#标签正向是否存在判断
				$hash{$tag}{$class}++;
			}else{#标签反向互补，并记录
				$tag=~tr/ATCG/TAGC/;
				$tag=reverse($tag);
				$hash{$tag}{$class}++;
			}
		}
		close EN;
	}
	close LIST;
	print STDOUT "###($level) Record tag complete, ",`date`;#输出日志

	print STDOUT "###($level) Identification label start, ",`date`;#输出日志
#	&CheckDir("$outdir/$hs_site2enzyme{$site}/database");
	open OU,"|gzip > $outdir/$hs_site2enzyme{$site}/$level.gz" or die "can not open $outdir/$hs_site2enzyme{$site}/$level.gz\n";
	print OU join("\t",@head[0..$hs_type{$level}]),"\tTags...\n";
	#计算每个基因组在各水平下unique标签数目，并输出每个水平每个基因组unique标签
	open LIST,"$list" or die "cannot open $list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);#去除注释行和空行
		chomp;
		my @tmp=split /\t/;
		if(defined($verbose)){#冗长模式
			print STDOUT "[(level) Identification label] Analyze $tmp[0]\n";
		}
		print OU join("\t",@tmp[0..$hs_type{$level}]);
		open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n";
		while(<EN>){
			chomp;
			next if(/^>/);
			my $detection_tag=$_;#待检测的标签
			unless(exists $hash{$detection_tag}){#判断标签是否正向存在，不存在则进行反向互补
				$detection_tag=~tr/ATCG/TAGC/;
				$detection_tag=reverse($detection_tag);
			}
			if(keys %{$hash{$detection_tag}}==1){
				print OU "\t$detection_tag";
			}
		}
		close EN;
		print OU "\n";
	}
	close LIST;
	close OU;
	undef %hash;
	print STDOUT "###($level) Identification label complete, ",`date`;#输出日志
	#去除冗余标签
	if($remove_redundant eq "yes"){
		print STDOUT "###($level) Remove redundant start, ",`date`;#输出日志
		open IN,"gzip -dc $outdir/$hs_site2enzyme{$site}/$level.gz|" or die "can not open $outdir/$hs_site2enzyme{$site}/$level.gz\n";#打开冗余标签
		open OU,"|gzip > $outdir/$hs_site2enzyme{$site}/$level\_norebundancy.gz" or die "cannot open $outdir/$hs_site2enzyme{$site}/$level\_norebundancy.gz\n"; #输出去除冗余标签后文件
		my $num_tag_col;
		while(<IN>){
			chomp;
			my $line=$_;
			my %hs_unique_tag;
			my @tmp=split /\t/,$line;
			if($.==1){
				$num_tag_col=$#tmp;
				print OU "$line\n";
				next;
			}
			for my $i($num_tag_col .. $#tmp){
				$hs_unique_tag{$tmp[$i]}++;
			}
			print OU join("\t",@tmp[0 .. $num_tag_col-1]);
			for my $i($num_tag_col .. $#tmp){
				if($hs_unique_tag{$tmp[$i]}==1){
					print OU "\t$tmp[$i]";
				}
			}
			print OU "\n";
			undef %hs_unique_tag;
		}
		close IN;
		close OU;
		print STDOUT "###($level) Remove redundant complete, ",`date`;#输出日志
	}
}
print STDOUT "###Record tag and Identification label complete, ",`date`;#输出日志

print STDOUT "###Data stat start, ",`date`;#输出日志
my %stat;
#统计每个基因组总标签数
open LIST,"$list" or die "cannot open $list\n";
while(<LIST>){
	next if(/^#/ || /^$/);#去除注释行和空行
	chomp;
	my @tmp=split /\t/;
	open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n";
	while(<EN>){
		chomp;
		$stat{$tmp[0]}{"all"}++ if(/^>/);
	}
	close EN;
#	system("gzip -f $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa"); #压缩酶切文件
}
close LIST;

#统计每个基因组各水平unique标签数
for my $i(keys %hs_type_database){
	next unless(-e "$outdir/$hs_site2enzyme{$site}/$i.gz");
	open IN,"gzip -dc $outdir/$hs_site2enzyme{$site}/$i.gz|" or die "cannot open $outdir/$hs_site2enzyme{$site}/$i.gz\n";
	while(<IN>){
		chomp;
		my $line=$_;
		my @tmp=split /\t/,$line;
		next if($.==1 && $line=~/^#/);
		$stat{$tmp[0]}{$i}=$#tmp-$hs_type_database{$i};
	}
	close IN;
}

#输出
open OU,">$outdir/$hs_site2enzyme{$site}/stat.xls" or die "cannot open $outdir/$hs_site2enzyme{$site}/stat.xls\n";
print OU "#Unique_Name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecie\tStrain\t";
print OU "All_Tag\tKingdom_Unique\tPhylum_Unique\tClass_Unique\tOrder_Unique\tFamily_Unique\tGenus_Unique\tSpecie_Unique\tStrain_Unique\n";
open LIST,"$list" or die "cannot open $list\n";
while(<LIST>){
	next if(/^#/ || /^$/);#去除注释行和空行
	chomp;
	my @tmp=split /\t/;
	$stat{$tmp[0]}{"all"}=0 unless(exists $stat{$tmp[0]}{"all"});
	print OU join("\t",@tmp[0..$#tmp-1]),"\t",$stat{$tmp[0]}{"all"};#输出表头
	for my $i(sort {$hs_type_database{$a} <=> $hs_type_database{$b}} keys %hs_type_database){
		if(exists $stat{$tmp[0]}{$i}){
			print OU "\t$stat{$tmp[0]}{$i}";
		}else{
			print OU "\t-";
		}
	}
	print OU "\n";
}
close LIST;
close OU;
undef %stat;
print STDOUT "###Data stat complete, ",`date`;#输出日志

print STDOUT "###Delete enzyme result start, ",`date`;#输出日志
system("rm -rf $outdir/$hs_site2enzyme{$site}/enzyme_result/");#删除酶切结果
print STDOUT "###Delete enzyme result complete, ",`date`;#输出日志

print STDOUT "ALL DONE, ",`date`;

sub execute{#打印出并执行命令
	my $cmd = shift;
	print "$cmd\n";
	system($cmd);
}
sub CheckDir{#创建目录
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit;}
		}
		return 1;
}
