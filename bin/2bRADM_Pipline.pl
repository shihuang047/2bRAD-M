#!/usr/bin/env perl
# Authors: Zheng Sun , Rongchao Zhang, Shi Huang
# Last update: 2020.12.22
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use Parallel::ForkManager;
use Cwd 'abs_path';

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="20201222";

#scripts path
my $Bin="$Bin/../scripts/";

#定性参数
my $qual ||="yes";#是否进行定性
my $site1 ||=5;#酶切位点
my $level1 ||="species";#水平

#定量参数
my $quan ||="yes";#是否进行定量
my $g_score_threshold ||=5;#对定性的合并结果，进行分类筛选 gscore阈值
my $GCF_threshold ||=1;#鉴定到某个基因组几个标签以上，该基因组才会被纳入定量建库
my $site2 ||=5;#酶切位点
my $level2 ||="species";#水平

#数据酶切，基因组酶切CPU
my $cpu1 ||=10;
#多酶定性cpu
my $cpu2 ||=8;

#数据质控参数
my $qc ||="yes";#是否进行质控
my $qc_n ||=0.08;
my $qc_q ||=30;
my $qc_p ||=80;
my $qc_b ||=33;

my $mock_sample="";
my $negative_control_sample="";

select STDOUT;$|=1;#标准输出清楚缓存

my($list,$type,$database,$outdir,$help);
GetOptions(
		"t:i"   => \$type,
		"l:s"   => \$list,
		"d:s"   => \$database,
		"o:s"   => \$outdir,
		#初步定性
		"p:s"   => \$qual,
		"s1:s"  => \$site1,
		"t1:s"   => \$level1,

		#精细定量
		"q:s"   => \$quan,
		"gsc:i"   => \$g_score_threshold,#筛选分类
		"gcf:i" => \$GCF_threshold,#筛选分类中的基因组
		"s2:s"  => \$site2,
		"t2:s"  => \$level2,

		#cpu
		"c1:i"  => \$cpu1,
		"c2:i"  => \$cpu2,
		
		#质控参数
		"qc:s"  => \$qc,
		"qcn:f"  => \$qc_n,
		"qcq:i"  => \$qc_q,
		"qcp:i"  => \$qc_p,
		"qcb:i"  => \$qc_b,

		#丰度结果过滤
		"ms:s"  => \$mock_sample,
		"ncs:s" => \$negative_control_sample,

		"h|help:s" => \$help,
		);

sub usage{#帮助
print STDERR "\e[;33;1m
	DESCRIPTION
    	We here provided a streamlined 2bRAD pipeline for analyzing microbial compositions from the 2bRAD/shotgun metagenomics data based on the species-specific 2bRAD markers.
	USAGE
	  perl $0
        PARAMETERS
          -t   <int>    The acceptable types of an input sequencing data file. The file path should be also listed in the sample list file (para -l)
                        [1] generic genome data in a fasta format
                        [2] shotgun metagenomic data in a fastq format(either SE or PE platform is accepted)
                        [3] 2bRAD data from a SE sequencing platform in a fastq format
                        [4] 2bRAD data from a PE sequencing platform in a fastq format
          -l   <file>   The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. Only four formats of a sample list file are accepted and should match with parameter -t: 
                        [1] sample<tab>sample.fa(.gz)
                        [2] sample<tab>shotgun.1.fq(.gz)(<tab>shotgun.2.fq.gz)
                        [3] sample<tab>2bsingle.fq(.gz or 2bsingle.1.fq.gz)
                        [4] sample1<tab>sample2<tab>sample3<tab>sample4<tab>sample5<tab>R1.fq(.gz)<tab>R2.fq(.gz)
          -d   <dir>    The working path of 2B-Tag-DB
          -o   <dir>    The output directory (if it doesn't exist, will be created automatically as 'outdir')
        OPTIONS of Qualitative Analysis
          -p   <str>   If qualitative analysis applies or not [default: $qual] (yes or no)
          -s1  <str>   The enzyme site(s) for the qualitative analysis. One or more sites can be specified(comma separated) [default: $site1]
                       It represents which enzyme(s) will be used for digital restriction digestion, the contruction of 2b-Tag-DB for the following qualitative analysis and quantitative analysis).
                       [1]CspCI  [5]BcgI  [9]BplI     [13]CjePI  [17]AllEnzyme
                       [2]AloI   [6]CjeI  [10]FalI    [14]Hin4I
                       [3]BsaXI  [7]PpiI  [11]Bsp24I  [15]AlfI
                       [4]BaeI   [8]PsrI  [12]HaeIV   [16]BslFI
          -t1  <str>   The taxonomic level for 2bRAD markers in the qualitative database, which should be one of the following: kingdom,phylum,class,order,family,genus,species,strain. [default: $level1]
        OPTIONS of Quantitative Analysis
          -q   <str>   If the quantitative analysis applies or not [default: $quan] (yes or no)
          -gsc <int>   G score threshold for identifying the condidate microbes present in a sample in qualitative analysis, which also determines the membership of sample-specific 2B-Tag-DB database in the quantitative analysis step. [default: $g_score_threshold, it means >$g_score_threshold]
          -gcf <int>   The threshold of the 2bRAD tag number for the presence of a microbial genome (i.e., GCF) in the qualitative analysis, which also determines the membership of sample-specific 2B-Tag-DB database in the quantitative analysis step. [default: $GCF_threshold, it means >$GCF_threshold]
          -s2  <str>   The enzyme site for the quantitative analysis. (refer to -s1) [default: $site2, must be included in para -s1]
          -t2  <str>   The taxonomic level for 2bRAD markers in the quantitative database, which should be one of the following: kingdom,phylum,class,order,family,genus,species,strain. [default: $level2]
        OPTIONS of CPU
          -c1  <int>   The number of CPUs used in the digital digestion step for multiple samples [default: $cpu1]
          -c2  <int>   The number of CPUs used for calculating the abundance for multiple samples [default: $cpu2] (each CPU needs about 15~65G of memory)
        OPTIONS of Quality Control
          -qc  <str>   If quality control apply or not [default: $qc] (yes or no)
          -qcn <float> The maximum ratio of base \"N\" [default: $qc_n]
          -qcq <int>   The minimum quality score to keep [default: $qc_q]
          -qcp <int>   The minimum percentage of bases that must have [-qcq] quality [default: $qc_p]
          -qcb <int>   The quality values of base [default: $qc_b]
        OPTIONS of Abundance Estimation
          -ms  <str>   The mock-community sample name(s) (separated by commas)
          -ncs <str>   The sample name(s) (separated by commas) of negative control that can be used for filtering potential contaminations
          -h|help   Print this help information.
	AUTHOR:  $author $time\e[0m\n";
}

if(defined($help)){
	&usage;
	exit 0;
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
#参数检测
unless($type && $list && $database && $outdir){
	&usage;
	print STDERR "Parameter -t -l -d or -o error.\n";
	exit 1;
}

#转化为绝对路径
$list=abs_path($list);
$database=abs_path($database);
$outdir=abs_path($outdir);

#输入数据类型检测
unless($type==1 || $type==2 || $type==3 || $type==4){
	&usage;
	print STDERR "Parameter -t is wrong.";
	exit 1;
}

#数据库文件检测
unless(-e "$database/abfh_classify_with_speciename.txt.gz"){
	print STDERR "incomplete database, $database/abfh_classify_with_speciename.txt.gz does not exists\n";
	exit 1;
}

#定性参数检测
unless($qual eq "no" || $qual eq "yes"){
	&usage;
	print STDERR "Parameter -p is wrong. Cannot get $qual\n";
	exit 1;
}
#定性鉴定水平检测
unless($level1 eq "kingdom" || $level1 eq "phylum" || $level1 eq "class" || $level1 eq "order" || $level1 eq "family" || $level1 eq "genus" || $level1 eq "species" || $level1 eq "strain"){
	&usage;
	print STDERR "Parameter -t1 is wrong. Cannot get $level1\n";
	exit 1;
}
#定性酶切位点检查及数据库检测
my %hs_site1;
if($site1=~/17/){
	$site1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16";
}
my @site1=split /,/,$site1;
for my $i(@site1){
	$hs_site1{$i}++;
	unless(exists $hs_site2enzyme{$i}){
		&usage;
		print STDERR "parameter -s1 is wrong, $i does not exists\n";
		exit 1;
	}
	#检测定性数据库并检测species数据库，用于定量(findgenome脚本)
	unless(-e "$database/$hs_site2enzyme{$i}.$level1.fa.gz" && -e "$database/$hs_site2enzyme{$i}.species.fa.gz"){
		&usage;
		print STDERR "incomplete database, $database/$hs_site2enzyme{$i}.$level1.fa.gz or $database/$hs_site2enzyme{$i}.species.fa.gz does not exists\n";
		exit 1;
	}
}

#定量参数检测
#定性鉴定水平检测
unless($level2 eq "kingdom" || $level2 eq "phylum" || $level2 eq "class" || $level2 eq "order" || $level2 eq "family" || $level2 eq "genus" || $level2 eq "species" || $level2 eq "strain"){
	&usage;
	print STDERR "Parameter -t2 is wrong. Cannot get $level2\n";
	exit 1;
}
#定性酶切位点检查
my @site2;
if($quan eq "yes"){
	if($site2=~/17/){
		$site2="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16";
	}
	@site2=split /,/,$site2;
	for my $i(@site2){
		unless(exists $hs_site2enzyme{$i}){#检测site2输入是否正确
			&usage;
			print STDERR "parameter -s2 is wrong, $i does not exists\n";
			exit 1;
		}
		unless(exists $hs_site1{$i}){#检测site2是否包含于site1
			&usage;
			print STDERR "parameter -s2 is wrong, $i is not included in para -s1\n";
			exit 1;
		}
	}
}elsif($quan eq "no"){
	;
}else{
	&usage;
	print STDERR "parameter -q is wrong\n";
	exit 1;
}

#质控参数检查
unless($qc eq "yes" || $qc eq "no"){
	&usage;
	print STDERR "Parameter -qc is wrong. Cannot get $qc\n";
}

print STDOUT "###COMMAND: perl $0 -t $type -l $list -d $database -o $outdir -p $qual -s1 $site1 -t1 $level1 -q $quan -gsc $g_score_threshold -gcf $GCF_threshold -s2 $site2 -t2 $level2 -c1 $cpu1 -c2 $cpu2 -qc $qc -qcn $qc_n -qcq $qc_q -qcp $qc_p -qcb $qc_b -ms $mock_sample -ncs $negative_control_sample\n";
&CheckDir($outdir);
#数据酶切
&CheckDir("$outdir/enzyme_result");
print STDOUT "###Electron digestion start, ",`date`;
open LIST,"$list" or die "cannot open $list\n";
my $pm=new Parallel::ForkManager($cpu1); #多线程
while(<LIST>){
	my $line=$_;
	if(/^#/ || /^$/){;}else{#去除注释行和空行
		chomp($line);
		my @tmp=split /\t/,$line;
		for my $i(@site1){
			my $pid=$pm->start and next;
			if($type!=4){#除2brad五标签之外其他数据处理
				unless(-e "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz"){
					&execute("perl $Bin/2bRADExtraction.pl -i @tmp[1..$#tmp] -t $type -s $i -od $outdir/enzyme_result -op $tmp[0] -gz yes -qc $qc -n $qc_n -q $qc_q -p $qc_p -b $qc_b -fm fa 1>/dev/null");
					`sleep 1s`;
				}else{
					print STDOUT "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz already exists, skip.\n";
					`sleep 1s`;
				}
			}else{#2brad五标签处理
	 			unless(-e "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[1].$hs_site2enzyme{$i}.fa.gz" &&
					   -e "$outdir/enzyme_result/$tmp[2].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[3].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[4].$hs_site2enzyme{$i}.fa.gz" ){
					&execute("perl $Bin/2bRADExtraction.pl -i $tmp[-2] $tmp[-1] -t $type -s $i -od $outdir/enzyme_result -op @tmp[0..4] -gz yes -qc $qc -n $qc_n -q $qc_q -p $qc_p -b $qc_b -fm fa 1>/dev/null");
				}else{
					print STDOUT "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[1].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[2].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[3].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[4].$hs_site2enzyme{$i}.fa.gz already exist, skip\n";
					`sleep 1s`;
				}
			}
			$pm->finish;
		}
	}
}
$pm->wait_all_children;
close LIST;
print STDOUT "###Electron digestion complete, ",`date`;


##整理列表
&CheckDir("$outdir/list");
for my $i(@site1){
	open LIST,"$list" or die "cannot open $list\n";
	open OU,">$outdir/list/$hs_site2enzyme{$i}.list" or die "cannot open $outdir/list/$hs_site2enzyme{$i}.list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);
		chomp;
		my @tmp=split /\t/;
		if($type!=4){#除2brad五标签之外其他数据处理
			print OU "$tmp[0]\t$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz\n";
		}else{#2brad五标签处理
			print OU "$tmp[0]\t$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz\n";
			print OU "$tmp[1]\t$outdir/enzyme_result/$tmp[1].$hs_site2enzyme{$i}.fa.gz\n";
			print OU "$tmp[2]\t$outdir/enzyme_result/$tmp[2].$hs_site2enzyme{$i}.fa.gz\n";
			print OU "$tmp[3]\t$outdir/enzyme_result/$tmp[3].$hs_site2enzyme{$i}.fa.gz\n";
			print OU "$tmp[4]\t$outdir/enzyme_result/$tmp[4].$hs_site2enzyme{$i}.fa.gz\n";
		}
	}
	close LIST;
	close OU;
}
if($type==4){#2brad五标签样品名格式行转化成列
	open LIST,"$list" or die "cannot open $list\n";
	open OULI,">$outdir/list/2brad_5tag.list" or die "cannot open $outdir/list/2brad_5tag.list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);#去除注释行和空行
		chomp;
		my @tmp=split /\t/;
		for my $i(0..4){
			print OULI "$tmp[$i]\t$outdir/quantitative/$tmp[$i]/$tmp[$i].combine.xls\n";
		}
	}
	close LIST;
	close OULI;
}else{
	open LIST,"$list" or die "cannot open $list\n";
	open OULI,">$outdir/list/Abundance_Stat.list" or die "cannot open $outdir/list/Abundance_Stat.list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);#去除注释行和空行
		chomp;
		my @tmp=split /\t/;
		print OULI "$tmp[0]\t$outdir/quantitative/$tmp[0]/$tmp[0].combine.xls\n";
	}
	close LIST;
	close OULI;
}

##多线程初步定性
if($qual eq "yes"){#是否需要定性
	print STDOUT "###Qualitative start, ",`date`;
	&CheckDir("$outdir/qualitative");
	$pm=new Parallel::ForkManager($cpu2);
	for my $i(@site1){
		my $pid=$pm->start and next;
		&execute("perl $Bin/CalculateRelativeAbundance_Single2bEnzyme.pl -l $outdir/list/$hs_site2enzyme{$i}.list -d $database/ -t $level1 -s $i -o $outdir/qualitative -g 0 -v yes 1> /dev/null");#未对G_score过滤
		`sleep 1s`;
		$pm->finish;
	}
	$pm->wait_all_children;
	
	##定性结果合并
	if($type!=4){#除2brad五标签之外其他数据处理
		&execute("perl $Bin/CalculateRelativeAbundance_Combined2bEnzymes.pl -l $list -s $site1 -io $outdir/qualitative -m combine -g 0");#未对G_score过滤
	}else{#2brad五标签处理
		&execute("perl $Bin/CalculateRelativeAbundance_Combined2bEnzymes.pl -l $outdir/list/2brad_5tag.list -s $site1 -io $outdir/qualitative -m combine -g 0");#未对G_score过滤
	}
	print STDOUT "###Qualitative complete, ",`date`;
}else{
	print STDOUT "All Done, ",`date`;
	exit 0;
}


if($quan eq "no"){
	print STDOUT "All Done, ",`date`;
	exit 0;
}
#精细定量
print STDOUT "###Quantitative start, ",`date`;

&CheckDir("$outdir/quantitative_sdb");
&CheckDir("$outdir/quantitative");

#FindGenome_ByQualitative
if($type!=4){#除2brad五标签之外其他数据处理
	&execute("perl $Bin/FindGenome_ByQualitative.pl -l $list -d $database -o $outdir/quantitative_sdb -qualdir $outdir/qualitative -gscore $g_score_threshold -gcf $GCF_threshold");
}else{#2brad五标签处理
	&execute("perl $Bin/FindGenome_ByQualitative.pl -l $outdir/list/2brad_5tag.list -d $database -o $outdir/quantitative_sdb -qualdir $outdir/qualitative -gscore $g_score_threshold -gcf $GCF_threshold");
}

if($type!=4){#除2brad五标签之外其他数据处理
	open LIST,"$list" or die "cannot open $list\n";
}else{#2brad五标签处理
	open LIST,"$outdir/list/2brad_5tag.list" or die "cannot open $outdir/list/2brad_5tag.list\n";
}
while(<LIST>){
	next if(/^#/ || /^$/);
	chomp;
	my $sample_name=(split /\t/)[0];
	print STDOUT "Analysis $sample_name, ",`date`;
	&CheckDir("$outdir/quantitative_sdb/$sample_name/database");
	&execute("cp $outdir/quantitative_sdb/$sample_name/sdb.list $outdir/quantitative_sdb/$sample_name/database/abfh_classify_with_speciename.txt && gzip -f $outdir/quantitative_sdb/$sample_name/database/abfh_classify_with_speciename.txt");
	#精细定量开始
	$pm=new Parallel::ForkManager($cpu2);
	for my $i(@site2){
		my $pid=$pm->start and next;
		open SA,">$outdir/quantitative_sdb/$sample_name/$hs_site2enzyme{$i}.list" or die "cannot open $outdir/quantitative_sdb/$sample_name/$hs_site2enzyme{$i}.list\n";
		print SA "$sample_name\t$outdir/enzyme_result/$sample_name.$hs_site2enzyme{$i}.fa.gz\n";
		close SA;
		&execute("perl $Bin/CreateQuanDatabase_2bRAD.pl -l $outdir/quantitative_sdb/$sample_name/sdb.list -e $database/$hs_site2enzyme{$i}.species.fa.gz -s $i -t $level2 -o $outdir/quantitative_sdb/$sample_name/database -r no 1> /dev/null");#建库
		&execute("perl $Bin/CalculateRelativeAbundance_Single2bEnzyme.pl -l $outdir/quantitative_sdb/$sample_name/$hs_site2enzyme{$i}.list -d $outdir/quantitative_sdb/$sample_name/database -t $level2 -s $i -o $outdir/quantitative -g 0 -v yes 1> /dev/null");#定量 不对gscore进行过滤
		$pm->finish;
	}
	$pm->wait_all_children;
	&execute("rm -rf $outdir/quantitative_sdb/$sample_name/database");#删除数据库
}
close LIST;

#精细定量结果合并
if($type!=4){#除2brad五标签之外其他数据处理
	&execute("perl $Bin/CalculateRelativeAbundance_Combined2bEnzymes.pl -l $list -s $site2 -io $outdir/quantitative -m combine -g 0 1> /dev/null");#不对gscore过滤
}else{#2brad五标签处理
	&execute("perl $Bin/CalculateRelativeAbundance_Combined2bEnzymes.pl -l $outdir/list/2brad_5tag.list -s $site2 -io $outdir/quantitative -m combine -g 0 1> /dev/null");#不对gscore过滤
}
print STDOUT "###Quantitative complete, ",`date`;

print STDOUT "###Abundance_Stat start, ",`date`;
if($type!=4){#除2brad五标签之外其他数据处理
	&execute("perl $Bin/Abundance_Stat.pl -l $outdir/list/Abundance_Stat.list -o $outdir/quantitative -p Abundance_Stat -m $mock_sample -c $negative_control_sample");
}else{#2brad五标签处理
	&execute("perl $Bin/Abundance_Stat.pl -l $outdir/list/2brad_5tag.list -o $outdir/quantitative -p Abundance_Stat -m $mock_sample -c $negative_control_sample");
}

print STDOUT "###Abundance_Stat complete, ",`date`;

print STDOUT "All Done, ",`date`;

sub CheckDir{#创建目录
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}


sub execute{#打出命令并执行
	my $cmd = shift;
	print STDOUT "$cmd\n";
	my $exit_code=system($cmd);
	if($exit_code!=0){
		print STDERR "Command $cmd failed with an exit code of $exit_code.\n";
		exit($exit_code >> 8);
	}
}
