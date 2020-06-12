#!/usr/bin/env perl
#Writer: Sunzheng , Zhangrongchao 2019.10.21
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use Parallel::ForkManager;

my $author="Sunzheng, Zhangrongchao";
my $time="20200301";

#定性参数
my $qual ||="yes";#是否进行定性
my $site1 ||=5;#酶切位点
my $level1 ||="specie";#水平

#定量参数
my $quan ||="yes";#是否进行定量
my $g_score_threshold ||=0;#对定性的合并结果，进行分类筛选 gscore阈值
my $GCF_threshold ||=1;#鉴定到某个基因组几个标签以上，该基因组才会被纳入定量建库
my $site2 ||=5;#酶切位点
my $level2 ||="specie";#水平

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


select STDOUT;$|=1;#标准输出清楚缓存

my($list,$type,$database,$outdir);
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
		"gscore:i"   => \$g_score_threshold,#筛选分类
		"gcf:i" => \$GCF_threshold,#筛选分类中的基因组
		"s2:s"  => \$site2,
		"t2:s"  => \$level2,

		#cpu
		"c1:i"  => \$cpu1,
		"c2:i"  => \$cpu2,
		
		#质控参数
		"qc:s"  => \$qc,
		"qc_n:f"  => \$qc_n,
		"qc_q:i"  => \$qc_q,
		"qc_p:i"  => \$qc_p,
		"qc_b:i"  => \$qc_b,
		);

sub usage{#帮助
print STDERR "\e[;33;1m
	DESCRIPTION
	  shotgun/2brad pipline
	USAGE
	  perl $0
	PARAMETERS
	  -t  <i> Type of Input File in sample list(para -l)
	          [1] Genome Data in Fasta Format
	          [2] Shotgun Data in Fastq Format(SE or PE)
	          [3] SE Platform Data in Fastq Format
	          [4] PE Platform Data in Fastq Format
	  -l  <s> sample list (the line which begins with # will be ignored)
	          [1] sample<tab>sample.fa(.gz)
	          [2] sample<tab>shotgun.1.fq(.gz)(<tab>shotgun.2.fq.gz)
	          [3] sample<tab>2bsingle.fq(.gz or 2bsingle.1.fq.gz)
	          [4] sample1<tab>sample2<tab>sample3<tab>sample4<tab>sample5<tab>R1.fq(.gz)<tab>R2.fq(.gz)
	  -d  <s> database path
	  -o  <s> outdir (if not exists,it will be created)
	OPTIONS of Qualitative Analysis
	  -p  <s> qualitative or not [$qual] (yes or no)
	  -s1 <s> qualitative enzyme site. One or more of site. (comma separated) [$site1]
	          It also represents enzymatic digestion or data splitting, and combining qualitative analysis results(for quantitative analysis).
	          [1]CspCI  [5]BcgI  [9]BplI      [13]CjePI  [17]AllEnzyme
	          [2]AloI   [6]CjeI  [10]FalI    [14]Hin4I
	          [3]BsaXI  [7]PpiI  [11]Bsp24I  [15]AlfI
	          [4]BaeI   [8]PsrI  [12]HaeIV   [16]BslFI
	  -t1 <s> qualitative database level. One of kingdom,phylum,class,order,family,genus,specie,strain. [$level1]
	OPTIONS of Quantitative Analysis
	  -q  <s> quantitative or not [$quan] (yes or no)
	  -gscore <i> G score threshold of classify in qualitative analysis, it decides quantitative database. [$g_score_threshold, it means >$g_score_threshold]
	  -gcf <i> detected tag threshold of GCF in qualitative analysis, it decides quantitative database. [$GCF_threshold, it means >$GCF_threshold]
	  -s2 <s> quantitative enzyme site (refer to -s1) [$site2, must be included in para -s1]
	  -t2 <s> quantitative database level. One of kingdom,phylum,class,order,family,genus,specie,strain. [$level2]
	OPTIONS of CPU
	  -c1 <i> enzyme cpu [$cpu1]
	  -c2 <i> calculate cpu [$cpu2] (each CPU needs about 15~65G of memory)
	OPTIONS of Quality Control
	  -qc   <s> quality control or not [yes] (yes or no)
	  -qc_n <f> Maximum Ratio of Base \"N\" [$qc_n]
	  -qc_q <i> Minimum Quality Score to Keep [$qc_q]
	  -qc_p <i> Minimum Percent of Bases that must have [-q] Quality [$qc_p]
	  -qc_b <i> Quality Values Base [$qc_b]
	AUTHOR:  $author $time\e[0m\n";
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
	exit;
}

#输入数据类型检测
unless($type==1 || $type==2 || $type==3 || $type==4){
	&usage;
	print STDERR "Parameter -t is wrong.";
	exit;
}

#数据库文件检测
unless(-d "$database/database" && -e "$database/abfh_classify.txt" && -d "$database/genome_ref"){
	print STDERR "incomplete database, please check\n";
	exit;
}

#定性参数检测
unless($qual eq "no" || $qual eq "yes"){
	&usage;
	print STDERR "Parameter -p is wrong. Cannot get $qual\n";
	exit;
}
#定性鉴定水平检测
unless($level1 eq "kingdom" || $level1 eq "phylum" || $level1 eq "class" || $level1 eq "order" || $level1 eq "family" || $level1 eq "genus" || $level1 eq "species" || $level1 eq "strain"){
	&usage;
	print STDERR "Parameter -t is wrong. Cannot get $level1\n";
	exit;
}
#定性酶切位点检查及数据库检测
if($site1=~/17/){
	$site1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16";
}
my @site1=split /,/,$site1;
for my $i(@site1){
	unless(exists $hs_site2enzyme{$i}){
		&usage;
		print STDERR "parameter -s1 is wrong, $i does not exists\n";
		exit;
	}
	unless(-e "$database/database/$hs_site2enzyme{$i}/$level1.gz"){
		&usage;
		print STDERR "incomplete database, $database/database/$hs_site2enzyme{$i}/$level1.gz does not exists";
		exit;
	}
}

#定量参数检测
#定性鉴定水平检测
unless($level2 eq "kingdom" || $level2 eq "phylum" || $level2 eq "class" || $level2 eq "order" || $level2 eq "family" || $level2 eq "genus" || $level2 eq "species" || $level2 eq "strain"){
	&usage;
	print STDERR "Parameter -t is wrong. Cannot get $level2\n";
	exit;
}
#定性酶切位点检查
my @site2;
if($quan eq "yes"){
	if($site2=~/17/){
		$site2="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16";
	}
	@site2=split /,/,$site2;
	for my $i(@site2){
		unless(exists $hs_site2enzyme{$i}){
			&usage;
			print STDERR "parameter -s2 is wrong, $i does not exists\n";
			exit;
		}
	}
}elsif($quan eq "no"){
	;
}else{
	&usage;
	print STDERR "parameter -q is wrong\n";
	exit
}

#质控参数检查
unless($qc eq "yes" || $qc eq "no"){
	&usage;
	print STDERR "Parameter -qc is wrong. Cannot get $qc\n";
}

print STDOUT "###COMMAND: perl $0 -t $type -l $list -d $database -o $outdir -p $qual -s1 $site1 -t1 $level1 -q $quan -gscore $g_score_threshold -gcf $GCF_threshold -s2 $site2 -t2 $level2 -c1 $cpu1 -c2 $cpu2 -qc $qc -qc_n $qc_n -qc_q $qc_q -qc_p $qc_p -qc_b $qc_b\n";
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
					&execute("$Bin/EeTt.pl -i @tmp[1..$#tmp] -t $type -s $i -od $outdir/enzyme_result -op $tmp[0] -gz yes -qc $qc -n $qc_n -q $qc_q -p $qc_p -b $qc_b -fm fa 1>/dev/null");
					`sleep 1s`;
				}else{
					print STDOUT "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz already exists\n";
					`sleep 1s`;
				}
			}else{#2brad五标签处理
	 			unless(-e "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[1].$hs_site2enzyme{$i}.fa.gz" &&
					   -e "$outdir/enzyme_result/$tmp[2].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[3].$hs_site2enzyme{$i}.fa.gz" && 
					   -e "$outdir/enzyme_result/$tmp[4].$hs_site2enzyme{$i}.fa.gz" ){
					&execute("$Bin/EeTt.pl -i $tmp[-2] $tmp[-1] -t $type -s $i -od $outdir/enzyme_result -op @tmp[0..4] -gz yes -qc $qc -n $qc_n -q $qc_q -p $qc_p -b $qc_b -fm fa 1>/dev/null");
				}else{
					print STDOUT "$outdir/enzyme_result/$tmp[0].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[1].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[2].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[3].$hs_site2enzyme{$i}.fa.gz && ";
					print STDOUT "$outdir/enzyme_result/$tmp[4].$hs_site2enzyme{$i}.fa.gz already exist\n";
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
			print OULI "$tmp[$i]\n";
		}
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
		&execute("$Bin/Calculate_Tag_Percent.pl -l $outdir/list/$hs_site2enzyme{$i}.list -d $database/ -t $level1 -s $i -o $outdir/qualitative -g 0 -v yes 1> /dev/null");#未对G_score过滤
		`sleep 1s`;
		$pm->finish;
	}
	$pm->wait_all_children;
	
	##定性结果合并
	if($type!=4){#除2brad五标签之外其他数据处理
		&execute("$Bin/calculate_combine.pl -l $list -s $site1 -io $outdir/qualitative -m combine -g 0");#未对G_score过滤
	}else{#2brad五标签处理
		&execute("$Bin/calculate_combine.pl -l $outdir/list/2brad_5tag.list -s $site1 -io $outdir/qualitative -m combine -g 0");#未对G_score过滤
	}
	print STDOUT "###Qualitative complete, ",`date`;
}

#数据库读取
my %gcf2classify_path;
open ABFH,"$database/abfh_classify.txt" or die "cannot open $database/abfh_classify.txt\n";
while(<ABFH>){
	chomp;
	next if(/^#unique_name/); #unique_name    kingdom phylum  class   order   family  genus   species  strain  genome_pat
	my @tmp=split /\t/;
	$gcf2classify_path{$tmp[0]}=$_;
}
close ABFH;

=pod
my %hs_anno;
open AN,"$database/species_ID_annotation.txt" or die "cannot open $database/species_ID_annotation.txt\n";
while(<AN>){
	chomp;
	next if(/^#/);
	my @tmp=split /\t/;
#	push @{$hs_anno{$tmp[2]}},"$tmp[0]\t$tmp[-1]";
	$hs_anno{$tmp[0]}=$tmp[-1];
}
close AN;
=cut;


if($quan eq "no"){
	print STDOUT "All Done, ",`date`;
	exit;
}
#精细定量
print STDOUT "###Quantitative start, ",`date`;
&CheckDir("$outdir/quantitative");
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
	#样品定量列表
	#单样品重建库list准备
	#记录通过G_score阈值的分类
	my (%hs_pass_Gscore_class);
	if(-e "$outdir/qualitative/$sample_name/$sample_name.combine.xls"){#如果没有定性合并结果，不进行后续分析
		open XI,"$outdir/qualitative/$sample_name/$sample_name.combine.xls" or die "cannot open $outdir/qualitative/$sample_name/$sample_name.combine.xls\n";
	}else{
		print STDERR "!!!$sample_name does not have $outdir/qualitative/$sample_name/$sample_name.combine.xls, can't do quantitative analysis\n";
		next;
	}
	while(<XI>){
		chomp;
		my @tmp=split /\t/;
		next if(/^#/);
		my $class=join("\t",@tmp[0..$#tmp-8]);
		$hs_pass_Gscore_class{$class}++ if($tmp[-1]>$g_score_threshold);#通过gscore阈值的分类
	}
	close XI;
	&CheckDir("$outdir/quantitative/$sample_name\_sdb");
	#输出通过GCF标签检测阈值且通过Gscore阈值的GCF基因组
	open OU,"|sort|uniq > $outdir/quantitative/$sample_name\_sdb/sdb.list" or die "cannot open $outdir/quantitative/$sample_name\_sdb/sdb.list\n";
	for my $i(@site1){
		if(-e "$outdir/qualitative/$sample_name/$sample_name.$hs_site2enzyme{$i}.GCF_detected.xls"){
			open QU,"$outdir/qualitative/$sample_name/$sample_name.$hs_site2enzyme{$i}.GCF_detected.xls" or die "cannot open $outdir/qualitative/$sample_name/$sample_name.$hs_site2enzyme{$i}.GCF_detected.xls\n";
		}else{
			print STDERR "warning: $sample_name does not have $outdir/qualitative/$sample_name/$sample_name.$hs_site2enzyme{$i}.GCF_detected.xls\n";
			next;
		}
		while(<QU>){
			chomp;
			my @tmp=split /\t/;
			my $class=join("\t",@tmp[0..$#tmp-4]);
			if(exists $hs_pass_Gscore_class{$class} && $tmp[-2]>$GCF_threshold){#基因组测到的标签种类数大于阈值，并且该基因组对应的分类通过了gscore阈值
				my @all_class=split /\t/,$gcf2classify_path{$tmp[-4]};
				print OU join("\t",@all_class[0..8]),"\t$database/genome_ref/$all_class[-1]\n";
			}
		}
		close QU;
	}
	close OU;
	
	#精细定量开始
	for my $i(@site2){
		open SA,">$outdir/quantitative/$sample_name\_sdb/$hs_site2enzyme{$i}.list" or die "cannot open $outdir/quantitative/$sample_name\_sdb/$hs_site2enzyme{$i}.list\n";
		print SA "$sample_name\t$outdir/enzyme_result/$sample_name.$hs_site2enzyme{$i}.fa.gz\n";
		close SA;
		&execute("$Bin/makedatabase.pl -l $outdir/quantitative/$sample_name\_sdb/sdb.list -s $i -t $level2 -o $outdir/quantitative/$sample_name\_sdb/database -c $cpu1 1> /dev/null");#建库
		&execute("$Bin/Calculate_Tag_Percent.pl -l $outdir/quantitative/$sample_name\_sdb/$hs_site2enzyme{$i}.list -d $outdir/quantitative/$sample_name\_sdb/ -t $level2 -s $i -o $outdir/quantitative -g 0 -v yes 1> /dev/null");#定量 不对gscore进行过滤
	}
	&execute("rm -rf $outdir/quantitative/$sample_name\_sdb/database");#删除数据库
}
close LIST;

#精细定量结果合并
if($type!=4){#除2brad五标签之外其他数据处理{
	&execute("$Bin/calculate_combine.pl -l $list -s $site2 -io $outdir/quantitative -m combine -g 0 1> /dev/null");#不对gscore过滤
}else{#2brad五标签处理
	&execute("$Bin/calculate_combine.pl -l $outdir/list/2brad_5tag.list -s $site2 -io $outdir/quantitative -m combine -g 0 1> /dev/null");#不对gscore过滤
}
print STDOUT "###Quantitative complete, ",`date`;

print STDOUT "All Done, ",`date`;

sub CheckDir{#创建目录
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit;}
	}
	return 1;
}

sub execute{
	my $cmd = shift;
	print "$cmd\n";
	system($cmd);
}
