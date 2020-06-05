#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="2020.06.03";

#set default parameters
my $g_score_threshold ||=0;
my $verbose ||="yes";

select STDOUT;$|=1;#Standard output for clearing cache

my ($list,$database,$site,$outdir,$level);
GetOptions(
		"l:s" => \$list,
		"d:s" => \$database,
		"t:s" => \$level,
		"s:s" => \$site,
		"o:s" => \$outdir,

		"g:i" => \$g_score_threshold,
		"v:s" => \$verbose,
		);


sub usage{# help information
print STDERR "\e[;33;1m
DESCRIPTION
    It computes the relative abundance of taxa identified from each of 2b-RAD samples using a precalcuated taxa-specific 2b-RAD reference database by a single type 2b restriction enzyme.
USAGE
    perl $0
Required:
    -l <str>    The path of the input filepath list (the lines that begin with # will be ignored) e.g: sample_name<tab>data_path(fa|fq)(.gz).
    -d <str>    The database filepath.
    -t <str>    The taxonomy level of the taxa-specific 2b-RAD database used. It should be one of the following: kingdom,phylum,class,order,family,genus,specie,strain.
    -s <str>    One of the type 2b restriction enzymes (sites).
                [1]CspCI  [9]BplI
                [2]AloI   [10]FalI
                [3]BsaXI  [11]Bsp24I
                [4]BaeI   [12]HaeIV
                [5]BcgI   [13]CjePI
                [6]CjeI   [14]Hin4I
                [7]PpiI   [15]AlfI
                [8]PsrI   [16]BslFI
    -o <str>    The output directory (automatically create if it does not exist)
Optional:
    -g <int>    The threshold of G score [$g_score_threshold, it means >=$g_score_threshold]. To control the false-positive in the species identification, G score was derived for each speciesidentified within a sample, which is a harmonious mean of read coverage of 2b-RAD markers belongs to a species and number of all possible 2b-RAD markers of this species. Therecommended/default threshold is 5.
    -v <str>    This specify if more detailed information will be shown [$verbose] (yes or no)
AUTHOR:  $author $time\e[0m\n";
}


my %hs_site2enzyme=(# the codes for all restriction enzymes
	'1'  =>  'CspCI',   '2'  =>  'AloI',
	'3'  =>  'BsaXI',   '4'  =>  'BaeI',
	'5'  =>  'BcgI',    '6'  =>  'CjeI',
	'7'  =>  'PpiI',    '8'  =>  'PsrI',
	'9'  =>  'BplI',    '10' =>  'FalI',
	'11' =>  'Bsp24I',  '12' =>  'HaeIV',
	'13' =>  'CjePI',   '14' =>  'Hin4I',
	'15' =>  'AlfI',    '16' =>  'BslFI',
	);

unless($list && $database && $level && $site && $outdir){
	&usage;
	exit;
}
# parameter checking
unless($verbose eq "yes" || $verbose eq "no"){
	&usage;
	print STDERR "Parameter -v is wrong\n";
	exit;
}
# check the taxonomic level of a 2b-RAD reference genome database
unless($level eq "kingdom" || $level eq "phylum" || $level eq "class" || $level eq "order" || $level eq "family" || $level eq "genus" || $level eq "specie" || $level eq "strain"){
	&usage;
	print STDERR "Parameter -t is wrong. Cannot get $level\n";
	exit;
}
# check the parameter -s and -e
unless(exists $hs_site2enzyme{$site}){
	&usage;
	print STDERR "Parameter -s $site is wrong\n";
	exit;
}
unless(-e "$database/database/$hs_site2enzyme{$site}/$level.gz"){
	&usage;
	print STDERR "$database/database/$hs_site2enzyme{$site}/$level.gz does not exist\n";
	exit;
}

print STDOUT "COMMAND: perl $0 -l $list -d $database -t $level -s $site -o $outdir -g $g_score_threshold -v $verbose\n";

&CheckDir($outdir);

# load the database
print STDOUT "### Loading the database, $database/database/$hs_site2enzyme{$site}/$level.gz, ",`date`;
my (%hs_tag2GCF,%hs_GCF2class,%hs_tag_theory_num,$head,$tag_col_start);
open IN,"gzip -dc $database/database/$hs_site2enzyme{$site}/$level.gz|" or die "cannot open $database/database/$hs_site2enzyme{$site}/$level.gz\n";
while(<IN>){
	chomp(my $line=$_);
	my @tmp=split /\t/,$line;
	if($.==1 && $line=~/^#Name/){
		for my $i(0..$#tmp){
			$tag_col_start=$i if($tmp[$i]=~/Tags\.\.\./);
		}
		$head=join("\t",@tmp[1..$tag_col_start-1]);
		next;
	}
	my $class=join("\t",@tmp[1..$tag_col_start-1]);#分类id
	$hs_GCF2class{$tmp[0]}=$class;#GCF对应的分类信息
	for my $i($tag_col_start .. $#tmp){
		push @{$hs_tag2GCF{$tmp[$i]}},$tmp[0];#标签对应的GCF
		$hs_tag_theory_num{$class}{$tmp[0]}{$tmp[$i]}++;#每个分类每个GCF每个标签的数目
	}
}
close IN;
print STDOUT "###Loading database completed, ",`date`;

#样品处理
open LI,"$list" or die "cannot open $list\n";
while(<LI>){
	next if(/^#/ || /^$/);#去除注释行和空行
	chomp;
	my ($sample_name,$sample_data)=split /\t/;
	print STDOUT "###($sample_name) Sample identification started, ",`date`;
	my (%hs_tag_num,%hs_detected_GCF_tag);
	#样品数据读取
	if($sample_data=~/\.gz$/){
		open IN,"gzip -dc $sample_data|" or die "cannot open $sample_data\n";
	}else{
		open IN,"$sample_data" or die "cannot open $sample_data\n";
	}
	while(<IN>){
		my $line=$_;
		if($line=~/^@/){#fastq
			$line .=<IN> . <IN> . <IN>;
		}elsif($line=~/^>/){#fasta
			$line .=<IN>;
		}
		my $tag=(split /\n/,$line)[1];
		if(exists $hs_tag2GCF{$tag}){
			my $class=$hs_GCF2class{$hs_tag2GCF{$tag}[0]};
			$hs_tag_num{$class}{$tag}++;#实际样品标签深度
			for my $i(0..$#{$hs_tag2GCF{$tag}}){
				$hs_detected_GCF_tag{$class}{$hs_tag2GCF{$tag}[$i]}{$tag}=$hs_tag_theory_num{$class}{$hs_tag2GCF{$tag}[$i]}{$tag};
			}
		}else{#反向互补
			$tag=~s/ATCG/TAGC/;
			$tag=reverse($tag);
			if(exists $hs_tag2GCF{$tag}){
				my $class=$hs_GCF2class{$hs_tag2GCF{$tag}[0]};
				$hs_tag_num{$class}{$tag}++;#实际样品标签深度
				for my $i(0..$#{$hs_tag2GCF{$tag}}){
					$hs_detected_GCF_tag{$class}{$hs_tag2GCF{$tag}[$i]}{$tag}=$hs_tag_theory_num{$class}{$hs_tag2GCF{$tag}[$i]}{$tag};
				}
			}
		}
	}
	close IN;

	if((keys %hs_tag_num)==0){#不存在微生物样品
		print STDERR "!!!($sample_name) Warning: $hs_site2enzyme{$site}  $level the number of 2b-RAD tags is zero\n";
		print STDOUT "###($sample_name) Sample idenfication completed, ",`date`;
		next;
	}
	&CheckDir("$outdir/$sample_name");#建立文件夹
	#输出各GCF测到的标签数种类数
	open DE,">$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.GCF_detected.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.GCF_detected.xls\n";
	for my $class(sort {$a cmp $b} keys %hs_detected_GCF_tag){
		for my $GCF(sort {$a cmp $b} keys %{$hs_detected_GCF_tag{$class}}){
			my $GCF_all_theory_num;
			my $detected_tag_num;
			$GCF_all_theory_num=keys %{$hs_tag_theory_num{$class}{$GCF}};#理论上 基因组unique标签种类数
			$detected_tag_num=keys %{$hs_detected_GCF_tag{$class}{$GCF}};#实际样品 测到的基因组unique标签种类数
			my $percent=sprintf "%.4f",$detected_tag_num/$GCF_all_theory_num;#测到的标签比例
			print DE "$class\t$GCF\t$GCF_all_theory_num\t$detected_tag_num\t$percent\n";
		}
	}
	close DE;
	# Output
	&CheckDir("$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}") if($verbose eq "yes");
	open OU,">$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.xls\n";
	print OU "$head\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\t";
	print OU "Sequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\tSequenced_Reads_Num/Sequenced_Tag_Num\tSequenced_Tag_Num(depth>1)\t";
	print OU "G_Score\n";
	for my $class(keys %hs_tag_num){
		if($verbose eq "yes"){
			my $output_name=(split /\t/,$class)[-1];
			open DETAIL,">$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}/$output_name.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}/$output_name.xls\n";
		}
		my ($Theoretical_Tag_Num,$Sequenced_Tag_Num,$Sequenced_Tag_Num_2,$Sequenced_Tag_Num2Theoretical_Tag_Num);
		my ($Sequenced_Reads_Num,$Sequenced_Reads_Num2Theoretical_Tag_Num,$Sequenced_Reads_Num2Sequenced_Tag_Num);
		my ($G_Score);
		$Sequenced_Tag_Num=$Sequenced_Reads_Num=$Sequenced_Tag_Num_2=0;
		for my $tag(keys %{$hs_tag_num{$class}}){#每个标签
			$Sequenced_Tag_Num++;
			$Sequenced_Tag_Num_2++ if($hs_tag_num{$class}{$tag}>1);#标签深度>1的标签数统计
			$Sequenced_Reads_Num+=$hs_tag_num{$class}{$tag};
			if($verbose eq "yes"){
				print DETAIL "$tag\t$hs_tag_num{$class}{$tag}\n";
			}
		}
		if($verbose eq "yes"){
			close DETAIL;
		}
		#平均理论标签数
		my $species_all_theory_num;
		for my $GCF(keys %{$hs_tag_theory_num{$class}}){
			for my $tag(keys %{$hs_tag_theory_num{$class}{$GCF}}){
				$species_all_theory_num+=$hs_tag_theory_num{$class}{$GCF}{$tag};
			}
		}
		$Theoretical_Tag_Num=$species_all_theory_num/(keys %{$hs_tag_theory_num{$class}});#平均理论标签数
		$Sequenced_Tag_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Tag_Num/$Theoretical_Tag_Num*100;#测到的标签占理论的百分比
		$Sequenced_Reads_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Theoretical_Tag_Num;#测到的标签深度/理论标签数
		$Sequenced_Reads_Num2Sequenced_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Sequenced_Tag_Num;#测到的标签平均深度
		$G_Score=sprintf "%.8f",sqrt($Sequenced_Tag_Num*$Sequenced_Reads_Num);#计算g_score
		next if ($G_Score<$g_score_threshold);#跳过阈值，不输出
		print OU "$class\t$Theoretical_Tag_Num\t$Sequenced_Tag_Num\t$Sequenced_Tag_Num2Theoretical_Tag_Num%\t";
		print OU "$Sequenced_Reads_Num\t$Sequenced_Reads_Num2Theoretical_Tag_Num\t$Sequenced_Reads_Num2Sequenced_Tag_Num\t$Sequenced_Tag_Num_2\t";
		print OU "$G_Score\n";
	}
	close OU;
	undef %hs_tag_num;
	undef %hs_detected_GCF_tag;
	print STDOUT "###($sample_name) Sample identification completed, ",`date`;
}
close LI;
		
# clean cache
print STDOUT "### Cleaning the cached objects started, ",`date`;
undef %hs_tag_theory_num;
undef %hs_GCF2class;
undef %hs_tag2GCF;
print STDOUT "###Cleaning the cached objects completed, ",`date`;

sub CheckDir{
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit;}
		}
		return 1;
}


