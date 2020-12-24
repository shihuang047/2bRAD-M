#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

my $author="Zheng Sun, Rongchao Zhang, Shi Huang";
my $time="2020.12.21";

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
    -l <str>    The path of the input filepath list (the line that begin with # will be ignored) e.g: sample_name<tab>data_path(fa|fq)(.gz).
    -d <str>    The database filepath.
    -t <str>    The taxonomy level of the taxa-specific 2b-RAD database used. It should be one of the following: kingdom,phylum,class,order,family,genus,species,strain.
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
    -g <int>    The threshold of G score [$g_score_threshold, it means >=$g_score_threshold]. To control the false-positive in the species identification, G score was derived for each speciesidentified within a sample, which is a harmonious mean of read coverage of 2b-RAD markers belongs to a species and number of all possible 2b-RAD markers of this species. Therecommended/default threshold is $g_score_threshold.
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

my @HEAD=(
	'Kingdom',
	'Phylum',
	'Class',
	'Order',
	'Family',
	'Genus',
	'Species',
	'Strain',
	);

unless($list && $database && $level && $site && $outdir){
	&usage;
	exit 1;
}

#转换绝对路径
$list=abs_path($list);
$database=abs_path($database);
$outdir=abs_path($outdir);


# parameter checking
unless($verbose eq "yes" || $verbose eq "no"){
	&usage;
	print STDERR "Parameter -v is wrong\n";
	exit 1;
}
# check the taxonomic level of a 2b-RAD reference genome database
unless($level eq "kingdom" || $level eq "phylum" || $level eq "class" || $level eq "order" || $level eq "family" || $level eq "genus" || $level eq "species" || $level eq "strain"){
	&usage;
	print STDERR "Parameter -t is wrong. Cannot get $level\n";
	exit 1;
}
# check the parameter -s and -d
unless(exists $hs_site2enzyme{$site}){
	&usage;
	print STDERR "Parameter -s $site is wrong\n";
	exit 1;
}
#检查库文件
unless(-e "$database/$hs_site2enzyme{$site}.$level.fa.gz" && -e "$database/abfh_classify_with_speciename.txt.gz"){
	&usage;
	print STDERR "Incomplete database, please check the parameter(-d).\n";
	exit 1;
}

print STDOUT "COMMAND: perl $0 -l $list -d $database -t $level -s $site -o $outdir -g $g_score_threshold -v $verbose\n";

&CheckDir($outdir);

my $head=join("\t",@HEAD[0..$hs_type_database{$level}-1]);
# load the database
print STDOUT "### Loading the database, $database/$hs_site2enzyme{$site}.$level.fa.gz, ",`date`;
my (%hs_tag2GCF,%hs_GCF2class,%hs_tag_theory_num);
my $all_genome_num;
open LI,"gzip -dc $database/abfh_classify_with_speciename.txt.gz|" or die "cannot open $database/abfh_classify_with_speciename.txt.gz\n";
while(<LI>){
	next if(/^#/ || /^$/);#去掉注释行和空行
	chomp;
	my @tmp=split /\t/;
	my $class=join("\t",@tmp[1..$hs_type_database{$level}]);# all taxonomic levels
	$hs_GCF2class{$tmp[0]}=$class;# record the corresponding taxonomy for each GCF
	$all_genome_num++;
}
close LI;

my (%hash_gcf_rank,%complete);
$/=">";
open IN,"gzip -dc $database/$hs_site2enzyme{$site}.$level.fa.gz|" or die "cannot open $database/$hs_site2enzyme{$site}.$level.fa.gz\n";
<IN>;
while(<IN>){
	chomp;
	my($id,$tag)=split /\n/;
	#GCF号|基因组内部标签排序|scaffoldid|startpos|正反向酶切|是否为指定水平下unique&&noredundancy标签
	my @tmp=split /\|/,$id;
	$hash_gcf_rank{$tmp[0]}++;
	next if($tmp[5]!=1);#跳过非unique标签
	my $class=$hs_GCF2class{$tmp[0]};# all taxonomic levels
	push @{$hs_tag2GCF{$tag}},$tmp[0];# record GCF for each 2b tag
	$hs_tag_theory_num{$class}{$tmp[0]}{$tag}++;# compute the number of 2b tags of each GCF under a given taxon and record the # of all 2b tags from the same taxa
	for (my $i=100;$i>0;$i=$i-1){ #每完成1%则输出进度
		if((keys %hash_gcf_rank)/$all_genome_num*100>=$i){
			print STDOUT "$i% " unless(exists $complete{$i});#仅没输出过的才会输出日志
			$complete{$i}++;
			last;
		}
	}
}
close IN;
$/="\n";
print STDOUT "###Loading database completed, ",`date`;

# process each sample in the list file
open LI,"$list" or die "cannot open $list\n";
while(<LI>){
	next if(/^#/ || /^$/); # 去除注释行和空行
	chomp;
	my ($sample_name,$sample_data)=split /\t/;
	$sample_data=abs_path($sample_data);#转为绝对路径
	print STDOUT "###($sample_name) Sample identification started, ",`date`;
	my (%hs_tag_num,%hs_detected_GCF_tag);
	# load a single sample
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
			$tag=~tr/ATCG/TAGC/;
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

	if((keys %hs_tag_num)==0){# go to the next sample if no 2b-RAD tag was detected in a sample
		print STDERR "!!!($sample_name) Warning: $hs_site2enzyme{$site}  $level the number of 2b-RAD tags for this sample is zero\n";
		print STDOUT "###($sample_name) Sample idenfication completed, ",`date`;
		next;
	}
	&CheckDir("$outdir/$sample_name");# create a filepath for each sample
	# compute the number of therotical and actual 2b-RAD tags for each GCF in each sample
	open DE,">$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.GCF_detected.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}.GCF_detected.xls\n";
	for my $class(sort {$a cmp $b} keys %hs_detected_GCF_tag){
		for my $GCF(sort {$a cmp $b} keys %{$hs_detected_GCF_tag{$class}}){
			my $GCF_all_theory_num;
			my $detected_tag_num;
			$GCF_all_theory_num=keys %{$hs_tag_theory_num{$class}{$GCF}};# theoratical, the number of taxa-specific 2b-RAD tags from each GCF
			$detected_tag_num=keys %{$hs_detected_GCF_tag{$class}{$GCF}};# in a real sample, the number of taxa-specific 2b-RAD tags from each GCF
			my $percent=sprintf "%.4f",$detected_tag_num/$GCF_all_theory_num;# the percentage of detected 2b-RAD tags from a GCF specific to a taxon
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
		if($verbose eq "yes"){ # output the detailed information on sequencing coverage of each 2b-RAD tag from a given genome detected in the real sample
			my $output_name=(split /\t/,$class)[-1];
			open DETAIL,">$outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}/$output_name.xls" or die "cannot open $outdir/$sample_name/$sample_name.$hs_site2enzyme{$site}/$output_name.xls\n";
		}
		my ($Theoretical_Tag_Num,$Sequenced_Tag_Num,$Sequenced_Tag_Num_2,$Sequenced_Tag_Num2Theoretical_Tag_Num);
		my ($Sequenced_Reads_Num,$Sequenced_Reads_Num2Theoretical_Tag_Num,$Sequenced_Reads_Num2Sequenced_Tag_Num);
		my ($G_Score);
		$Sequenced_Tag_Num=$Sequenced_Reads_Num=$Sequenced_Tag_Num_2=0;
		for my $tag(keys %{$hs_tag_num{$class}}){# iterate each 2b-RAD tag
			$Sequenced_Tag_Num++;
			$Sequenced_Tag_Num_2++ if($hs_tag_num{$class}{$tag}>1);# compute the number of 2b-RAD tags that have the sequencing coverage >1
			$Sequenced_Reads_Num+=$hs_tag_num{$class}{$tag}; # the number of reads detected
			if($verbose eq "yes"){
				print DETAIL "$tag\t$hs_tag_num{$class}{$tag}\n";
			}
		}
		if($verbose eq "yes"){
			close DETAIL;
		}
		# average number of theoretical 2b-RAD tags for each taxon
		my $species_all_theory_num;
		for my $GCF(keys %{$hs_tag_theory_num{$class}}){
			for my $tag(keys %{$hs_tag_theory_num{$class}{$GCF}}){
				$species_all_theory_num+=$hs_tag_theory_num{$class}{$GCF}{$tag};
			}
		}
		$Theoretical_Tag_Num=$species_all_theory_num/(keys %{$hs_tag_theory_num{$class}});# average number of theoretical 2b-RAD tags for each taxon
        # statistical summmary
		$Sequenced_Tag_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Tag_Num/$Theoretical_Tag_Num*100;# 测到的标签占理论的百分比
		$Sequenced_Reads_Num2Theoretical_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Theoretical_Tag_Num;# 测到的标签深度/理论标签数
		$Sequenced_Reads_Num2Sequenced_Tag_Num=sprintf "%.8f",$Sequenced_Reads_Num/$Sequenced_Tag_Num;# 测到的标签平均深度
		$G_Score=sprintf "%.8f",sqrt($Sequenced_Tag_Num*$Sequenced_Reads_Num);#compute the g_score for each taxon
		next if ($G_Score<$g_score_threshold);# filter taxa that have g_score < $g_score_threshold
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
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
		}
		return 1;
}


