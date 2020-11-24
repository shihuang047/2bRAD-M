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

select STDOUT;$|=1;# cache cleaning

my($list,$site,$type,$outdir,$verbose,$help);
GetOptions(
		"l:s"  => \$list,
		"s:i"  => \$site,
		"t:s"  => \$type,
		"o:s"  => \$outdir,

		"c:i"  => \$cpu,
		"v:s"  => \$verbose,
		"r:s"  => \$remove_redundant,# yes or no, default value is "no"

		"h|help:s" => \$help,
		);

sub usage{# help
print STDERR "\e[;33;1m
DESCRIPTION
  It constructs the taxa-specific 2b-RAD reference genome database from a whole-genome reference database.
USAGE
  perl $0
PARAMETERS
  -l <str> genome classification list (the line which begins with # will be ignored)
           eg:GCF<tab>kingdom<tab>phylum<tab>class<tab>order<tab>family<tab>genus<tab>species<tab>strain<tab>genome_path
  -s <int> One or multiple type 2b restriction enzymes (sites).
           [1]CspCI  [9]BplI
           [2]AloI   [10]FalI
           [3]BsaXI  [11]Bsp24I
           [4]BaeI   [12]HaeIV
           [5]BcgI   [13]CjePI
           [6]CjeI   [14]Hin4I
           [7]PpiI   [15]AlfI
           [8]PsrI   [16]BslFI
  -t <str> The database level. One or more taxonomy level of the 2b-RAD reference database can be specified: kingdom,phylum,class,order,family,genus,species,strain. Use 'all' for any levels. (comma separated).
  -o <str> outdir (if not exists,it will be created)
OPTION
  -c <int> cpu [$cpu]
  -v <str> verbose
  -r <str> whether to delete redundant tags within the genome (yes or no) [$remove_redundant]
  -h|help  print this help
Author:  $author
Last update:  $time\e[0m\n";
}

if(defined($help)){
	&usage;
	exit 0;
}

unless($list && $site && $type && $outdir){
	&usage;
	print STDERR "para -l -s -t or -o error.\n";
	exit 1;
}

my %hs_site2enzyme=(# enzymes
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
			'species'  => '7',
			'strain'  => '8',
			);
my @head=('#Name','Kingdom','Phylum','Class','Order','Family','Genus','Species','Strain');

# check the parameter -t: specify the taxonomic level of 2b-RAD genome database
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
			exit 1;
		}
	}
}
# check the parameter -s:
unless(exists $hs_site2enzyme{$site}){
	&usage;
	print STDERR "-s parameter error: cannot find $site\n";
	exit 1;
}

#check the parameter -r: using default value "no"
unless($remove_redundant eq "yes" || $remove_redundant eq "no"){
	&usage;
	print STDERR "-r parameter error: cannot find $remove_redundant\n";
	exit 1;
}

&CheckDir($outdir);# create the output directory
&CheckDir("$outdir/$hs_site2enzyme{$site}");#within the output directory, create sub-directories for enzymes.

# digital digestion
&CheckDir("$outdir/$hs_site2enzyme{$site}/enzyme_result");# for each enzyme, create a directory for the fasta files with 2b-RAD reads by digital digestion from the whole genome

print STDOUT "###Electron digestion start, ",`date`;# output the log file
open LIST,"$list" or die "cannot open $list\n";
my $pm=new Parallel::ForkManager($cpu); # multi-threads computation
while(<LIST>){
	my $pid=$pm->start and next;
	my $line=$_;
	if(/^#/ || /^$/){;}else{# remove blank lines or lines starting with #
		chomp($line);
		my @tmp=split /\t/,$line;
		if(-e $tmp[-1]){# check the availability of a genome fasta file
			if(defined($verbose)){# verbose model
				print STDOUT "[Electron digestion] Analyze $tmp[0]\n";
			}
			system("perl $Bin/2bRADExtraction.pl -i $tmp[-1] -t 1 -s $site -od $outdir/$hs_site2enzyme{$site}/enzyme_result/ -op $tmp[0] -gz no 1> /dev/null; mv $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].$hs_site2enzyme{$site}.fa $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa");# digital digestion
		}else{
			print STDERR "[ERROR] $tmp[-1] does not exist,please check your genome file\n";
			exit 1;
		}
	}
	$pm->finish;
}
$pm->wait_all_children;
close LIST;
print STDOUT "###Electron digestion complete, ",`date`;# STDOUT

# Record the taxonomies of each 2b-RAD tag and identify taxa-specifc 2b-RAD tags
print STDOUT "###Record the taxonomies of each 2b-RAD tag and identification of taxa-specifc 2b-RAD tags -- start, ",`date`;#STDOUT
for my $level(sort {$hs_type{$a}<=>$hs_type{$b}} keys %hs_type){ #iterate all taxonomic levels of 2b-RAD database
	#Record the taxonomies of each 2b-RAD tag
	print STDOUT "###($level) Record the taxonomies of each 2b-RAD tag -- start, ",`date`;#STDOUT
	my %hash;
	open LIST,"$list" or die "cannot open $list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);# remove blank lines or lines starting with #
		my $line=$_;
		chomp($line);
		my @tmp=split /\t/,$line; # separate list items by \t
		my $class=join("\t",@tmp[1..$hs_type{$level}]); # concatenate the full taxonomic annotation
		if(defined($verbose)){# verbose model
			print STDOUT "[($level) Record tag] Analyze $tmp[0]\n";
		}
		open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n"; # open the 2b-RAD genome file by digital digestion
		while(<EN>){
			chomp;
			next if(/^>/);
			my $tag=$_;
			if(exists $hash{$tag}){# read in all 2b-RAD tags
				$hash{$tag}{$class}++; # record the corresponding taxonomies of each 2b-RAD tag
			}else{# check the reverse compliment
				$tag=~tr/ATCG/TAGC/;
				$tag=reverse($tag);
				$hash{$tag}{$class}++;
			}
		}
		close EN;
	}
	close LIST;
	print STDOUT "###($level) Record the taxonomies of each 2b-RAD tag -- complete, ",`date`;# STDOUT

	print STDOUT "###($level) Identification of taxa-specifc 2b-RAD tags -- start, ",`date`;# STDOUT
#	&CheckDir("$outdir/$hs_site2enzyme{$site}/database");
	open OU,"|gzip > $outdir/$hs_site2enzyme{$site}/$level.gz" or die "can not open $outdir/$hs_site2enzyme{$site}/$level.gz\n";
	print OU join("\t",@head[0..$hs_type{$level}]),"\tTags...\n"; # output headers including taxonomies of each 2b-RAD tag
	# compute the number of 2b-RAD tags specific to a taxon for a genome
    # output all taxa-specific 2b-RAD tags of a genome into the filepath specified
	open LIST,"$list" or die "cannot open $list\n";
	while(<LIST>){
		next if(/^#/ || /^$/);# remove blank lines and lines starting with #
		chomp;
		my @tmp=split /\t/;
		if(defined($verbose)){# verbose model
			print STDOUT "[(level) Identification label] Analyze $tmp[0]\n";
		}
		print OU join("\t",@tmp[0..$hs_type{$level}]); # concatenate the full taxonomy annotation
		open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n";
		while(<EN>){
			chomp;
			next if(/^>/);
			my $detection_tag=$_;# 2b-RAD tags to be examined
			unless(exists $hash{$detection_tag}){# check if a 2b-RAD tag exist in the hash table, otherwise it should exist in the hash table in the reverse-compliment format
				$detection_tag=~tr/ATCG/TAGC/;
				$detection_tag=reverse($detection_tag);
			}
			if(keys %{$hash{$detection_tag}}==1){ # examine if a 2b-RAD tag specific to a given taxon
				print OU "\t$detection_tag"; # output the 2b-RAD tag sequence into database file (last column)
			}
		}
		close EN;
		print OU "\n";
	}
	close LIST;
	close OU;
	undef %hash;
	print STDOUT "###($level) Identification of taxa-specifc 2b-RAD tags -- complete, ",`date`; #STDOUT
	# remove the redundant 2b-RAD markers within a given genome
    # The definition of a "redundant" 2bRAD marker: a 2b-RAD tag can be identified as a taxa-specific marker but occurs more than once in the original genome.
	if($remove_redundant eq "yes"){ # The default value is "no" not displaying to users in the help information
		print STDOUT "###($level) Remove redundant 2b-RAD markers -- start, ",`date`;# STDOUT
		open IN,"gzip -dc $outdir/$hs_site2enzyme{$site}/$level.gz|" or die "can not open $outdir/$hs_site2enzyme{$site}/$level.gz\n";# open the database file
		open OU,"|gzip > $outdir/$hs_site2enzyme{$site}/$level\_norebundancy.gz" or die "cannot open $outdir/$hs_site2enzyme{$site}/$level\_norebundancy.gz\n"; #output the de-redundant database file
		my $num_tag_col;
		while(<IN>){
			chomp;
			my $line=$_;
			my %hs_unique_tag;
			my @tmp=split /\t/,$line;
			if($.==1){ # line number ==1 or header
				$num_tag_col=$#tmp;
				print OU "$line\n";
				next;
			}
			for my $i($num_tag_col .. $#tmp){ # iterate columns for each line
				$hs_unique_tag{$tmp[$i]}++;
			}
			print OU join("\t",@tmp[0 .. $num_tag_col-1]);
			for my $i($num_tag_col .. $#tmp){ # check how many times a 2b-RAD tag show up within a genome
				if($hs_unique_tag{$tmp[$i]}==1){
					print OU "\t$tmp[$i]";
				}
			}
			print OU "\n";
			undef %hs_unique_tag;
		}
		close IN;
		close OU;
		print STDOUT "###($level) Remove redundant 2b-RAD markers -- complete, ",`date`;# STDOUT
	}
}
print STDOUT "###Record the taxonomies of each 2b-RAD tag and identification of taxa-specifc 2b-RAD tags -- complete, ",`date`;#STDOUT

print STDOUT "###The statistical summary of 2b-RAD markers -- start, ",`date`;# STDOUT
my %stat;
# compute the total number of 2b-RAD tags within a genome
open LIST,"$list" or die "cannot open $list\n"; # open each genome in the genome classification list
while(<LIST>){
	next if(/^#/ || /^$/);# remove blank lines and lines starting with #
	chomp;
	my @tmp=split /\t/;
	open EN,"$outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa" or die "cannot open $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa\n";
	while(<EN>){
		chomp;
		$stat{$tmp[0]}{"all"}++ if(/^>/); # the total number of 2b-RAD tags in the genome
	}
	close EN;
#	system("gzip -f $outdir/$hs_site2enzyme{$site}/enzyme_result/$tmp[0].fa"); # compress
}
close LIST;

# compute the total number of 2b-RAD taxa-specific markers within a genome
for my $i(keys %hs_type_database){ #iterate databases at different taxonomic levels
	next unless(-e "$outdir/$hs_site2enzyme{$site}/$i.gz");
	open IN,"gzip -dc $outdir/$hs_site2enzyme{$site}/$i.gz|" or die "cannot open $outdir/$hs_site2enzyme{$site}/$i.gz\n"; # open the file that have redundant 2b-RAD markers within a genome
	while(<IN>){
		chomp;
		my $line=$_;
		my @tmp=split /\t/,$line;
		next if($.==1 && $line=~/^#/);
		$stat{$tmp[0]}{$i}=$#tmp-$hs_type_database{$i}; # the total number of 2b-RAD taxa-specific markers in the genome
	}
	close IN;
}

# output
open OU,">$outdir/$hs_site2enzyme{$site}/stat.xls" or die "cannot open $outdir/$hs_site2enzyme{$site}/stat.xls\n";
print OU "#Unique_Name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\t";
print OU "All_Tag\tKingdom_Unique\tPhylum_Unique\tClass_Unique\tOrder_Unique\tFamily_Unique\tGenus_Unique\tSpecies_Unique\tStrain_Unique\n";
open LIST,"$list" or die "cannot open $list\n";
while(<LIST>){
	next if(/^#/ || /^$/);# remove blank lines and lines starting with #
	chomp;
	my @tmp=split /\t/;
	$stat{$tmp[0]}{"all"}=0 unless(exists $stat{$tmp[0]}{"all"});
	print OU join("\t",@tmp[0..$#tmp-1]),"\t",$stat{$tmp[0]}{"all"};# headers
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

#统计去冗余后的unique标签数
if($remove_redundant eq "yes"){
	print STDOUT "###Count the number of tags after deletion of redundant tags within the genome start, ",`date`;
	for my $GCF(keys %stat){#%stat变量 除基因组电子标签外，其他数据初始化
		for my $shuiping(keys %{$stat{$GCF}}){
			next if($shuiping eq "all");
			$stat{$GCF}{$shuiping}="-";
		}
	}

	#统计每个基因组各水平unique标签数
	for my $i(keys %hs_type_database){
		next unless(-e "$outdir/$hs_site2enzyme{$site}/$i\_norebundancy.gz");
		open IN,"gzip -dc $outdir/$hs_site2enzyme{$site}/$i\_norebundancy.gz|" or die "cannot open $outdir/$hs_site2enzyme{$site}/$i\_norebundancy.gz\n";
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
	open OU,">$outdir/$hs_site2enzyme{$site}/stat_norebundancy.xls" or die "cannot open $outdir/$hs_site2enzyme{$site}/stat_norebundancy.xls\n";
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
	print STDOUT "###Count the number of tags after deletion of redundant tags within the genome complete, ",`date`;
}
undef %stat;
print STDOUT "###The statistical summary of 2b-RAD markers -- complete, ",`date`;# STDOUT

print STDOUT "###Delete enzyme result start, ",`date`;# STDOUT
system("rm -rf $outdir/$hs_site2enzyme{$site}/enzyme_result/");# delete the 2b-RAD genome fasta files generated from the whole genome
print STDOUT "###Delete enzyme result complete, ",`date`;# STDOUT

print STDOUT "ALL DONE, ",`date`;

sub execute{# print the command line and execute
	my $cmd = shift;
	print STDOUT "$cmd\n";
	my $exit_code=system($cmd);
	if($exit_code!=0){
		print STDERR "Command $cmd failed with an exit code of $exit_code.\n";
		exit($exit_code >> 8);
	}
}
sub CheckDir{# create the directory
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
		}
		return 1;
}
