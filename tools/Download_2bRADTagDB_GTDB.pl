#!/usr/bin/perl
#Author:zhangrongchao, zhangrongchaoxx@163.com
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

$ARGV[0] ||="2B-RAD-M-ref_db_GTDB";

#if($#ARGV!=0){
#	print STDERR "perl $0 outdir\n";
#	exit 1;
#}

my $outdir=$ARGV[0];#下载目录

$outdir=abs_path($outdir);
&CheckDir("$outdir");


my @a=('abfh_classify','MSA1002','simulate_50');#分类表，实际数据，模拟数据
#my @b=('BcgI.species');#需要下载的库文件
my @b=('BcgI.species','CjePI.species');#需要下载的库文件

my %hash_path=(
	'abfh_classify'=>['https://figshare.com/ndownloader/files/31653170/abfh_classify_with_speciename.txt.gz',],

	'MSA1002'      =>['https://figshare.com/ndownloader/files/25623566/MSA1002_R1.fq.gz',],

#	'simulate_50'  =>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25621832/simulate_50.fa.gz',],
	'simulate_50'  =>['https://figshare.com/ndownloader/files/25915428/simulate_50.BcgI.fq.gz',],

	'BcgI.species' =>['https://figshare.com/ndownloader/files/31653911/BcgI.species.fa.gz0',
	                  'https://figshare.com/ndownloader/files/31659299/BcgI.species.fa.gz1',
	                  'https://figshare.com/ndownloader/files/31653614/BcgI.species.fa.gz2',],

	'CjePI.species'=>['https://figshare.com/ndownloader/files/31660241/CjePI.species.fa.gz0',
	                  'https://figshare.com/ndownloader/files/31660358/CjePI.species.fa.gz1',
	                  'https://figshare.com/ndownloader/files/31662320/CjePI.species.fa.gz2',
	                  'https://figshare.com/ndownloader/files/31662794/CjePI.species.fa.gz3',
	                  'https://figshare.com/ndownloader/files/31659818/CjePI.species.fa.gz4',],
	);

my %hash_md5=(
	'abfh_classify'=>['c2faa9ae97b704b3d0705709cf22ecb4',],

	'MSA1002'      =>['bc2b189213975f6d6c0833a4ba726239',],

#	'simulate_50'  =>['9defe990462d3fef8eb69a2c359d72da',],
	'simulate_50'  =>['04cafca5b5c23c48774e9d515dde42a8',],

	'BcgI.species' =>['a1b70d0de71093a0bb9bedbadab641b0',
	                  '383fd8c85a23aee4a48d48aa41845f17',
	                  'd19a5ce115fac8708fb0919f619ddf19',],
	
	'CjePI.species'=>['8b1c62c80bdf3b05182f2fe47d0f0751',
	                  '4662c85ef0e12a749d8b9284302e2a18',
	                  'ed3d3a27df05b7c0eb97140f78f54a75',
	                  '063b3c362f41889037b3bb15d8a0617f',
	                  '021a06a6e926b4ba91acba0c398877d7',]
	);

#合并后文件md5
my %complete_md5=(
	'BcgI.species' =>'eea6b5ec34b00a749d45199a91fd3e34',
	'CjePI.species'=>'3d9913da22ac340357d4e708a7506de8',
	);

#download abfh_classify && MSA1002 && simulate_50
for my $i(@a){
	my @tmp=split /\//,$hash_path{$i}[0];
	my $url=join("/",@tmp[0..$#tmp-1]);
	my $name=$tmp[-1];
	my $file_md5;#下载的文件的MD5值
	while(1){
		if(-e "$outdir/$name"){
			chomp($file_md5=`md5sum $outdir/$name`);
			$file_md5=(split /\s+/,$file_md5)[0];
		}
		if(-e "$outdir/$name" && $file_md5 eq $hash_md5{$i}[0]){
			print STDOUT "File $name has been downloaded.\n";
			last;
		}else{
			`wget -t 0 -O $outdir/$name $url`;
		}
	}
}
#example list
open OU,">$outdir/list_mock" or die "cannot open $outdir/list_mock\n";
print OU "MSA1002\t$outdir/MSA1002_R1.fq.gz\n";
close OU;

open OU,">$outdir/list_simulation" or die "cannot open $outdir/list_simulation\n";
#print OU "simulate_50\t$outdir/simulate_50.fa.gz\n";
print OU "simulate_50\t$outdir/simulate_50.BcgI.fq.gz\n";
close OU;



#下载数据库文件
for my $i(@b){
	my $cat="";
	while(1){
		my $md5;
		if(-e "$outdir/$i.fa.gz"){#存在完成文件
			chomp($md5=`md5sum $outdir/$i.fa.gz`);
			$md5=(split /\s+/,$md5)[0];
		}
		if(-e "$outdir/$i.fa.gz" && $md5 eq $complete_md5{$i}){
			print STDOUT "File $i.fa.gz hash been downloaded.\n";
			`rm -rf $cat`;
			last;
		}else{
			for my $j(0..$#{$hash_path{$i}}){#循环每个文件
				my @tmp=split /\//,$hash_path{$i}[$j];
				my $url=join("/",@tmp[0..$#tmp-1]);
				my $name=$tmp[-1];
				my $file_md5;#下载的文件的MD5值
				while(1){
					if(-e "$outdir/$name"){
						chomp($file_md5=`md5sum $outdir/$name`);
						$file_md5=(split /\s+/,$file_md5)[0];
					}
					if(-e "$outdir/$name" && $file_md5 eq $hash_md5{$i}[$j]){
						print STDOUT "File $name has been downloaded.\n";
						$cat .=" $outdir/$name";
						last;
					}else{
						`wget -t 0 -O $outdir/$name $url`;
					}
				}
			}
			`cat $cat > $outdir/$i.fa.gz`;
		}
	}
}

print STDOUT "Congratulations! All databases have been downloaded.\n";

sub CheckDir{
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}






