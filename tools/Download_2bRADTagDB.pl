#!/usr/bin/perl
#Author:zhangrongchao, zhangrongchaoxx@163.com
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

$ARGV[0] ||="2B-RAD-M-ref_db";

#if($#ARGV!=0){
#	print STDERR "perl $0 outdir\n";
#	exit 1;
#}

my $outdir=$ARGV[0];#下载目录

$outdir=abs_path($outdir);
&CheckDir("$outdir");


my @a=('abfh_classify','MSA1002','simulate_50');#分类表，实际数据，模拟数据
my @b=('BcgI.species');#需要下载的库文件
#my @b=('BcgI.species','CjePI.species');#需要下载的库文件

my %hash_path=(
	'abfh_classify'=>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25889157/abfh_classify_with_speciename.txt.gz',],

	'MSA1002'      =>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25623566/MSA1002_R1.fq.gz',],

#	'simulate_50'  =>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25621832/simulate_50.fa.gz',],
	'simulate_50'  =>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25915428/simulate_50.BcgI.fq.gz',],

	'BcgI.species' =>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25889544/BcgI.species.fa.gz0',
	                  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25889658/BcgI.species.fa.gz1',],

	'CjePI.species'=>['https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25891653/CjePI.species.fa.gz0',
	                  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25890987/CjePI.species.fa.gz1',
	                  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25890996/CjePI.species.fa.gz2',
	                  'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/25891002/CjePI.species.fa.gz3',],
	);

my %hash_md5=(
	'abfh_classify'=>['25f3a20babb56fd9f2a61eeddb82151a',],

	'MSA1002'      =>['bc2b189213975f6d6c0833a4ba726239',],

#	'simulate_50'  =>['9defe990462d3fef8eb69a2c359d72da',],
	'simulate_50'  =>['04cafca5b5c23c48774e9d515dde42a8',],

	'BcgI.species' =>['b36cc8e85fb68f1b3cc5301c49cafe98',
	                  '071b711730ce87e6c1f85f29319a5979',],
	
	'CjePI.species'=>['a32c1998d0d800fe336d9f03756b8409',
	                  '1eb528474f89a6550f69c160d0885dd8',
	                  'b803ea6b0e2bca1c6381b2a15a76876d',
	                  'a3d5f018fb3410b507759f2eabee4d04',]
	);

#合并后文件md5
my %complete_md5=(
	'BcgI.species' =>'75171aabcb754e827e5824ae755d06af',
	'CjePI.species'=>'bcfdef3722dfc763e09fd185f580198d',
	);

#download abfh_classify && MSA1002 && simulate_50
for my $i(@a){
	my $name=(split /\//,$hash_path{$i}[0])[-1];
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
			`wget -t 0 -O $outdir/$name $hash_path{$i}[0]`;
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
				my $name=(split /\//,$hash_path{$i}[$j])[-1];
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
						`wget -t 0 -O $outdir/$name $hash_path{$i}[$j]`;
					}
				}
			}
			`cat $cat > $outdir/$i.fa.gz`;
		}
	}
}



sub CheckDir{
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}






