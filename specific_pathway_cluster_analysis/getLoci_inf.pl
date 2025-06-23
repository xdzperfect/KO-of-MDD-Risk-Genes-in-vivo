#!/usr/bin/perl
use strict;
#use warnings;
my $Author="Xuzhzh";
my $writeTime="2018-1-21";
#Any question please contact:xdzperfect@163.com;
my $version=2.00;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i1=s","i2=s","o=s","cf=s","h");
if (!(defined $opts{i1} &&$opts{i2} &&defined $opts{o}) || defined $opts{h}) {
&usage;
}
my $in1 = $opts{'i1'};
my $in2 = $opts{'i2'};
my $out = $opts{'o'};
my $convertFilter=defined $opts{'cf'} ? $opts{'cf'} : "no";
sub usage{
print <<"USAGE";
Version: $version
Author: $Author
WriteTime: $writeTime
Usage:
$0 -i1 -i2 -o -cf -h
options:
-i1 input file for loci1;
-i2 input file for include loc1 inf;
-o output file for get loci1 inf;
-cf filter give list,defualt is no;
-h help
USAGE
exit(1);
}
#*************************************************************************************************************************
=pod
=head input file format1(determined.fastq)
cg22461835
cg17689707s
=head input file format2(/share_bio/unisvx3/ciwm_group/xuzhzh/software/wgs_annotation/annovar/humandb2/hg19_cytoBandSplit/chr1.txt)
Keys TCGA-CH-5788-01.txt TCGA-EJ-5505-01.txt
cg13332474 0.898283202914206 0.167695921735896
cg17027195 0.0410937375033316 0.0368024965434248
cg22461835 0.131571403575738 0.37056097217209
cg04244851 0.854719963666444 0.923094049134669
cg19669385 0.944744945352528 0.755723542652723
cg04244855 0.91824429616199 0.91966346024429
cg17689707 0.0194723375302965 0.0932638904396283
cg02434381 0.0525068199121555 0.0875269828118246
cg05777492 0.0242416176183741 0.0324280629030931
▒▒
=head output file format
cg22461835 0.131571403575738 0.37056097217209
cg17689707 0.0194723375302965 0.0932638904396283
=cut
#***********************************************************************************************************************
#attached cytoband loci information if not coveraged
#version2: add parameter "cf", default is to get given list, if -cf yes, then script will filter the given list
#***********************************************************************************************************************
my (@key,%hash);
open A,"$in1" or warn "open failed:$!";
while(<A>)
{
        chomp;#p36.33  0       2300000 146     1       146
        my @array=split;
        push @key,$array[0];
}
close A;
for(my $i=0;$i<=$#key;$i++)
{
        $hash{$key[$i]}=0;
}
open B,"$in2" or warn "open failed:$!";
open C,">>$out" or warn "open failed:$!";
while(<B>)
{
	chomp;
	my @array=split;#p36.33  chr1    0       2300000 p36.33
	if($convertFilter eq "no")
	{
		if(exists $hash{$array[0]})
  	{
    	print C "$_\n" or die "can't write into the file";
  	}else
  	{
    	next;
  	}
  }else
  {
		if(exists $hash{$array[0]})
  	{
    	next;
  	}else
  	{
    	print C "$_\n" or die "can't write into the file";
  	}  
  }
}
close B;
close C;
