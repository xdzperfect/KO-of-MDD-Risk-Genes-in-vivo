#!/usr/bin/perl
use strict;
#use warnings;
my $Author="Xuzhzh";
my $writeTime="2018-1-5";
my $contact="xdzperfect@163.com";
my $version=3.00; 
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","f=s","o=s","h");		
if (!(defined $opts{i} &&defined $opts{o}) || defined $opts{h}) {
&usage;
}
my $in = $opts{'i'};
my $out = $opts{'o'};
my $type = defined $opts{'f'} ? $opts{'f'}:"\t";
sub usage{
print <<"USAGE";
Version: $version
Author: $Author
WriteTime: $writeTime
Contact: $contact
Usage:
$0 -i -o -f -h
options:
-i input file waitFor convertion Raw and Colum;
-o output file finished convertion newRaw and newColum;
-f the split formate, default Tab;
-h help
USAGE
exit(1);
}
#*************************************************************************************************************************
=pod
=head input file format(dataFram or Matirx file)
Key	A	B	C
K1	a1	b1	c1
K2	a2	b2	c2
K3	a3	b3	c3
����
=head output file format(convertedMatrix file)
Key	K1	K2	K3
A	a1	a2	a3
B	b1	b2	b3
C	c1	c2	c3
����
=cut
#***********************************************************************************************************************************
#�������ܣ��Ծ���������������������ݻ���С���ݶ����á�����С���ݿ���ֱ����Excel����ת�ã����߽������ݲ�ֳ�С����Ȼ����Excelת�á�
#***********************************************************************************************************************************
open A,"$in" or warn "open failed:$!";
my @row;
while(<A>)
{
	chomp;
	push @row,$_;	
}
close A;
my @newRow=split '\s',$row[0];                     #ȡ��һ�е�������Ϊ�µ�������#Key	A	B	C
open B,">$out" or warn "open failed:$!";
for(my $i=0;$i<=$#newRow;$i++)                     #����ȥÿ�ж�Ӧ����ֵ���кϲ���һ���ַ�����Ȼ�������
{	                                                 #��дÿ���ж�Ӧ�е��ַ������顣
	for(my $i2=1;$i2<=$#row;$i2++)
	{
		my @array=split '\s',$row[$i2];
		$newRow[$i]= $type eq "\t" ? $newRow[$i]."\t".$array[$i]: $newRow[$i]."$type".$array[$i];
		#$newRow[$i]=$newRow[$i]." ".$array[$i];
	}	                                               #��������������������Ӧ����ֵ��
	print B "$newRow[$i]\n" or die "can't write into the file";
	$newRow[$i]=0;  #updata in version2 for giving up the string, so that the mem will become smaller
}
close B;
