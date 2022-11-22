#!/usr/bin/env perl
use strict;
#use warnings;
use File::Basename;
use Getopt::Std;
use POSIX; #计算时间

####################################################
# 
# Usage of run_dat2fa_swiss_prot.pl # 2022.11.22 update
# 
####################################################
sub usage{
    die(
        qq!
Usage:     perl $0 -i <dat file(required)> -o <fasta file>
Function:  to convert dat file from uniprot(Swiss-Prot) to fasta file.
Command:   -i input file (dat file) from uniprot(Swiss-Prot)(required). A gzip file with '.gz' suffix is also acceptable.
           -o output file with fasta format. default: related to input file. 
           -h get the help tips.
Author:    Shimao Zheng, zhengshimao007@163.com
version:   v0.2
Update:    2022.11.22
Notes:     Follow my WeChat public account 'The_Elder_Student'
\n! 
    )
}
#&usage; #调用帮助测试
my $dat_file;
my $fasta_file; 
my $fasta_info; #2022.11.22 add
####################################################
#
# 命令行参数的定义和获取，记录程序初始时间，设置参数默认值.
# Get the parameter and provide the usage.
#
####################################################
my %opts;
getopts( 'i:o:h', \%opts );
&usage && exit 1 unless ( exists $opts{i} );
&usage if ( exists $opts{h} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#output file
my ($name, $path, $suff)=fileparse("$opts{i}", qr/\.\w+$/);

$dat_file = $opts{i};

#不知怎么把输出文件处理部分写重复了 #2022.11.22 update
if($suff eq "\.gz"){
    if($opts{o}){
        $fasta_file = $opts{o};
        $fasta_info = "$opts{o}.info";
    }else{
        $fasta_file = "$name.fasta";
        $fasta_info = "$name.info";
    }
}else{
    #$fasta_file = $opts{o}||"$opts{i}.fasta";
    if($opts{o}){
        $fasta_file = $opts{o};
        $fasta_info = "$opts{o}.info";
    }else{
        $fasta_file = "$opts{i}.fasta";
        $fasta_info = "$opts{i}.info";
    } 
}

die "\nERROR: $fasta_file exist. Give another file name please!\n" if (-e $fasta_file);
# die "\nERROR: $fasta_info exist. Give another file name please!\n" if (-e $fasta_info); # 2022.11.22 add

print "Input file is $dat_file\nOutput file is $fasta_file , $fasta_info\n\n"; # 2022.11.22 update


####################################################
#
# To convert dat file from uniprot(Swiss-Prot)  to fasta file.
#
####################################################
if($suff eq "\.gz"){
    open DAT,"zcat $opts{i} | " || die "ERROR:can't open $opts{i} !";
}else{
    open DAT,"<","$opts{i}" || die "ERROR:can't open $opts{i} !";
}

open FA,">","$fasta_file" || die "ERROR:can't open $fasta_file!";
#number 44064;
#number 44180; #2022.10.22

open INFO,">","$fasta_info" || die "ERROR:can't open $fasta_info!"; #2022.11.22 add

    my $id; #UniProt ID
    my $ac; #accession
    my $sv; #sequence version
    my $full; #full name;
    my $os; #Organism
    my $ox; #Organism Taxonomy
    my $pe; #Protein Existence
    my $sq; #Sequence

$/="//\n"; # don't use '//',because http links are included in dat files;
print INFO "AC\tID\tFULL\tOS\tOX\tGN\tPE\tSV\n"; #2022.11.22 add
while(<DAT>){

    my $gn; #gene name #必须命名在循环内，因为部分序列无GN。如果定义不在循环内，则会造成下一个模块内GN未被重置。
   
    chomp;
    #print ;

    #if(m/^ID\s+(\S+)\s+/){ $id = $1 }; #first method
    $id = $1 if(m/^ID\s+(\S+)\s+/); #second method

    $ac = $1 if(m/AC\s+(\w+);/); #Accession #部分含有多个ID，这里只选了第一个。
    #print "$ac\n";

    $sv = $1 if( m/DT.*sequence version\s+(\d+)/); #why can't use "^DT"?;
    #if(/DT.*sequence version\s+(\d+)/){print "$1\n"};

    $full = $1 if(m/DE\s+RecName: Full=(.*);/); #匹配的第一个DE，后面的DE有同样格式的，被忽略了。
    #去掉可能含有的{ECO.*}
    $full =~ s/\s+\{.*\}//g;#{}从5.22开始，正则中的{}要加转义符 #报错解决https://www.codercto.com/a/92293.html

    $os = $1 if(m/OS\s+(.*)\./);
    $os =~ s/\s+\(.*\)//; #NOTE:'\'; #去掉物种名后面可能的其它信息

    $ox = $1 if(m/OX\s+NCBI_TaxID=(\d+)/); #数字之后不一定是 ';'，如：11S1_CARIL

    $gn = $1 if(m/GN\s+Name=(.*?);/); #部分序列无GN,如108_SOLLC; 部分序列还有多个GN; 部分GN含有多个 ';'; 
    #print "$gn\n"; #报错无法解决，去掉 'use warnings;' 才可以运行。但是部分可能为空值。
    #尝试写GN是否存在，如果不存在，$gn重置为空

    $pe = $1 if(m/PE\s+(\d+):/);

    #print "$id\$ac\t$sv\t$full\t$os\t$ox\t$gn\t$pe\n";

    # 打印fasta格式的第一行
#=pod
    if($gn){
        print FA ">sp|$ac|$id $full OS=$os OX=$ox GN=$gn PE=$pe SV=$sv";
    }else{
        print FA ">sp|$ac|$id $full OS=$os OX=$ox PE=$pe SV=$sv";
    }
#=cut
    print INFO "$ac\t$id\t$full\t$os\t$ox\t$gn\t$pe\t$sv\n"; # fasta_info文件

    if(m/SQ   SEQUENCE\s+.*;/){  #学习捕获中的$`, $&, $'的含义。
=pod
    $`表示匹配起始位置之前的字符串
    $&表示匹配的内容，即//内的内容
    $'表示匹配终结位置之后的内容
=cut 
        #$sq = $';
        #$sq =~ s/\s+//g; #记得加/g 进行全局替换。
        #赋值加替换，单行写法
        ($sq = $') =~ s/ //g;#记得加/g 进行全局替换，删除空格。 
        print FA "$sq";
    }
}
close DAT;
close FA;
close INFO; #2022.11.22 add
####################################################
#
# Record the program running time!
# 输出程序运行时间
#
####################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\n";
print "Follow my WeChat public account 'The_Elder_Student'\n";
print "OK!\n";
