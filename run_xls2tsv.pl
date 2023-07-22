#!/usr/bin/perl -w
 
use strict;
use Getopt::Long;
use File::Basename;
use Spreadsheet::ParseExcel; # 查看版本 perl -MSpreadsheet::ParseExcel -e  "print Spreadsheet::ParseExcel->VERSION"
use Spreadsheet::ParseExcel::FmtUnicode; # Unicode::Map


####################################################
# 
# Usage of run_xls2tsv.pl
# 
####################################################

sub USAGE {
	my $usage=<<"USAGE";
Usage:     perl $0  
Function:  To stat results of RNA seq.
Command:   --xls	-x	Excel file.
           --out	-o	Output file.
           --help	-h	Get the help tips.
Author:    Shimao Zheng, zhengshimao007\@163.com
version:   v0.1
Create:    2023.07.22
Update:    -
USAGE
	print $usage;
	exit;
}

# &USAGE; # 测试调用打印帮助信息

####################################################
# 
# Default parameters⭐
# 
####################################################

my $xls;
my $output_txt;

# 获取输出文件名
#my $filename = basename($xls, ".xls");
#my $dir = dirname($xls);
#my $output_txt = $dir."/".$filename.".txt";

####################################################
# 
# GetOptions
# 
####################################################
#Getopt::Long::Configure("bundling");

GetOptions(
				"help|h" => \&USAGE,
				"xls|x:s" => \$xls,
                "out|o:s" => \$output_txt
				) or &USAGE;
&USAGE unless ($xls && $output_txt); # 使用--help，-h或者在不设置$sample_list时均可调用&USAGE


# 使其支持中文
my $code = "CP936"; 
my $oFmtJ = Spreadsheet::ParseExcel::FmtUnicode->new(Unicode_Map => $code);

my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse($xls,$oFmtJ) or die "Can\'t open $xls $!";

=pod
if ( !defined $workbook ) {
    die $parser->error(), ".\n";
}
=cut

my $value;
my $txt;
my $worksheet = $workbook -> worksheet("Sheet1");

    my ( $row_min, $row_max ) = $worksheet->row_range();
    my ( $col_min, $col_max ) = $worksheet->col_range();
 
    for my $row ( $row_min .. $row_max ) {
        for my $col ( $col_min .. $col_max ) {
 
            my $cell = $worksheet->get_cell( $row, $col );
=pod            
            $value = "" unless $cell;
            $value = $cell -> value();
=cut
            unless($cell){
                $value = ""; # 未定义值记为"";

            }else{
                $value = $cell -> value();
            }

            $txt .= $value;
            if($col != $col_max){
                $txt .= "\t";
            }else{
                $txt .= "\n";
            }
        }
    }
    #print $row_min, "\n",$row_max, "\n",$col_min,"\n",$col_max,"\n";

open (OUT, ">$output_txt");
# print $txt;
print OUT  $txt;
close OUT;
