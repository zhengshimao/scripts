#!/usr/bin/bash

function usage(){
	echo "USAGE: sh $0 *.fastq.gz"
}

if [ $# -eq 0 ];then
	usage
fi

function check_fq (){
	file=$(basename $1) # print filename only
	zcat  $1 |\
 head -n 10000 | \
  awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
  awk -v file=$file 'BEGIN{min=100;max=0;} \
      {for(i=1;i<=NF;i++) \
          {if($i>max) max=$i; \
               if($i<min) min=$i;}}END \
          {if(max<=74 && min<59) \
                     print file"\tPhred+33"; \
           else \
           if(max>73 && min>=64) \
                     print file"\tPhred+64"; \
           else \
           if(min>=59 && min<64 && max>73) \
                     print file"\tSolexa+64"; else print "Unknown score encoding!";}'
}

for i in $@;do
	check_fq $i
done
