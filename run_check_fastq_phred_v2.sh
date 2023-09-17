#!/usr/bin/bash

## ---------- 定义脚本的使用方法 ----------------------------------
usage() {
cat << EOF
Usage: $0 -i <fastq.gz>
	ls *.fastq.gz |xargs -i bash $0 -i {}

Options:
	-i, --input Input file.
	-o, --ouput Output file.
	-n, --read_nums Number of reads to test [default: 2500].
	-f, --force Force to rerun it.
	-h, --help Print the tips.
EOF
    exit 1
}

## ---------- 脚本选项设置 ----------------------------------
OPTIONS=$(getopt -o hi:o::n:: --long help,input:,output::,read_nums:: -n "ERROR:$0" -- "$@")

#if [ $? != 0]; then 
#	echo "Terminating..." >&2 ; 
#	exit 1; 
#fi

if [ $? -ne 0 ]; then 
	usage && exit 1
fi

# Note the quotes around `$OPTIONS': they are essential!
eval set -- "$OPTIONS"

## ---------- 设置默认值 ----------------------------------
FILE1=""
FILE2=""
READ_NUMS=2500
FORCE=false
## ---------- 处理解析后的命令行参数 ----------------------------------
while true; do
  case "$1" in
    -h|--help)
      usage
      ;;
	-f|--force)
	  FORCE=true
	  shift
	  ;;
    -i|--input)
      FILE1="$2" # 赋值
      shift 2
      ;;
	-o|--output)
      FILE2="$2" # 赋值
      shift 2
      ;;
	-n|--read_nums)
      READ_NUMS="$2" # 赋值
      shift 2
      ;;  
    --)
      shift
      break
      ;;
    *)
      echo "Cannot recognize: $1"
      usage
      ;;
  esac
done

## ---------- 流程接续用文件 ----------------------------------
<<EOF
this_scirpt=$(readlink -f "$0")
#echo "$this_scirpt" # 测试脚本所在位置。
end_log="${this_scirpt}.end_log"
EOF

## ---------- 参数的处理 ----------------------------------
num_lines=$((READ_NUMS * 4))


## ---------- 定义函数 ----------------------------------
### 处理压缩文件
function get_fq_lines(){ # $1 fastq.gz file; $2 number of lines of fastq

	if [[ "$1" =~ \.gz$ ]]; then
		zcat "$1" |head -n "$2"
	else
		cat "$1" |head -n "$2"
	fi
}

#get_fq_lines "$FILE1" "$num_lines" # 简单测试下这个函数

### 检查fastq文件格式
function check_fq (){ # $1 fastq.gz file; $2 number of lines of fastq
	file=$(basename $1) # print filename only
	#zcat  $1 |\
	#head -n 10000 | \
	get_fq_lines $1 $2 | \
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
                     print file"\tSolexa+64"; else print file"\tUnknown score encoding!";}'
}

## ---------- 检查与运行 ----------------------------------
## 检查输入文件是否存在
if [ ! -f "$FILE1" ]; then
    echo "ERROR: $i not exits."
	usage
    exit 1
fi

#echo "---------- step 1 start: `date`----------------------------------"
check_fq $FILE1 $num_lines  # check phred score
#echo "---------- step 1 done: `date`----------------------------------"

## ---------- 流程接续 ----------------------------------
<<EOF

if [[ $FORCE ]] ;then rm -rf $end_log; fi

#[[ -f $end_log ]] || ( echo "$FILE1" && touch ${end_log} )

if [[ ! -f $end_log ]] ;then
	for i in $(ls $FILE1); do
		check_fq $i $num_lines
	done         # 主体命令
fi

if [ $? -ne 0 ]; then
		usage && exit 1
	else
		touch ${end_log}
fi
EOF

#echo "---------- DONE: `date`----------------------------------"