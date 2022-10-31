#!/usr/bin/env bash

#检查$1 # $#表示参数格式  $?表示上一个命令的返回状态
if [ $# -eq 0 ];then
	echo -e "Usage:bash $0 PRJNA881764\nplease give a parameter of project number,e.g. PRJNA275632";
	exit 1;
fi
echo $1 |grep -q 'PRJNA'
if [ $? -ne 0 ];then
	echo -e "Usage:bash $0 PRJNA881764\nplease give a parameter of project number,e.g. PRJNA275632";
	exit 1;
fi

mkdir -p ena_info_sample && cd ena_info_sample

project=$1 #PRJNA275632
url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$project&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_md5,fastq_ftp,fastq_aspera,submitted_ftp,sra_md5,sra_ftp,sample_title&format=tsv&download=true&limit=0"
if [ -f progect_infor.txt ];then
	next
else
	wget  "$url" -O progect_infor.txt -q;
fi

if [ $? -eq 0 ];then
    echo "Successfully! progect_infor.txt downloaded!";
else
    echo "ERROR:Download of the progect_infor.txt failed!"
    exit 1;
fi

# rename sample name #处理空格部分略显啰嗦了2022.10.7
sed '1d' progect_infor.txt| awk -F "\t" '{print "rename s/"$4"/"$13"/g "$4"*.gz"}' |sed -e "s/rename s/rename 's/g" -e "s/\/g/\/g'/g" -e 's/-/_/g' -e 's/ /_/g' |sed -e 's/rename_/rename /g' -e 's/_SRR/ SRR/g' > run_rename_sample.sh

# get md5 of fastq
cat progect_infor.txt |awk -F"\t"  '{print $7}' |sed 's/;/\n/g' > fq_md5
cat progect_infor.txt |awk -F"\t"  '{print $8}' |sed 's/;/\n/g'  | sed 's/ftp.*SRR/SRR/g' > fq_file
paste -d "  " fq_md5 fq_file |sed '1d'> fastq_md5.txt
rm fq_file fq_md5

# run_downloaded_fastq_using_aspera.sh
cat progect_infor.txt |awk -F "\t" '{print $9 }'| sed -e '1d' -e 's/;/\n/g' |awk '{print "ascp  -vQT -l 500m -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@"$1" ./"}' > run_downloaded_fastq_using_aspera.sh
# example
# ascp  -vQT -l 500m -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR122/079/SRR12207279/SRR12207279_1.fastq.gz  ./

# run_downloaded_sra_using_wget|axel.sh

cat progect_infor.txt |sed '1d'|awk -F "\t" '{print "wget -c -O "$4".sra "$12 }' > run_wget_sra_from_ena.sh
cat progect_infor.txt |sed '1d'|awk -F "\t" '{print "axel -n10 -o "$4".sra "$12 }' > run_axel_sra_from_ena.sh

# get md5 of sra
cat progect_infor.txt |sed '1d'|awk -F "\t" '{print $11"  "$4".sra"}' >sra_md5.txt

# fasterq-dump + pigz #首选
cat progect_infor.txt |sed '1d'|awk -F "\t" '{print "fasterq-dump --threads 8 --split-3 -O ./ "$4".sra;pigz -p 8 "$4"*.fastq"}' > run_pigz_fasterq_dump.sh #328837
# parallel-fastq-dump不加--gzip是未压缩的fastq文件,太占空间，舍
#cat progect_infor.txt |sed '1d'|awk -F "\t" '{print "parallel-fastq-dump --threads 8 --outdir ./ --split-3 -s "$4".sra"}' > parallel-fastq-dump.sh #280122
# parallel-fastq-dump的--gzip模式很慢
cat progect_infor.txt |sed '1d'|awk -F "\t" '{print "parallel-fastq-dump --threads 8 --outdir ./ --split-3  --gzip -s "$4".sra"}' > run_gzip_parallel-fastq-dump.sh #282962

# 赋予所有bash脚本可执行权限
chmod u+x *.sh
# 回到原文件夹
cd ..
# 赋予该bash脚本权限
chmod u+x $0


# 常用Kingfisher下载
echo "###################################################"
echo -e "#The 1st way"
echo "conda activate kingfisher"
echo "kingfisher get -p $1 -m ena-ftp aws-http aws-cp prefetch --download-threads 8 -f fastq.gz"
echo "md5sum -c fastq_md5.txt > md5.res"
echo "#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩"
echo "bash run_rename_sample.sh"
echo
# ascp从ENA下载fastq文件
echo "###################################################"
echo -e "#The 2nd way"
echo "ParaFly -c run_downloaded_fastq_using_aspera.sh -CPU 5 -failed_cmds failed_run_downloaded_fastq_using_aspera.txt"
echo "md5sum -c fastq_md5.txt > md5.res"
echo "#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩"
echo "bash run_rename_sample.sh"
echo
# ascp从ENA下载sra文件 #使用axel或者wget
echo "###################################################"
echo -e "#The 3rd way"
echo "ParaFly -c run_wget_sra_from_ena.sh -CPU 10 -failed_cmds failed_wget_sra_from_ena.txt #wget下载sra"
echo "ParaFly -c run_axel_sra_from_ena.sh -CPU 1 -failed_cmds failed_axel_sra_from_ena.txt #axel下载sra #设置多CPU可能报错,原因未知"
echo "md5sum -c sra_md5.txt > md5.res #校验"
echo "bash rename_sample.sh"
echo "ParaFly -c run_pigz_fasterq_dump.sh -CPU 2 -failed_cmds failed_pigz_fasterq_dump.txt #首选"
parallel-fastq-dump -h 1>/dev/null 2>&1
if [ $? -eq 0 ];then
    echo 'ParaFly -c run_gzip_parallel-fastq-dump.sh -CPU 5 -failed_cmds failed_gzip_parallel-fastq-dump.txt #次选'
else
    echo "conda install -c bioconda parallel-fastq-dump"
    echo 'ParaFly -c run_gzip_parallel-fastq-dump.sh -CPU 5 -failed_cmds failed_gzip_parallel-fastq-dump #次选'
fi
echo "#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩"
echo "bash run_rename_sample.sh"
exit 0
