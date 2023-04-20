# scripts

 我的常用脚本。记录下小工具们。能运行不代表一定正确，脚本都不难，有语言基础的话建议看下其中代码。

使用建议：将脚本直接放在已经加入环境变量的路径下`~/scripts` 

博客[《学生信的大叔》](https://zhengshimao.netlify.app/)

公众号：《学生信的大叔》

## 脚本介绍

### rd

`rd` == `real directory`

这是一个小`bash脚本`。

根据文件(夹)的路径获取其绝对路径。支持`多个参数`输入，支持`管道符`。

详情见博客 [封个bash小脚本 | 由文件（夹）的相对路径获取绝对路径](https://zhengshimao.netlify.app/post/get_pwd.html)

使用示例：

#### 获取单个文件绝对路径

```sh
$ rd 1.qc.sh
/home/data/vip13t28/workspace/pig/1.qc.sh

```

#### 获取单个文件夹绝对路径

```sh
$ rd 2.clean_reads/
/home/data/vip13t28/workspace/pig/1.qc.sh
/home/data/vip13t28/workspace/pig/2.clean_reads
```

#### 获取多个文件（夹）

```sh
$ rd 4.rseqc/ 5.stringtie_gtf  nohup.out 
/home/data/vip13t28/workspace/pig/4.rseqc
/home/data/vip13t28/workspace/pig/5.stringtie_gtf
/home/data/vip13t28/workspace/pig/nohup.out

```

#### 获取多个文件绝对路径

```sh
$ ls data/SRR18059* |rd
/home/data/vip13t28/workspace/pig/data/SRR1805929_1.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805929_2.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805930_1.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805930_2.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805931_1.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805931_2.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805932_1.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805932_2.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805933_1.fastq.gz
/home/data/vip13t28/workspace/pig/data/SRR1805933_2.fastq.gz

```

### run_merge_fc_counts_normalization.R

这是一个R脚本。

该脚本用于将转录组中`linux`版本`featurecounts`定量结果合并为`count`矩阵，并根据`featurecounts`提供的"有效基因长度"计算`FPKM`和`TPM`。

博客地址：[封个R脚本 | featurecounts 结果合并、计算FPKM/TPM一步到位](https://zhengshimao.netlify.app/run_merge_fc_counts_normalization.r)

#### 脚本依赖

```sh
conda install -c conda-forge r-argparser #0.7
conda install -c conda-forge r-tidyverse
conda install -c bioconda bioconductor-edger #用了rpkm()函数算的FPKM
```

#### 脚本参数介绍与使用

- `-i, --input_path`  指定含有featurecounts结果的文件夹
- `-p, --pattern` 用R中正则指定文件pattern。建议直接将featurecounts的结果文件后缀写为`.count` ，类似`sample.count` 格式，这个选项用默认的就可以了。

- `-o, --output_path` 输出文件结果文件夹指定
- `-f, --prefix ` 输出结果文件名前缀，默认为"`my`"

命令行：

```sh
./run_merge_fc_counts_normalization.R -i ./lncrna  -o ./
```

然后在当前目录有三个文件生成：`my_genes.counts`，`my_genes.fpkm`，`my_genes.tpm`

分别是基因的raw count矩阵，FPKM矩阵和TPM矩阵文件，所有结果文件均未进行基因筛选。

`特别注意`：过程提示会有“`n  genes were not expressed in all samples!`” ，运行时`n`会是具体数值，这句提示的意思是有`n` 个基因在所有样本中的count值均为0。但是，这些基因并未过滤掉！

### run_get_bioinfo_sample.sh

这是个`bash`脚本。

博客地址：[封个Bash脚本 | 脚本根据Project 编号获取3.5种方式从ENA获取数据的脚本与执行命令](https://zhengshimao.netlify.app/post/run_get_bioinfo_sample.sh.html)

该脚本用于根据测序数据的`Project编号`生成从`ENA`下载测序项目数据的`命令`和`脚本`。

每种方式都包含了`数据下载` 、 `数据校验` 和`重命名`等命令，最终的结果文件都是`fastq.gz`。

关于 `run_rename_sample.sh` 的脚本，由于上传数据的人什么奇葩都有，这一步中的脚本必须在运行前进行检查是否符合我们常用的命名规范以及是否能按照区分不同样品。主要是改个名字，做下游分析方便点，不至于眼花。

适用于`单端测序`与`双端测序`。

由于所有脚本的输入和输出都是`在当前目录` 中，所以，屏幕提示的3种下载方式为`主体命令`，需要你在使用时根据需要进行简单修改。

脚本内所用数据`PRJNA881764` 是一个miRNA的单端测序项目，6个数据文件。`PRJNA275632`是一个猪的普通转录组双端测序项目。

#### 脚本依赖

每个方式单独说依赖，不需要全部安装，根据自己选择的方式进行下载。

第一种方式：`kingfisher` 【我一直在用的下载方法】

第二种方式：`ascp` 

第三种方式：`axel` 、`sra-tools` 、 `parallel-fastq-dump` 和`pigz`

屏幕输出的代码中还涉及到了`ParaFly` ，不装也行，直接`bash `相应脚本，就是慢点而已。

#### 脚本使用

项目测序编号`PRJNA881764`

```sh
bash run_get_bioinfo_sample.sh PRJNA881764 #建议将屏幕输出重定向到一个文件内。屏幕输出即为你运行所产生bash脚本的主体命令
```

不知道怎么用了就运行

```sh
bash run_get_bioinfo_sample.sh
```

会获取如下提示：

```sh
Usage:bash run_get_bioinfo_sample.sh PRJNA881764
please give a parameter of project number,e.g. PRJNA275632
```

正确运行会获取一个`ena_info_sample`文件夹

```sh
tree ena_info_sample/
ena_info_sample/
├── fastq_md5.txt
├── progect_infor.txt
├── run_axel_sra_from_ena.sh
├── run_downloaded_fastq_using_aspera.sh
├── run_gzip_parallel-fastq-dump.sh
├── run_pigz_fasterq_dump.sh
├── run_rename_sample.sh
├── run_wget_sra_from_ena.sh
└── sra_md5.txt
```

还会有屏幕提示三种下载方式，如下：

```sh
###################################################
#The 1st method
conda activate kingfisher
kingfisher get -p PRJNA881764 -m ena-ftp aws-http aws-cp prefetch --download-threads 8 -f fastq.gz
md5sum -c fastq_md5.txt > md5.res
#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩
bash run_rename_sample.sh

###################################################
#The 2nd method

ParaFly -c run_downloaded_fastq_using_aspera.sh -CPU 5 -failed_cmds failed_run_downloaded_fastq_using_aspera.txt
md5sum -c fastq_md5.txt > md5.res
#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩
bash run_rename_sample.sh

###################################################
#The 3rd method
ParaFly -c run_wget_sra_from_ena.sh -CPU 10 -failed_cmds failed_wget_sra_from_ena.txt #wget下载sra
ParaFly -c run_axel_sra_from_ena.sh -CPU 1 -failed_cmds failed_axel_sra_from_ena.txt #axel下载sra #设置多CPU可能报错,原因未知
md5sum -c sra_md5.txt > md5.res #校验
ParaFly -c run_pigz_fasterq_dump.sh -CPU 2 -failed_cmds failed_pigz_fasterq_dump.txt #首选
ParaFly -c run_gzip_parallel-fastq-dump.sh -CPU 5 -failed_cmds failed_gzip_parallel-fastq-dump.txt #次选
#在进行下一步之前先检查重命名后的名字是否规范,有的项目信息有奇葩
bash run_rename_sample.sh

```

### run_dat2fa_swiss_prot.pl

这是一个perl脚本。

Uniprot数据下目录 [Taxonomic divisions](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/) 存放不同物种大类的`dat` 文件，但是没有`fasta` 文件，这里写了个脚本将该目录下的文件转为fasta格式。

脚本依赖的都是基础模块。如果输入文件是`.gz` 格式，在win上必须解压或者使用gitbash才能处理（因为调用了`zcat` 命令）。

博客地址 [Uniprot数据库dat文件转换为fasta文件 - 学生信的大叔](https://zhengshimao.netlify.app/post/uniprot_dat2fasta.html)

#### 版本升级

##### 2022.12.11 v0.3 【已转为非公kai脚本】

- 添加脚本`run_dat2fa_uniprot.pl` ，原`run_dat2fa_swiss_prot.pl` 未动。`run_dat2fa_uniprot.pl`在原脚本基础上做出如下修改
-  添加了`-f` 选项，即使结果文件已存在，仍然执行，使其覆盖。
- 对每条序列`ID`行添加检查`Reviewed/Unreviewed` （分别对应fasta结果文件序列ID行中的`sp/tr`），使其适用`swiss-prot`与`trembl`数据库`dat`文件转换。
- 思考：每条序列都添加检查可能会使速度稍慢，这可能对于`trembl` 库的dat文件转换才比较明显。如何在指定数据库来源的同时减少内部判断增加的运行时间呢？分块编程重新写？太麻烦，暂时不改了。

##### 2022.11.22 v0.2 【仍公开】

- 将脚本名称由原来的`dat2fa_swiss_prot.pl` 改为了`run_dat2fa_swiss_prot.pl` 
- 新增了fasta结果信息文件`.info` 。将fasta序列标题行信息整理成了tab分隔的文件，方便处理。不过得多说一句，`.info` 文件只是提取了部分`.dat` 文件信息。
- 修改了输出结果文件命名处理部分。原来的处理部分重复了；处理`.gz` 文件信息时在命名上也会带上`.gz`让我觉得别扭，也修改了。
- 所有修改行在对应行标记了修改日期`2022.11.22` 以及修改的方式。

#### FASTA序列信息注解

https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

```sh
>sp|P31237|ACCO_ACTDE 1-aminocyclopropane-1-carboxylate oxidase OS=Actinidia deliciosa OX=3627 GN=ACO PE=2 SV=1
MEAFPVIDMEKLNGEERAPTMEKIKDACENWGFFELVNHGISHELMDTVERLTKEHYNKC
MEQRFKEMVATKGLEAVQSEINDLDWESTFFLRHLPVSNISEIPDLEQDHRKAMKEFAEK
LEKLAEQLLDLLCENVGLEKGYLKKAFYGSKGPTFGTKVSNYPPCPRPELIKGLRAHTDA
GGIILLFQDNKVSGLQLLKDGEWIDVPPMKHSIVINIGDQLEVITNGKYKSVMHRVIAQP
DGNRMSIASFYNPGSDAVMYPAPALVDKEEDQQKQVYPKFVFEDYMKLYAGLKFQAKEPR
FEAMKAMENAVNLGPIATI
```

\> 后的注释信息进行简单说明。

- **sp**：Swiss-Prot数据库的简称，也就是上面说的验证后的蛋白数据库。

那Trembl数据库呢？未下载过，没看过。

- **P31237**：UniProt ID号   对应dat文件中的`AC  P31237;`  <u>dat文件中的编号可能有多个，并以`;` 分割，而在fasta序列中只包含第一个ID号。</u>

- **ACCO_ACTDE 1-aminocyclopropane-1-carboxylate oxidase**：蛋白质名称。`ACCO_ACTDE`对应dat文件中的`ID`，`1-aminocyclopropane-1-carboxylate oxidase`是dat文件中的`DE  RecName: Full=1-aminocyclopropane-1-carboxylate oxidase;`。

- **OS=Actinidia deliciosa** ：OS是Organism简称，Actinidia deliciosa是美味猕猴桃的拉丁文名称，说明该蛋白是来自美味猕猴桃。对应dat文件中的`OS ` 。
- **OX=3627**：Organism Taxonomy，也就是物种分类数据库Taxonom y ID。对应dat文件中的`OX  NCBI_TaxID=3627;`
- **GN=ACO** ：Gene name，基因名为ACO 。对应dat文件中的`GN  Name=ACO;` 。部分序列无`GN项`，如`108_SOLLC` (一种植物蛋白)。<u>部分蛋白中可能没有`GN` 项</u>。
- **PE=2** ：Protein Existence，蛋白质可靠性，对应5个数字，数字越小越可靠：对应dat文件中的`PE  2: Evidence at transcript level;` 
  1：Experimental evidence at protein level
  2：Experimental evidence at tranlevel
  3：Protein inferred from homology
  4：Protein predicted
  5：Protein uncertain
- **SV=1**：Sequence Version，序列版本号。对应dat文件中的`DT  01-JUL-1993, sequence version 1.` 

#### 脚本使用方法

```sh
perl dat2fa_swiss_prot.pl -i uniprot_sprot_plants.dat.gz -o output.fasta
```

- -i 输入文件，是必需参数。允许解压后的dat文件，也允许gzip压缩的`.gz` 后缀文件。
- -o 指定输出文件。如果未指定，则根据输入文件进行定义，且输出路径与输入文件相同。如果结果文件已经存在，则无法运行，必须更换或者删除输出文件名称。

注意：脚本有部分内容（序列ID中的'sp|'）写死，只能用于[Taxonomic divisions](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/) 下的uniprot_sprot文件`uniprot_sprot*.gz` ，不能用于 uniprot_trembl 文件【其实有网友测试了能转格式，应该能用于注释，但是被我写死的`'sp|'` 部分确实不妥，你必须要保证自己备注好真实的数据库来源】。
