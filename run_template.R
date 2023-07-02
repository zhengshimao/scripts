#!/usr/bin/env Rscript

Sys.setenv(LANGUAGE = "en") # 设置提示信息为英文
suppressPackageStartupMessages(library(magrittr, quietly = T))


# # 所需安装包，仅作展示查看，不做运行。
# {
#   library(argparse) # 该模块依赖于Python的同名模块。
#   library(prettyunits)
#   library(lubridate)
#   library(stringr)
#   library(vroom)
#   library(dplyr) 
# }

# conda安装所需模块
##mamba install -c conda-forge r-argparse  r-prettyunits r-lubridate # r-stringr r-vroom r-dplyr

# 设置脚本选项参数
{
  parser <- argparse::ArgumentParser(description = "Get bed file of ssr sequence！")
  # 添加选项
  parser$add_argument("-i","--input",action = "store",type = "character", 
                      help = "The input.")
  parser$add_argument("-t","--thread",action = "store",type = "integer", default = 10,
                      help = "Threads. [default: 10]")
  
  # type: 'logical', 'integer', 'double' or 'character'
  
  # 解析参数为一个list
  args <- parser$parse_args()
  # args #打印查看参数内容
  parser$print_help() #打印帮助文档,用于调试 # 测试完毕需要注释掉这一行。
}


# 时间设置：起始时间记录。
{
  start_time <- lubridate::now(tzone = 'Asia/Shanghai') #记录开始时间
  cli::cat_line("Start time: ",cli::col_red(start_time),"!")
}


# 将选项转换为常规变量
{
  input <- args$input #输入文件
}

# 调试用参数：封脚本时需要改为F。
if(F){
  input <- "./Trinity.fa.misa"
}

# 检查文件
{
  # 输入文件：两个，其中misa文件为显，fasta文件也必须包含。
  ## 直接检查misa文件
  if(file.exists(misa)){
    cli::cli_inform(message = paste0("OK: ",cli::col_red(input)))
  }else{
    cli::cli_abort(message = c("{cli::col_red(input)} does not exist!",
                               "x"="You've supplied {cli::col_red(input)} in current directory!"))
  }
  
  ## 检查fasta文件
  fasta_file <- stringr::str_remove(input,pattern = "\\.misa")
}

# 定义函数
{
  # WGCNA::collectGarbage # 不讲武德地直接复制过来代码了,因为WGCNA安装常常不是很顺。
  wgcna_collectGarbage <- function(){
    while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1, 4]) {
    }
  }
  
  # 读入fasta为df
  read_bio_fasta <- function(fasta){
    fa <- vroom::vroom_lines(file = fasta, skip_empty_rows = TRUE)
    # a vector of fasta id
    id <- fa[stringr::str_detect(fa, "^>")] %>% stringr::str_remove(pattern = "^>")
    # a vector of fasta sequence
    seq <- stringr::str_replace_all(fa,pattern = "^>.*",replacement = ">") %>% paste0(collapse = "") %>% stringr::str_split(pattern = ">")
    seq <- seq[[1]][-1]
    
    return(data.frame(id = id, sequence = seq))
  }
}

# 脚本代码⭐⭐⭐
{
  # 写入你的内容
  
}


# 回收内存：根据需要使用
{
  # 可以使用简单的invisible(gc())
  # invisible(gc())
  # 如果安装了R包WGCNA，还可以调用WGCNA::collectGarbage 
  wgcna_collectGarbage()
}

# 计算运行用时
{
  # 时间设置：运行时间记录
  end_time <- lubridate::now(tzone = 'Asia/Shanghai') #记录结束时间
  cli::cat_line("End time: ",cli::col_red(end_time),"!") # 打印结束时间
  
  # 计算时间差
  runtime <- lubridate::interval(start = start_time, end = end_time) %>% 
    lubridate::time_length(unit = "second")
  
  cli::cat_line("Time consuming: ", cli::col_red(prettyunits::pretty_sec(runtime)), "!")# 使用prettyunits美化时间差。
}

