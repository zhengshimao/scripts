#!/usr/bin/env Rscript

Sys.setenv(LANGUAGE = "en") # 设置提示信息为英文
options(stringsAsFactors=F)
args<-commandArgs(T) 
# class(args) # [1] "character"  


{
  USAGE <- "
  =================================================
  Descripthon: the function of this script!
  Usage:       run_template_base.r <input> <output>
                <input> input file
                <output> output file
  Author:      Shimao Zheng, zhengshimao007@163.com
  Version:     v0.1
  Create:      2023.07.23
  Update:      -
  NOTE:        Follow my WeChat public account \"The_Elder_Student!<学生信的大叔>\"
  =================================================
  "
}

{
  input <- args[1] #<- "test.txt" #输入文件
  output <- args[2] #<- "out.txt"
  if(length(args) < 2) stop(USAGE)
}

# 调试用参数：封脚本时需要改为F。
if(F){
  input <- "./test.txt"
  output <- "./out.txt"
}

# 检查文件
{
  ## 检查输入文件
  if(file.exists(input)){
    message("OK: ", input)
  }else{
    stop("File Not Exists: ",input,"!\n",USAGE)
  }
  
  
}

# 时间设置：起始时间记录。
{
  start_time <- Sys.time() #记录开始时间
  message("Start time: ",start_time,"!")
}


# 定义函数
{
  # WGCNA::collectGarbage # 不讲武德地直接复制过来代码了,因为WGCNA安装常常不是很顺。
  wgcna_collectGarbage <- function(){
    while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1, 4]) {
    }
  }
  
  # 本脚本函数
}


# 脚本代码⭐⭐⭐
{
  # 写入你的需求
  print(input)
  print(output)
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
  end_time <- Sys.time() #记录结束时间
  message("End time: ",end_time,"!") # 打印结束时间
  
  # 计算时间差
  runtime <- difftime(time1 = end_time, time2 = start_time, units = "auto")
  message("==>Time consuming: ",round(runtime, digits = 2)," ",units(runtime),"!") # units(runtime)获取单位
}












