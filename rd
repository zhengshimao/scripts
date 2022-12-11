#!/usr/bin/env bash

function get_pwd(){
    here=$PWD;
    if test -d $1 ;then
    cd $1;
    dir=$PWD;
    echo $dir;
    cd $here ;
elif test -f $1; then
    #here=$PWD;
    tmp_file=`basename $1`;
    tmp_dir=`dirname $1`;
    cd $tmp_dir;
    dir=$PWD;
    echo "${dir}/${tmp_file}"
    cd $here ;
else 
    echo "$1 : No such file or dirctory!"
fi
}
#get_pwd $1
#<<EOF

if [[ $# != 0 ]];then
    for i in $@ ;do
        get_pwd $i
    done
else
    while read line;
    do
        get_pwd $line
    done
fi
#EOF
