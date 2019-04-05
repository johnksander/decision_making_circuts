#!/usr/bin/env bash

DIR=$1
#DIR=diag_EtoIfixed
#DIR=diag_ItoEfixed
OUT=figdirs_$DIR
mkdir $OUT

for d in $DIR/*/
do
    d=${d%*/}
    #echo ${d##*/}
    targ_dir=$DIR/${d##*/}
    echo "copying $targ_dir"
    cp -r $targ_dir $OUT
done
