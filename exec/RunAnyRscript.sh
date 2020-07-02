#!/bin/bash
#this script run any R script and create output folder
# it should be run thoruhgt qsub thus its arguments need to by passed throught wsub call in this way:
# qsub -v "suf=RANDOM,env=ngrip2,fun=fun,script=exec/focusLast.R" -l nodes=1:ppn=16 exec/RunAnyRscript.sh 
echo $PWD 
cd  /home/share/simon/pleistoclimate/
node=$HOSTNAME 
name=${env}${fun}_Sls${sls}${suf}
nodefoldname=${name}/${name}_${node}

if [ ! -d $name ] 
then 
    echo newdir
    mkdir -p $name
fi


Rscript $script 16 1 ${nodefoldname} ${env} ${fun} ${sls} > ${nodefoldname}.log 2> ${nodefoldname}.err 
