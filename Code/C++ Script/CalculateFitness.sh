#该脚本用于计算fitness
#!/bin/bash

thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir1=`dirname $thisfold`
mainDir=`dirname $mainDir1`

cat $mainDir/TaskList.tsv | awk 'NR>1' | while read line
do
echo $line
input1=${line// /_}
input2=`echo ${input1//_/+}`
input3=`echo ${input2// /_}`
input4=`echo ${input3//-/}`
echo $input1
echo $input4
$mainDir/SelectiveCoEfficient/Script/HaploCaseCalculator.out $mainDir/SelectiveCoEfficient/Data $mainDir/SelectiveCoEfficient/Result $input1
Rscript $mainDir/SelectiveCoEfficient/Script/EpiEstim.R $mainDir/SelectiveCoEfficient/Result/$input4.DailyCases
$mainDir/SelectiveCoEfficient/Script/SelectCOEFCalculator.out $mainDir/SelectiveCoEfficient/Result/$input4.LaterProp
done
