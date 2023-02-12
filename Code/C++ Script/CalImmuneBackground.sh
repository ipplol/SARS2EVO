#该脚本用于计算ImmuneBackground
#!/bin/bash

thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir1=`dirname $thisfold`
mainDir=`dirname $mainDir1`
echo $mainDir

k=0
cat $mainDir/TaskList.tsv | awk 'NR>1' | while read line
do
let "k++"
#echo $line
input1=${line// /_}
input2=`echo ${input1//_/+}`
input3=`echo ${input2// /_}`
input4=`echo ${input3//-/}`
echo $input1
array=(${input4//_/ })
echo ${array[0]}
echo ${array[2]}
k=$((k%16))

#$mainDir/AntibodySpectrum/Script/CalInfectBackground.out $mainDir/SelectiveCoEfficient/Data $mainDir/AntibodySpectrum/Result $input1 &
$mainDir/AntibodySpectrum/Script/CalAntibodyPressure.out $mainDir/AntibodySpectrum/Result/${array[0]}_20191201_${array[2]}.InfectionList $mainDir/AntibodySpectrum/Data &
#16 threads
if [ $k == 0 ]; then
wait
fi

done

wait
