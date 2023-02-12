#该脚本用于计算ImmuneBackground
#!/bin/bash

thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir1=`dirname $thisfold`
mainDir=`dirname $mainDir1`
echo $mainDir

#compile
cd $mainDir/AntibodySpectrum/Script
g++ $mainDir/AntibodySpectrum/Script/CalEscapeScoreMatrix/main.cpp -o CalEscapeScoreMatrix.out -std=c++11 -lpthread
g++ $mainDir/AntibodySpectrum/Script/CalAntibodyPressure/main.cpp -o CalAntibodyPressure.out -std=c++11 -lpthread
g++ $mainDir/AntibodySpectrum/Script/CalInfectBackground/main.cpp -o CalInfectBackground.out -std=c++11 -lpthread
g++ $mainDir/AntibodySpectrum/Script/CorrelationBuild_PressureFitness/main.cpp -o CorrelationBuild_PressureFitness.out -std=c++11 -lpthread

cd $mainDir

#transfer escape score list to matrix
mkdir $mainDir/AntibodySpectrum/Data/EscapeScore
#$mainDir/Script/CalEscapeScoreMatrix.out $mainDir/Data/use_abs_res.csv $mainDir/Data/EscapeScore

#plot heatmap
#for f in `find $mainDir/AntibodySpectrum/Data/EscapeScore "*.EscapeMatrix"`
#do 
#Rscript $mainDir/AntibodySpectrum/Script/ESSMplot.R $f
#done
