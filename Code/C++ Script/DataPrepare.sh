#该脚本用于编译c++代码，准备数据
#!/bin/bash

echo "Compiling [Some warnings may show up which are normal, please relax.]:"
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#编译脚本
cd `dirname "thisfold"`
g++ $thisfold/HaploCaseCalculator/main.cpp -o HaploCaseCalculator.out -std=c++11 -lpthread
g++ $thisfold/Mut2Haplo/main.cpp -o Mut2Haplo.out -std=c++11 -lpthread
g++ $thisfold/SelectCOEFCalculator/main.cpp -o SelectCOEFCalculator.out -std=c++11 -lpthread
g++ $thisfold/SeqFindMeta/main.cpp -o SeqFindMeta.out -std=c++11 -lpthread
g++ $thisfold/CasePerMillion/main.cpp -o CasePerMillion.out -std=c++11 -lpthread

#更新数据
mainDir1=`dirname $thisfold`
mainDir=`dirname $mainDir1`
cd "$mainDir/SelectiveCoEfficient/Data"
#wget下载

#预处理
echo "Data processing. It may take a while."
#$mainDir/SelectiveCoEfficient/Script/SeqFindMeta.out $mainDir/SelectiveCoEfficient/Data
#$mainDir/SelectiveCoEfficient/Script/Mut2Haplo.out $mainDir/SelectiveCoEfficient/Data/sample_mut_loc_time.tsv
#$mainDir/SelectiveCoEfficient/Script/CasePerMillion.out $mainDir/SelectiveCoEfficient/Data

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`