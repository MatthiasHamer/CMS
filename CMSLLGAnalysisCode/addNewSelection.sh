#!/bin/bash

if [ ! -e cutflow${1}.cpp ]; then 
echo "#include \"LLGAnalysis.h\"" >> cutflow${1}.cpp
echo "void LLGAnalysis::Setup${1}() {" >> cutflow${1}.cpp
echo  >> cutflow${1}.cpp
echo "    // setup the cutflow">> cutflow${1}.cpp
echo  >> cutflow${1}.cpp
echo "    // and the histograms" >> cutflow${1}.cpp
echo  >> cutflow${1}.cpp
echo "    return;">> cutflow${1}.cpp
echo "}">> cutflow${1}.cpp
echo  >> cutflow${1}.cpp
echo "void LLGAnalysis::${1}Selection() {">> cutflow${1}.cpp
echo  >> cutflow${1}.cpp
echo "    return;">> cutflow${1}.cpp
echo "}">> cutflow${1}.cpp
fi

sed -ie '/INSERT YOUR SELECTION HERE/ i\        void Setup'${1}'();' LLGAnalysis.h
sed -ie '/INSERT YOUR SELECTION HERE/ i\        void '${1}'Selection();' LLGAnalysis.h
sed -ie '/SETUP YOUR SELECTION HERE/ i\    else if( SELECTION == "'${1}'" ) Setup'${1}'();' LLGAnalysis.cpp
sed -ie '/CALL YOUR SELECTION HERE/ i\        else if( SELECTION == "'${1}'" ) '${1}'Selection();' LLGAnalysis.cpp
sed -ie 's/cutflowSignalRegion.o/cutflowSignalRegion.o cutflowWJetsCR.o/' Makefile
