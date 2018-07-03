echo "start Iter Temp Fit"
rm -r IterationProgress
rm -r MCTemplatesAnData
rm -r MixedBGComp
mkdir IterationProgress
mkdir MCTemplatesAnData
mkdir MixedBGComp

root -l -b -q Iter_Temp_Fit.C++

echo "start BGComparison"

root -l -b -q MixedBGComp.C++

#going to GammaCalo-All_503_normal_and_extra and making everything there
cd ..
cd GammaCalo-All_503_normal_and_extra

cp -u ../GammaCalo-MC_503/Iter_Temp_Fit.C IterTemp.C
cp -u ../GammaCalo-MC_503/MixedBGComp.C MixedBGComp.C
cp -u ../GammaCalo-MC_503/CommonHeader.h CommonHeader.h


echo "start Iter Temp Fit"
rm -r IterationProgress
rm -r MCTemplatesAnData
rm -r MixedBGComp
mkdir IterationProgress
mkdir MCTemplatesAnData
mkdir MixedBGComp

root -l -b -q Iter_Temp_Fit.C++

echo "start BGComparison"

root -l -b -q MixedBGComp.C++


#going to GammaCalo-All_503 and making everything there
cd ..
cd GammaCalo-All_503

cp -u ../GammaCalo-MC_503/Iter_Temp_Fit.C IterTemp.C
cp -u ../GammaCalo-MC_503/MixedBGComp.C MixedBGComp.C
cp -u ../GammaCalo-MC_503/CommonHeader.h CommonHeader.h

echo "start Iter Temp Fit"
rm -r IterationProgress
rm -r MCTemplatesAnData
rm -r MixedBGComp
mkdir IterationProgress
mkdir MCTemplatesAnData
mkdir MixedBGComp

root -l -b -q Iter_Temp_Fit.C++

echo "start BGComparison"

root -l -b -q MixedBGComp.C++

#going to GammaCalo_503_extra and making everything there
cd ..
cd GammaCalo_503_extra

cp -u ../GammaCalo-MC_503/Iter_Temp_Fit.C IterTemp.C
cp -u ../GammaCalo-MC_503/MixedBGComp.C MixedBGComp.C
cp -u ../GammaCalo-MC_503/CommonHeader.h CommonHeader.h

echo "start Iter Temp Fit"
rm -r IterationProgress
rm -r MCTemplatesAnData
rm -r MixedBGComp
mkdir IterationProgress
mkdir MCTemplatesAnData
mkdir MixedBGComp

root -l -b -q Iter_Temp_Fit.C++

echo "start BGComparison"

root -l -b -q MixedBGComp.C++
echo ""
echo "DONE!"
echo ""

cd ..
cd GammaCalo-MC_503
