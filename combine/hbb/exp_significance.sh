year="2022"

echo "GGF SIGNIFICANCE"
combine -M Significance -m 125 --signif output/testModel_${year}/model_combined.root --cminDefaultMinimizerStrategy 0 -t -1 --redefineSignalPOI rggF --setParameters rggF=1

