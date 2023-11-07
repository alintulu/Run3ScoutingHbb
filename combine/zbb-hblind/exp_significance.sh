year="2022"

combine -M Significance -m 125 --signif output/testModel_${year}/model_combined.root --cminDefaultMinimizerStrategy 0 -t -1 --redefineSignalPOI rZbb --setParameters rZbb=1

