year="2022"

combine -M FitDiagnostics -m 125 output/testModel_${year}/model_combined.root --setParameters rZbb=1 -t -1 --saveShapes --saveWithUncertainties --cminDefaultMinimizerStrategy 0 --robustFit=1
