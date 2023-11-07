year="2022"

combine -M FitDiagnostics -m 125 output/testModel_${year}/model_combined.root --robustFit=1 --saveShapes --saveWithUncertainties --saveOverallShapes --setParameters rggF=1 --cminDefaultMinimizerStrategy 0 --ignoreCovWarning
