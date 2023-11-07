year="2022"

#combine -M FitDiagnostics -m 125 output/testModel_${year}/model_combined.root --saveShapes --saveWithUncertainties --robustFit=1 --robustHesse=1 --cminDefaultMinimizerStrategy=0 --setParameters rZbb=1 --verbose 9

combine -M FitDiagnostics -d output/testModel_2022/model_combined.root --robustFit=1 --setRobustFitAlgo=Minuit2,Migrad --setRobustFitStrategy=1 --rMin -15 --plots --saveShapes --saveWithUncertainties --setParameters rZbb=1
