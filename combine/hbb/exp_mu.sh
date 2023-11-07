year="2022"

npoints=100

combine -M MultiDimFit -m 125 output/testModel_${year}/model_combined.root --setParameters rggF=1 -t -1 --cminDefaultMinimizerStrategy 0 --algo grid --points ${npoints} --redefineSignalPOI rggF --saveWorkspace -n rggF
