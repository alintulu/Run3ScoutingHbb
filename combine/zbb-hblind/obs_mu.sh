year="2022"
frozen="rgx{CMS_.*},rgx{QCDscale_.*},rgx{UEPS_.*},rgx{pdf_.*},rgx{.*mcstat},rgx{qcd.*}"
npoints=1000

#combine -M MultiDimFit -m 125 output/testModel_${year}/model_combined.root --setParameters rggF=1,rZbb=1 --cminDefaultMinimizerStrategy 0 --algo grid --points ${npoints} --redefineSignalPOI rZbb --saveWorkspace -n rZbb --freezeParameters rggF

combine -M MultiDimFit -m 125 --setParameters rggF=1,rZbb=1 --cminDefaultMinimizerStrategy 0 --algo grid --points ${npoints} --redefineSignalPOI rZbb --saveWorkspace -n rZbbStatOnly -d higgsCombinerZbb.MultiDimFit.mH125.root -w w --snapshotName "MultiDimFit" --freezeParameters rggF,${frozen}
