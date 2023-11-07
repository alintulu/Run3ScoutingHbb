year="2022"

cd output/testModel_${year}/

. build.sh

text2workspace.py model_combined.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=.*/ZJetsbb:rZbb[1,0,3]' --PO verbose

cd ../..
