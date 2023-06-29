#!/bin/sh

python ExtractEdges.py --emFile data_230410/All_Original_data.txt --annFile data_230410/All_Original_meta.txt --coreNum 14 --out data_230410/All_Original_data
python ExtractEdges.py --emFile data_230410/All_PhaseX_data.txt --annFile data_230410/All_PhaseX_meta.txt --coreNum 14 --out data_230410/All_PhaseX_data
python DiffEdges.py --refFolder data_230410/All_Original_data --targetFolder data_230410/All_PhaseX_data --out data_230410/All_PhaseX-Original_mean
