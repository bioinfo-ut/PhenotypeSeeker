#!/bin/bash

python3 -m pip install --user --upgrade pip

cp bin/* ~/.local/bin/

sed -i "s:\"mash :\"~/.local/bin/mash :g" PhenotypeSeeker/modeling.py
sed -i "s:\"glistmaker:\"~/.local/bin/glistmaker:g" PhenotypeSeeker/modeling.py
sed -i "s:\"glistquery:\"~/.local/bin/glistquery:g" PhenotypeSeeker/modeling.py
sed -i "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py

sed -i "s:^phenotypeseeker:~/.local/bin/phenotypeseeker:" example/test_PS_modeling.sh
sed -i "s:^phenotypeseeker:~/.local/bin/phenotypeseeker:" example/test_PS_prediction.sh

python3 -m pip install --user .