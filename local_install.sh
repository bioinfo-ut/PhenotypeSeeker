#!/bin/bash

pip install --user --upgrade pip
pip install --user --upgrade setuptools
pip install --user --upgrade numpy

cp bin/* ~/.local/bin/
pip install --user .

sed -n "s:\"mash\":\"~/.local/bin/mash\":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistmaker \":\"~/.local/bin/glistmaker \":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistquery\":\"~/.local/bin/glistquery\":g" PhenotypeSeeker/modeling.py
sed -n "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py
