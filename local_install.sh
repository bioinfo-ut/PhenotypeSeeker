#!/bin/bash

pip3 install --user --upgrade pip
pip3 install --user --upgrade setuptools
pip3 install --user --upgrade --ignore-installed numpy

cp bin/* ~/.local/bin/
pip3 install --user .

sed -n "s:\"mash\":\"~/.local/bin/mash\":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistmaker \":\"~/.local/bin/glistmaker \":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistquery\":\"~/.local/bin/glistquery\":g" PhenotypeSeeker/modeling.py
sed -n "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py
