#!/bin/bash

pip2.7 install --user --upgrade pip
pip2.7 install --user --upgrade setuptools
pip2.7 install --user --upgrade numpy

cp bin/* ~/.local/bin/
pip2.7 install --user .

sed -n "s:\"mash\":\"~/.local/bin/mash\":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistmaker \":\"~/.local/bin/glistmaker \":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistquery\":\"~/.local/bin/glistquery\":g" PhenotypeSeeker/modeling.py
sed -n "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py
