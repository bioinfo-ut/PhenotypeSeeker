#!/bin/bash

python2.7 -m pip install --user --upgrade pip
python2.7 -m pip install --user --upgrade setuptools
python2.7 -m pip install --user --upgrade numpy

cp bin/* ~/.local/bin/
python2.7 -m pip install --user .

sed -n "s:\"mash\":\"~/.local/bin/mash\":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistmaker \":\"~/.local/bin/glistmaker \":g" PhenotypeSeeker/modeling.py
sed -n "s:\"glistquery\":\"~/.local/bin/glistquery\":g" PhenotypeSeeker/modeling.py
sed -n "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py
