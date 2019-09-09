#!/bin/bash

python3 -m pip install --user --upgrade pip

cp bin/* ~/.local/bin/

sed -ni "s:\"mash :\"~/.local/bin/mash :g" PhenotypeSeeker/modeling.py
sed -ni "s:\"mash\":\"~/.local/bin/mash\":g" PhenotypeSeeker/modeling.py
sed -ni "s:\"glistmaker:\"~/.local/bin/glistmaker:g" PhenotypeSeeker/modeling.py
sed -ni "s:\"glistquery\":\"~/.local/bin/glistquery\":g" PhenotypeSeeker/modeling.py
sed -ni "s:gmer_counter:~/.local/bin/gmer_counter:g" PhenotypeSeeker/prediction.py

python3 -m pip install --user .