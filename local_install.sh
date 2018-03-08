#!/bin/bash

pip install --user --upgrade pip
sudo pip install --user --upgrade setuptools
sudo pip install --user --upgrade numpy

cp bin/* ~/.local/bin/
pip install --user .