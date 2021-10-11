#!/bin/bash

sudo apt-get update
sudo apt-get install python3 python3-pip virtualenv

#python3-setuptools
#sudo python3 -m pip install --upgrade pip

virtualenv -p python3 .PSenv
cp bin/* ./PSenv/bin
python3 -m pip install -r requirements.txt
python3 -m pip install .

