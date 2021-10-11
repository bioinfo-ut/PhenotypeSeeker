#!/bin/bash

sudo apt-get update
sudo apt-get install python3 python3-pip virtualenv

virtualenv -p python3 .PSenv
source .PSenv/bin/activate

cp bin/* ./PSenv/bin
python3 -m pip install -r requirements.txt
python3 -m pip install .

deactivate
