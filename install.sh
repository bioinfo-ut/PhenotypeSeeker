#!/bin/bash

sudo apt-get update --assume-yes
sudo apt-get install --assume-yes python3 python3-pip python3-virtualenv

virtualenv -p python3 .PSenv --system-site-packages
source .PSenv/bin/activate

cp bin/* ./.PSenv/bin
yes | python3 -m pip install -r pkgs_to_inst.txt
