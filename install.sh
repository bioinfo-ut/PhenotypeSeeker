#!/bin/bash

if [ "$1" != "--user" ]; then
  sudo apt-get update --assume-yes
  sudo apt-get install --assume-yes python3 python3-pip virtualenv
fi

virtualenv -p python3 .PSenv --system-site-packages
source .PSenv/bin/activate

cp bin/* ./.PSenv/bin
python3 -m pip install setuptools==58.2.0
yes | python3 -m pip install -r pkgs_to_inst.txt
