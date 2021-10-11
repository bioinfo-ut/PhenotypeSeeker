#!/bin/bash

sudo apt-get update --assume-yes
sudo apt-get install python3 python3-pip virtualenv

virtualenv -p python3 .PSenv
source .PSenv/bin/activate

cp bin/* ./.PSenv/bin
yes | python3 -m pip install -r requirements.txt
