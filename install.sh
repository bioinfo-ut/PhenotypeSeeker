#!/bin/bash

sudo apt-get install software-properties-common
sudo apt-add-repository universe
sudo apt-get update
sudo apt-get install python2.7 python-pip python-dev

sudo pip2 install --upgrade pip
sudo pip2.7 install --upgrade setuptools
sudo pip2.7 install --upgrade --ignore-installed numpy

sudo cp bin/* /usr/bin
sudo pip2.7 install .

