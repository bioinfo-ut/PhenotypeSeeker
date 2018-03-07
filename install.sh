#!/bin/bash

sudo apt-get update
sudo apt-get install python2.7 python-pip python-dev

pip install --upgrade pip
sudo pip install --upgrade setuptools
sudo pip install --upgrade numpy

sudo cp bin/* /usr/bin
sudo pip install .

echo
phenotypeseeker --version
