#!/bin/bash

sudo apt-get update
sudo apt-get install python2.7 python-pip python-dev

sudo python2.7 -m pip install --upgrade pip
sudo python2.7 -m pip install --upgrade setuptools
sudo python2.7 -m pip install --upgrade numpy

sudo cp bin/* /usr/bin
sudo python2.7 -m pip install .

