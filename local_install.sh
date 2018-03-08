#!/bin/bash

pip install --user --upgrade pip
sudo pip install --user --upgrade setuptools
sudo pip install --user --upgrade numpy

sudo cp bin/* /usr/bin
sudo pip install .