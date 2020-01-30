#!/bin/bash

sudo apt-get update
sudo apt-get install python3 python3-pip

sudo python3 -m pip install --upgrade pip

sudo cp bin/* /usr/local/bin
sudo python3.6 -m pip install .

