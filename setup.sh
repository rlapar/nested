#!/bin/sh

#TODO ADD CONFIG PATH

python3 setup.py install

#remove unnecessary dirs
rm -rf build dist nested.egg-info

#move config file
mkdir -p /etc/nested
cp config.yml /etc/nested/.
cp nested/config/gt.style /etc/nested/.