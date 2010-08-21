#!/bin/sh -e

set -e

make clean
rm -rf hp15c.app
make
#rm hp15c.app.zip
#zip -r hp15c.app.zip hp15c.app
macdeployqt hp15c.app -dmg
