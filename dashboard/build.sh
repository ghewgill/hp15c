#!/bin/sh -e

set -e

DEST=HP15C.wdgt

rm -rf $DEST
mkdir $DEST
cp ../images/15.jpg $DEST/
cp ../images/*.png $DEST/
cp ../common/sprintf-0.6.js ../common/hp15c.js $DEST/
mkdir $DEST/jsmat
cp ../common/jsmat/matrix.js $DEST/jsmat/
cp Default.png $DEST/
cp Info.plist $DEST/
cp hp15c.html $DEST/
rm $DEST.zip
zip -r $DEST.zip $DEST
