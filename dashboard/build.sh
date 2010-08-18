#!/bin/sh -e

DEST=HP15C.wdgt

rm -rf $DEST
mkdir $DEST
cp ../15.jpg $DEST/
cp ../*.png $DEST/
cp ../sprintf-0.6.js ../hp15c.js $DEST/
mkdir $DEST/jsmat
cp ../jsmat/matrix.js $DEST/jsmat/
cp Default.png $DEST/
cp Info.plist $DEST/
cp hp15c.html $DEST/
rm $DEST.zip
zip -r $DEST.zip $DEST
