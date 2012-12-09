#!/bin/sh

convert -size 480x320 canvas:black -page +0+20 ../../../images/calc.jpg -rotate 270 -layers flatten Default.png
convert -size 960x640 canvas:black -page +0+40 ../../../images/calc@2x.jpg -rotate 270 -layers flatten Default@2x.png
convert -size 1136x640 canvas:black -page +88+40 ../../../images/calc@2x.jpg -rotate 270 -layers flatten Default-568h@2x.png

convert ../../../images/calc@2x.jpg \
    -page +174+40 ../../../images/ios/1@2x.png \
    -page +208+80 ../../../images/ios/decimal@2x.png \
    -page +220+40 ../../../images/ios/6@2x.png \
    -page +266+40 ../../../images/ios/1@2x.png \
    -page +312+40 ../../../images/ios/8@2x.png \
    -page +358+40 ../../../images/ios/0@2x.png \
    -page +404+40 ../../../images/ios/3@2x.png \
    -page +450+40 ../../../images/ios/3@2x.png \
    -page +496+40 ../../../images/ios/9@2x.png \
    -page +542+40 ../../../images/ios/8@2x.png \
    -page +588+40 ../../../images/ios/9@2x.png \
    -layers flatten screenshot-front.jpg
convert -size 1136x600 canvas:black -page +88+0 screenshot-front.jpg -layers flatten screenshot-front-568h.jpg

convert -size 960x600 canvas:black -page +0+55 ../../../images/back@2x.jpg -geometry 960 -layers flatten screenshot-back.jpg
convert -size 1136x600 canvas:black ../../../images/back@2x.jpg -geometry 1136 -layers flatten screenshot-back-568h.jpg
