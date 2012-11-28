#!/bin/sh

convert -size 480x320 canvas:black -page +0+20 ../../../images/calc.jpg -rotate 270 -layers flatten Default.png
convert -size 960x640 canvas:black -page +0+40 ../../../images/calc@2x.jpg -rotate 270 -layers flatten Default@2x.png
convert -size 1136x640 canvas:black -page +88+40 ../../../images/calc@2x.jpg -rotate 270 -layers flatten Default-568h@2x.png
