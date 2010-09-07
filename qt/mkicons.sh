#!/bin/sh

convert 15.png -geometry 128 -bordercolor none -border 0x100 -crop 128x128+0+77 15-128.png
convert 15.png -geometry 48 -bordercolor none -border 0x100 -crop 48x48+0+92 15-48.png
convert 15.png -geometry 32 -bordercolor none -border 0x100 -crop 32x32+0+95 15-32.png
convert 15.png -geometry 16 -bordercolor none -border 0x100 -crop 16x16+0+97 15-16.png
