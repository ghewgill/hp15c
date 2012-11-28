#!/bin/sh -e

go() {
    convert -geometry $1x$1 iTunesArtwork.png $2
}

convert -background none -geometry 300x450 1.svg 1.tmp.png
convert -background none -geometry 300x450 5.svg 5.tmp.png
convert -background none -geometry 300x450 C.svg C.tmp.png

convert 1.tmp.png \( +clone -background black -shadow 60x3+20+30 \) +swap -background none -layers merge +repage 1.tmp.png
convert 5.tmp.png \( +clone -background black -shadow 60x3+20+30 \) +swap -background none -layers merge +repage 5.tmp.png
convert C.tmp.png \( +clone -background black -shadow 60x3+20+30 \) +swap -background none -layers merge +repage C.tmp.png

convert icon-base.jpg \
    1.tmp.png -geometry -100+300 -composite \
    5.tmp.png -geometry +300+300 -composite \
    C.tmp.png -geometry +700+300 -composite \
    iTunesArtwork.png
rm *.tmp.png
#go 512
go 144 Icon-72@2x.png
go 114 Icon@2x.png
go 72 Icon-72.png
go 57 Icon.png
