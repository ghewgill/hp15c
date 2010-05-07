import os

Digits = [
    [' _ ', '   ', ' _ ', ' _ ', '   ', ' _ ', ' _ ', ' _ ', ' _ ', ' _ '],
    ['| |', '  |', ' _|', ' _|', '|_|', '|_ ', '|_ ', '  |', '|_|', '|_|'],
    ['|_|', '  |', '|_ ', ' _|', '  |', ' _|', '|_|', '  |', '|_|', ' _|'],
]

Segments = [
    "M 3 1 L 17 1 14 4 6 4 Z",
    "M 2 3 L 5 6 5 11 2 14 Z",
    "M 18 3 L 18 14 15 11 15 6 Z",
    "M 6 13 L 14 13 16 14.5 14 16 6 16 4 14.5 Z",
    "M 2 15 L 5 18 5 23 2 26 Z",
    "M 18 15 L 18 26 15 23 15 18 Z",
    "M 3 28 L 17 28 14 25 6 25 Z",
]

for i in range(len(Digits[0])):
    name = "0123456789"[i]
    f = open(name + ".svg", "w")
    print >>f, """<?xml version="1.0"?>
        <svg width="20" height="30">
            <g fill="black">"""
    if Digits[0][i][1] == '_': print >>f, """<path d="%s" />""" % Segments[0]
    if Digits[1][i][0] == '|': print >>f, """<path d="%s" />""" % Segments[1]
    if Digits[1][i][2] == '|': print >>f, """<path d="%s" />""" % Segments[2]
    if Digits[1][i][1] == '_': print >>f, """<path d="%s" />""" % Segments[3]
    if Digits[2][i][0] == '|': print >>f, """<path d="%s" />""" % Segments[4]
    if Digits[2][i][2] == '|': print >>f, """<path d="%s" />""" % Segments[5]
    if Digits[2][i][1] == '_': print >>f, """<path d="%s" />""" % Segments[6]
    print >>f, """ </g>
        </svg>"""
    f.close()
    os.system("convert -background none %s.svg %s.png" % (name, name))

f = open("decimal.svg", "w")
print >>f, """<?xml version="1.0"?>
    <svg width="3" height="3">
        <rect width="100%" height="100%" fill="black" />
    </svg>"""
f.close()
os.system("convert -background none decimal.svg decimal.png")

f = open("neg.svg", "w")
print >>f, """<?xml version="1.0"?>
    <svg width="12" height="3">
        <rect width="100%" height="100%" fill="black" />
    </svg>"""
f.close()
os.system("convert -background none neg.svg neg.png")
