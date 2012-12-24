nmake -f Makefile.release clean
nmake -f Makefile.release
"c:\program files\inno setup 5\iscc.exe" hp15c.iss

set QT=c:\qt\4.8.4\bin
erase HP15C.zip
zip -j -X -9 HP15C.zip release\HP15C.exe msvcr100.dll msvcp100.dll %QT%\QtCore4.dll %QT%\QtGui4.dll %QT%\QtScript4.dll
