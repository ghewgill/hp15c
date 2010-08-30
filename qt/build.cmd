pause "build release version"
"c:\program files\inno setup 5\iscc.exe" hp15c.iss

set QT=c:\qt\2010.02.1\qt\bin
erase HP15C.zip
zip -j -X -9 HP15C.zip release\HP15C.exe %QT%\libgcc_s_dw2-1.dll %QT%\mingwm10.dll %QT%\QtCore4.dll %QT%\QtGui4.dll %QT%\QtScript4.dll
