[Setup]
AppId={{9280DCF0-6A26-4CC1-96E8-ECB1E355B32B}
AppName=HP15C
AppVersion=1.0
;AppVerName=HP15C 1.0
AppPublisher=hewgill.com
AppPublisherURL=http://hewgill.com
AppSupportURL=http://hewgill.com
AppUpdatesURL=http://hewgill.com
DefaultDirName={pf}\HP15C
DefaultGroupName=HP15C
DisableProgramGroupPage=yes
OutputBaseFilename=HP15C-win32-install
OutputDir=.
Compression=lzma
SolidCompression=yes
PrivilegesRequired=lowest

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "release\hp15c.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "c:\qt\2010.02.1\qt\bin\libgcc_s_dw2-1.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "c:\qt\2010.02.1\qt\bin\mingwm10.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "c:\qt\2010.02.1\qt\bin\QtCore4.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "c:\qt\2010.02.1\qt\bin\QtGui4.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "c:\qt\2010.02.1\qt\bin\QtScript4.dll"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\HP15C"; Filename: "{app}\hp15c.exe"
Name: "{group}\{cm:UninstallProgram,HP15C}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\HP15C"; Filename: "{app}\hp15c.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\hp15c.exe"; Description: "{cm:LaunchProgram,HP15C}"; Flags: nowait postinstall skipifsilent

