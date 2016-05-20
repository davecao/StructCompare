REM A batch script to call java bilab-structure-1.0-jar-with-dependencies.jar

@echo off
(
REM Set JAVA_PATH if java is not available on system path
REM set JAVA_PATH="C:\Program Files (x86)\Java\jre1.8.0_91\bin"
REM set PATH=%PATH%;%JAVA_PATH%

set @FWDIR=%cd%
set @jarfile=%@FWDIR%\..\jars\bilab-structure-1.0-jar-with-dependencies.jar
)

java.exe -jar %@jarfile% -u
