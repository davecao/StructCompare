REM A batch script to call java bilab-structure-1.0-jar-with-dependencies.jar

@echo off
(
set JAVA_PATH=""C:\Program Files\Java\jdk1.7.0_17\bin";"
set PATH=%PATH%;%JAVA_PATH%

set @FWDIR=%CD%
set @jarfile=%@FWDIR%\..\jars\bilab-structure-1.0-jar-with-dependencies.jar
)

java.exe -jar %@jarfile% -u
