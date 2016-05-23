TITLE bilab structure comparison tool
REM A batch script to call java bilab-structure-1.0-jar-with-dependencies.jar

@echo off
for /f "tokens=2 delims=:" %%a in ('systeminfo ^| find "OS Name"') do set OS_Name=%%a
for /f "tokens=* delims= " %%a in ("%OS_Name%") do set OS_Name=%%a
for /f "tokens=3 delims= " %%a in ("%OS_Name%") do set OS_Name=%%a
if "%os_name%"=="XP" set version=XP
if "%os_name%"=="7" set version=7

:CheckOS
IF EXIST "%PROGRAMFILES(X86)%" (GOTO 64BIT) ELSE (GOTO 32BIT)

:64BIT
echo 64-bit...
:: GOTO END

:32BIT
echo 32-bit...
:: GOTO END


:: -------------------------------------
:: Check Windows Version
:: 5.0 = W2K
:: 5.1 = XP
:: 5.2 = Server 2K3
:: 6.0 = Vista or Server 2K8
:: 6.1 = Win7 or Server 2K8R2
:: 6.2 = Win8 or Server 2K12
:: 6.3 = Win8.1 or Server 2K12R2
:: 0.0 = Unknown or Unable to determine
:: --------------------------------------
echo OS Detection: Windows %%version%% on "%PROCESSOR_ARCHITECTURE%"
:: ver | findstr /i "5\.0\."
:: if %ERRORLEVEL% EQU 0 (
:: echo  OS = Windows 2000
:: )
:: ver | findstr /i "5\.1\."
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Windows XP
:: )
:: ver | findstr /i "5\.2\."
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Server 2003
:: )
:: ver | findstr /i "6\.0\." > nul
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Vista / Server 2008
:: )
:: ver | findstr /i "6\.1\." > nul
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Windows 7 / Server 2008R2
:: )
:: ver | findstr /i "6\.2\." > nul
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Windows 8 / Server 2012
:: )
:: ver | findstr /i "6\.3\." > nul
:: if %ERRORLEVEL% EQU 0 (
:: echo OS = Windows 8.1 / Server 2012R2
:: )

@echo off
(
REM Set JAVA_PATH if java is not available on system path
REM set JAVA_PATH="C:\Program Files (x86)\Java\jre1.8.0_91\bin"
REM set PATH=%PATH%;%JAVA_PATH%

set @FWDIR=%~dp0\..
set @jarfile=%@FWDIR%\jars\bilab-structure-1.0-jar-with-dependencies.jar
)

java -jar %@jarfile% -u
