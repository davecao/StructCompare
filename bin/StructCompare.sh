#!/usr/bin/env bash

FWDIR=$(cd `dirname $0`;pwd)

exec_jar=$FWDIR/../jars/bilab-structure-1.0-jar-with-dependencies.jar
java -jar $exec_jar -u
