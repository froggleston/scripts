#!/bin/bash

## Generate JSON output from an ENA accession
## Author: Rob Davey, The Genome Analysis Centre (TGAC), UK
## http://www.ebi.ac.uk/ena/data/warehouse/usage

## project, study, sample, experiment, or run accession
[ $# -ge 1 ] && PROJ="$1" || read PROJ

## returned result type. defaults to "read_run"
RESULT=$2

if [ -z "$PROJ" ] ; then
  echo "No accession supplied"
  exit 1
fi

if [ -z "$2" ] ; then
  RESULT="read_run"
fi

OUT=`curl --silent "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$PROJ&download=text&result=$RESULT"`

i=1
HEADERS=

echo "["

while read -r line ; do
  ## substitute tabs for commas. IFS doesn't do nice multi-whitespace separation.
  line=${line//$'\t'/,}

  if [ $i -eq 1 ] ; then
    # parse headers
    IFS=$',' read -r -a HEADERS <<< "$line"
  else
    # parse values
    echo "{"
    while IFS=$',' read -r -a VALUES ; do
      for j in "${!HEADERS[@]}" ; do
        echo "\"${HEADERS[j]}\":\"${VALUES[j]}\","
      done
    done <<< "$line" | sed '$s/,//'
    echo "},"
  fi
  i=$((i + 1));
done <<< "$OUT" | sed '$s/,//'

echo "]"
