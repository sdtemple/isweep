#!/bin/bash

# Zero step: make folder structure
# Seth D. Temple, sdtemple@uw.edu
# April 26, 2023

WHERE=$1 # name of study analysis
LOGS=$2 # folder to put -o and -e from jobs in

mkdir -p $WHERE
mkdir -p ${WHERE}/vcfs
mkdir -p ${WHERE}/maps
mkdir -p ${WHERE}/plots
mkdir -p ${WHERE}/ibdsegs
mkdir -p ${WHERE}/ibdsegs/hapibd
mkdir -p ${WHERE}/ibdsegs/ibdends
mkdir -p ${WHERE}/ibdsegs/ibdends/modified
mkdir -p ${WHERE}/ibdsegs/ibdends/modified/mom
mkdir -p ${WHERE}/ibdsegs/ibdends/modified/scan
mkdir -p ${LOGS}
