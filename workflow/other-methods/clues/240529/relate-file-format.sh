#!/bin/bash
relatefolder=$1
mainfolder=$2
prefix=$3
$relatefolder/bin/RelateFileFormats --mode ConvertFromVcf --haps $mainfolder/$prefix.haps --sample $mainfolder/$prefix.sample --input $mainfolder/$prefix
