#!/bin/bash
mainfolder=$1
vcffile=$2
tabix -fp vcf $mainfolder/$vcffile
