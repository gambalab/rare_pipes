#!/bin/bash

mkdir tmp_dbsnfp
cd tmp_dbsnfp

# Download
aria2c -x 16 -s 16 https://dbnsfp.s3.amazonaws.com/dbNSFP4.9a.zip

# Uncompress
unzip dbNSFP4.9a.zip
rm dbNSFP4.9a.zip


# Create a single file version
(zcat dbNSFP4.9a_variant.chr1.gz | head -n 1 | cut -f1,2,3,4,83,122,125,138,285,362,440,441,442,443,444,445,446,447,448 ; zcat dbNSFP4.9a_variant.chr*.gz | cut -f1,2,3,4,83,122,125,138,285,362,440,441,442,443,444,445,446,447,448 | grep -v "^#" ) > dbNSFP4.9_small.txt

# Compress using block-gzip algorithm
bgzip -@ 16 dbNSFP4.9_small.txt

# Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFP4.9_small.txt.gz

cd ..
mkdir data
mv tmp_dbsnfp/dbNSFP4.9_small.txt.gz data/
mv tmp_dbsnfp/dbNSFP4.9_small.txt.gz.tbi data/

rm -rf tmp_dbsnfp