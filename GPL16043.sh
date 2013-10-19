#!/usr/bin/bash

#first make temp directory

#get the zip of the cel files and the name of the zip
CEL_ZIP=$1
NAME=$2

#CEL_ZIP=/ark/home/cl512/ressrv19/raw/expression/JB20130926st.zip
#NAME=AFFY
ID=$RANDOM
#ID=1234
TEMP_DIR_ROOT=/ark/temp/
TEMP_DIR=$TEMP_DIR_ROOT$NAME\_tmp\_$ID
INITIAL_DIR=`pwd`

#making the temp directory
mkdir $TEMP_DIR

#copy and unzip stuff
cp $CEL_ZIP $TEMP_DIR/$NAME\_tmp_$ID.zip
cd $TEMP_DIR
unzip $TEMP_DIR/$NAME\_tmp_$ID.zip

#make an analysis output directory
mkdir $TEMP_DIR/output

#run the spikey normy
R --no-save $TEMP_DIR/ $NAME < $INITIAL_DIR/GPL16043.r

#run the GPL gene level script
python $INITIAL_DIR/GPL16043.py -i $TEMP_DIR/output/$NAME\_all_mas5_probe_exprs_raw.txt
python $INITIAL_DIR/GPL16043.py -i $TEMP_DIR/output/$NAME\_all_mas5_probe_exprs_norm.txt

#zip up the output
cd $TEMP_DIR
zip -r $NAME\_output.zip output

echo "WROTE OUTPUT TO: " $TEMP_DIR/$NAME\_output.zip
