#!/bin/bash
set -x

# this script assumes you have $TMP configured with a few GBs of space
# available

# create a tmpdir for easy cleanup
TMPDIR=$(mktemp -d -t bamliquidator-XXXXXXXX)
cd $TMPDIR

# download the data
 wget https://www.dropbox.com/s/bu75ojqr2ibkf57/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam
 wget https://www.dropbox.com/s/a71ngagu2k8pgiv/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam.bai
 wget https://www.dropbox.com/s/g7rcde76jya11y0/04032013_D1L57ACXX_4.TTAGGC.hg18.summary_chr1.tab

# run the test
bamliquidator_batch --flatten 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam

# check the results
diff -qs 04032013_D1L57ACXX_4.TTAGGC.hg18.summary_chr1.tab output/summary_chr1.tab

ECODE=$?

# cleanup 
rm -r $TMPDIR

if [ $ECODE -eq 0 ]; then
    echo -e "\033[0;32mOK: \033[0mtest produces expected results"
else
    echo -e "\033[0;31mFAIL: \033[0mtest script failed or produced unexpected results"
fi

