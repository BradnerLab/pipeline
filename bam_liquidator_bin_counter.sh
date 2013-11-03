#!/usr/bin/env bash

##################################################################################
# The MIT License (MIT)
#
# Copyright (c) 2013 John DiMatteo 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
##################################################################################

# todo: re-write this as a C++ program, and break out the relevant bamliquidator call to a function in
#       a separate .c file.  use poor man's profiler to understand why this takes so long (I think it is 
#       disk I/O bound) and maybe speed it up

# I intend for major version to go up in increments of 100, and minor versions to go up in increments of 1.  So for 
# example if the version is 205, and a minor change is made (e.g. to speed up performance), then the version can be
# changed to 206.  However if a major change is made (e.g. to fix a bug that had resulted in prior counts being 
# incorrect), then the version should be changed to 300.  Note that a change that fixes a bug that had previously 
# prevented a bin on a chromosome from being counted may be incremented by 1 (since the prior successful counts will
# be unchanged).  This will allow me to run different versions simultaneously and/or compare performance/correctness 
# of different versions.
version=100

# todo: set from arguments
force=false
file_path="/Users/jdimatteo/DanaFarber/copied_from_tod/05012013_C22WBACXX_3.AGTTCC.hg18.bwt.sorted.bam"

# todo: set by parsing file_path
file_name="04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam"
parent_directory="copied_from_tod"

# todo: rename this
database_name=bradnerlab

# todo: don't use root mysql user

# todo: select/lock the given file (or the next file in some directory if no file specified) 
# 
# we probably want the run table to be transactional, so that we can figure out the next file and
# insert a record for it, without worrying about anyone else stealing our file


gff="HG18_100KB_UNIQUENESS.gff"
both_strands="."
one_summary=1
zero_extension=0

select_chromosomes_sql="SELECT DISTINCT chromosome FROM bins WHERE gff_name = '$gff'
                        ORDER BY LENGTH(chromosome), chromosome";
while read chromosome
do
  echo ---------------------------------------
  echo counting chromosome $chromosome - `date`

  select_bins_sql="SELECT start, end, bin_number FROM bins
                   WHERE chromosome = '$chromosome' AND gff_name = '$gff' ORDER BY bin_number;"

  while read start end bin
  do
    # todo: Ask Charles: should I subtract one from the end so that I don't count the boundaries twice?

    count=`./bamliquidator $file_path $chromosome $start $end $both_strands $one_summary $zero_extension`
    count_status=$?
    if [ $count_status -ne 0 ]; then
      echo ^^^^^ error detected in bamliquidator run ^^^^
      echo ^^^^^ count failed with return code $count_status:
      echo ^^^^^ stdout: $count
      echo ^^^^^ chromosome=$chromosome, start=$start, end=$end, bin=$bin
      echo ^^^^^ skipping the rest of this chromosome
      echo 
      # todo: record some sort of error somewhere... maybe there should be an error log table?
      # todo: Ask Charles: should I record 0 for the bins where there is no chromosome?
      break 
    fi
    #echo status=$count_status
    #echo count=$count

    insert_count_sql="INSERT DELAYED INTO COUNTS (parent_directory,file_name,chromosome,bin,count,counter_version)
                      VALUES ('$parent_directory','$file_name', '$chromosome', $bin, $count,$version);"

    # todo: don't use root mysql user
    mysql -uroot $database_name -e "$insert_count_sql"
    insert_status=$?
    if [ $insert_status -ne 0 ]; then
        # note that this may not error in most conditions since insert delayed doesn't wait for return code
        echo "insert failed: $1"
        echo sql: $insert_count_sql
        exit $insert_status
    fi
  done < <(mysql -uroot $database_name -ss -e "$select_bins_sql")

  echo counting successful 
  echo ---------------------------------------
  echo
done < <(mysql -uroot $database_name -ss -e "$select_chromosomes_sql") 

