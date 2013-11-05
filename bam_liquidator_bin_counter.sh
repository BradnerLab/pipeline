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

usage_string="usage: $0 [options] path 

Records bamliquidator count on each bin in the given bam file.  If the path is
a directory instead of a bam file, then counts the next uncounted file in the
provided directory (can be recursive).

The counts are stored in a mysql database.

OPTIONS:
   -c      Check results against baseline version (non-zero rcode if any checks fail)
"

usage()
{
  echo "$usage_string" 
}

baseline_check=0

while getopts “:c” OPTION
do
  case $OPTION in
    c)
      baseline_check=1 
    ;;
    ?)
      usage
      exit 1
    ;;
  esac
done

shift $((OPTIND-1))
if [ -z "$@" ]; then
  usage
  exit 1
fi

path="$@"

if [ -d "$path" ]; then
  echo "todo: support directory as argument"
  exit 1
elif [ -f "$path" ]; then
  file_path="$path"
else
  echo "$path is not a valid path"
  usage
  exit 1
fi

echo  starting with parameters:
echo "  baseline_check: $baseline_check"
echo "  file_path: $file_path"

# I intend for major version to go up in increments of 100, and minor versions to go up in increments of 1.  So for 
# example if the version is 205, and a minor change is made (e.g. to speed up performance), then the version can be
# changed to 206.  However if a major change is made (e.g. to fix a bug that had resulted in prior counts being 
# incorrect), then the version should be changed to 300.  Note that a change that fixes a bug that had previously 
# prevented a bin on a chromosome from being counted may be incremented by 1 (since the prior successful counts will
# be unchanged).  This will allow me to run different versions simultaneously and/or compare performance/correctness 
# of different versions.
version=200
baseline_version=100 # used to verify counts haven't changed if baseline checking is enabled

file_name=`basename $file_path`
parent_directory=$(basename $(dirname $file_path))

database_name=meta_analysis
mysql_user=counter

# todo: select/lock the given file (or the next file in some directory if no file specified)
# 
# we probably want the run table to be transactional, so that we can figure out the next file and
# insert a record for it, without worrying about anyone else stealing our file

gff="HG18_100KB_UNIQUENESS.gff"
both_strands="."
one_summary=1
zero_extension=0

select_chromosomes_sql="SELECT DISTINCT chromosome FROM bins
                         WHERE gff_name = '$gff'
                           AND skip = false 
                         ORDER BY LENGTH(chromosome), chromosome";

return_code=0

while read chromosome
do
  echo --------------------------------------------------------------------------------
  echo counting chromosome $chromosome - `date`
  echo

  select_bins_sql="SELECT start, end - 1, bin_number FROM bins
                   WHERE chromosome = '$chromosome' AND gff_name = '$gff' ORDER BY bin_number;"

  while read start end bin
  do
    count=`./bamliquidator $file_path $chromosome $start $end $both_strands $one_summary $zero_extension`
    count_status=$?
    if [ $count_status -ne 0 ]; then
      echo error detected in bamliquidator run
      echo "   count failed with return code $count_status:"
      echo "   stdout: $count"
      echo "   chromosome=$chromosome, start=$start, end=$end, bin=$bin"
      echo "   skipping the rest of this chromosome"
      echo 
      # todo: record some sort of error somewhere... maybe there should be an error log table?
      break 
    fi
    #echo status=$count_status
    #echo count=$count

    insert_count_sql="INSERT DELAYED INTO COUNTS (parent_directory,file_name,chromosome,bin,count,counter_version)
                      VALUES ('$parent_directory','$file_name', '$chromosome', $bin, $count, $version);"

    mysql -u$mysql_user $database_name -e "$insert_count_sql"
    insert_status=$?
    if [ $insert_status -ne 0 ]; then
        # note that this may not error in most conditions since insert delayed doesn't wait for return code
        echo "insert failed: $1"
        echo sql: $insert_count_sql
        exit $insert_status
    fi

    if [ $baseline_check -ne 0 ]; then
      select_baseline_count_sql="SELECT count FROM counts
                                  WHERE counter_version=$baseline_version
                                    AND bin=$bin AND chromosome='$chromosome'
                                    AND file_name='$file_name'"

      baseline_count=`mysql -u$mysql_user $database_name -ss -e "$select_baseline_count_sql"`
      if [ $count -ne $baseline_count ]; then
        echo baseline check: current count does not match baseline
        echo "   baseline count: $baseline_count"
        echo "   current count:  $count"
        echo "   baseline select sql: " $select_baseline_count_sql
        echo
        ((return_code++))
      fi
    fi

  done < <(mysql -u$mysql_user $database_name -ss -e "$select_bins_sql")

done < <(mysql -u$mysql_user $database_name -ss -e "$select_chromosomes_sql") 


if [ $baseline_check -ne 0 ]; then
  select_baseline_number_of_counts_sql="SELECT count(*) FROM counts WHERE counter_version=$baseline_version;"
  select_current_number_of_counts_sql="SELECT count(*) FROM counts WHERE counter_version=$version;"
  
  baseline_number_of_counts=`mysql -u$mysql_user $database_name -ss -e "$select_baseline_number_of_counts_sql"`
  current_number_of_counts=`mysql -u$mysql_user $database_name -ss -e "$select_current_number_of_counts_sql"`

  if [ $baseline_number_of_counts -ne $current_number_of_counts ]; then
      echo baseline check: current count does not match baseline
      echo "   baseline number of counts: $baseline_number_of_counts"
      echo "   current number of counts:  $current_number_of_counts"
      echo "   baseline select number of counts sql: " $select_baseline_number_of_counts_sql
      echo
      ((return_code++))
  fi
fi

exit $return_code
