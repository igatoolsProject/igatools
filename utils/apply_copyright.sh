#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
#
# This file is part of the igatools library.
#
# The igatools library is free software: you can use it, redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-+--------------------------------------------------------------------

#!/usr/bin/bash
#
# Script using sed to apply the copyright header 
# written in the file utils/copyright_header.txt
# to all igatools files.
# 
# This script should be run from the root of igatools
# sh ./utils/apply_copyright.sh
#

# Function to apply the copyrigth
# copyright (comment_character, copyright text, files_to_apply_copyright) 
function copyright {
    com_char=$1
    copyrigth_file=$2
    
    # prepare a temporary header copyright for the give type of file
    hline="-+--------------------------------------------------------------------"
    header_file=tmp_copyright_header.txt
    echo "$com_char$hline" > $header_file
    while read -r line
    do
	if [[ -z "$line" ]]; then
	    echo "$com_char" >> $header_file
	else
    	    echo "$com_char $line" >> $header_file
	fi
    done < $copyrigth_file
    echo "$com_char$hline" >> $header_file

    
    files=$3
    # echo $files
    for file in $files
    do
	last_line=`grep -n $com_char-+ $file |tail -1 |cut -d: -f1`
	if [[ -z "$last_line" ]]; then
	    echo "No copyright found in $file"
	    cp  $file $file.1
	else
    	    sed  '1,'$last_line'd' <$file >$file.1 
	fi
	
    	cat $header_file $file.1 > $file
	
    	rm $file.1   
    done
    
    rm $header_file
}

files=`find . -name "*.h"`
files="$files `find . -name "*.cpp"`"
files="$files `find . -name "*.h.in"`"
copyright "//" "./utils/copyright_header.txt" "$files"

files="`find . -name "*.py"`  ./tutorial/cmake_template.txt"
files="$files `find ./utils -name "*.sh"`"
files="$files `find . -name "*.cmake"`"
files="$files `find . -name "CMakeLists.txt"`"
copyright "#" "./utils/copyright_header.txt" "$files"

files=`find ./utils/logo -name "*.tex"`
copyright "%" "./utils/copyright_header.txt" "$files"

# apply copyright to dealii copied files
files="include/igatools/base/exceptions.h 
include/igatools/base/logstream.h 
source/base/exceptions.cpp 
source/base/logstream.cpp 
./source/contrib/table_handler.cpp 
./include/igatools/contrib/table_indices.h 
./include/igatools/contrib/table.h 
./include/igatools/contrib/table_handler.h"
copyright "//" "./utils/dealii_copyright_header.txt" "$files"
