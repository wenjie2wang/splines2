#!/bin/bash

# Note: this script is should be sourced from the project root directory

set -e

if [ "$(uname)" != "Linux" ]; then
    printf "Remeber to update date and version number.\n"
else
    printf "Updating date, version, and copyright year.\n"

    # define some variables
    yr=$(date +%Y)
    dt=$(date +%Y-%m-%d)
    cprt_R=misc/copyright.R
    cprt_cpp=misc/copyright.cpp
    citation=inst/CITATION
    version=$(grep "Version" DESCRIPTION | awk '{print $NF}')

    # update copyright year in the template headers
    regexp1="s/Copyright \(C\) 2016-[0-9]+/Copyright \(C\) 2016-$yr/"
    sed -i -E "$regexp1" $cprt_R
    sed "s_#_/_g" $cprt_R > $cprt_cpp

    # update copyright year in all R scripts
    for Rfile in R/*.R
    do
        if ! grep -q 'Copyright (C)' $Rfile; then
            cat $cprt_R $Rfile > .tmp
            mv .tmp $Rfile
        fi
        sed -i -E "$regexp1" $Rfile
    done

    # update copyright year in all C++ scripts
    for cppfile in src/*.cpp inst/include/*.h inst/include/*/*.h
    do
        if ! grep -q 'Copyright (C)' $cppfile; then
            cat $cprt_cpp $cppfile > .tmp
            mv .tmp $cppfile
        fi
        sed -i -E "$regexp1" $cppfile
    done
    rm $cprt_cpp

    # update date in DESCRIPTION
    regexp2="s/Date: [0-9]{4}-[0-9]{1,2}-[0-9]{1,2}/Date: $dt/"
    sed -i -E "$regexp2" DESCRIPTION

    # update version and year in citation
    regexp3="s/version ([0-9]+\.*)+/version $version/"
    sed -i -E "$regexp3" $citation
    # restrict the search and only update the year of package
    regexp4="/splines2-package/,/^\)$/ s/20[0-9]{2}/$yr/"
    sed -i -E "$regexp4" $citation

    # done
    printf "All updated.\n"
fi
