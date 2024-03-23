#!/bin/bash

set -e

## run locally
PKG=$(grep "Package" DESCRIPTION | awk '{print $NF}')
BUILD_DIR=$(pwd)
DOCS_REPO=$HOME/git/wwenjie.org
TARGET_DIR=$DOCS_REPO/static/$PKG
GIT_COMMIT=$(git rev-parse --short HEAD)

## update docs by pkgdown
make pkgdown

# go to the repository for wwenjie.org
cd $DOCS_REPO
git checkout -f
git checkout main
git pull gitlab main
mkdir -p $TARGET_DIR
cp -r $BUILD_DIR/docs/* $TARGET_DIR
git add static/$PKG/
if [ -n "$(git diff --cached --exit-code)" ]
then
    git commit -m "deploy pkgdown for $PKG $GIT_COMMIT"
    git push gitlab main
else
    printf "The docs was not updated.\n"
fi
