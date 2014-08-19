#!/usr/bin/env bash

# this script is intended to be run by the makefile in this directory
# it used to be part of the makefile, but it resulted in git being a dependency of running "make"
# which was not permitted by the ppa build

git diff --quiet
if [ $? -eq 0 ]; then
CLEAN=clean
else
CLEAN=unclean
fi

set -ex

CHANGELOG="bamliquidator ($VERSION-0ppa1~%s) %s; urgency=low

  * Auto generated from makefile
  * $(git config --get remote.origin.url)
  * $(git rev-parse HEAD) ($CLEAN) 

 -- $UPLOADER <$UPLOADER_EMAIL>  $(date +"%a, %d %b %G %H:%M:%S %z")

"

mv bamliquidator-$VERSION.tar.gz bamliquidator_$VERSION.orig.tar.gz 
tar xf bamliquidator_$VERSION.orig.tar.gz
cp -R debian bamliquidator-$VERSION
cp -R bamliquidator-$VERSION bamliquidator-$VERSION-precise
mv bamliquidator-$VERSION bamliquidator-$VERSION-trusty
printf "$CHANGELOG" precise precise > bamliquidator-$VERSION-precise/debian/changelog 
printf "$CHANGELOG" trusty trusty > bamliquidator-$VERSION-trusty/debian/changelog 

pushd bamliquidator-$VERSION-precise/debian
debuild -S
popd

pushd bamliquidator-$VERSION-trusty/debian
debuild -S
popd
