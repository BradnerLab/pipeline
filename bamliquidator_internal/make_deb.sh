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

CHANGELOG="%s ($VERSION-0ppa1~%s) %s; urgency=low

  * Auto generated from makefile
  * $(git config --get remote.origin.url)
  * $(git rev-parse HEAD) ($CLEAN) 

 -- $UPLOADER <$UPLOADER_EMAIL>  $(date +"%a, %d %b %G %H:%M:%S %z")

"

py2dsc -m "$UPLOADER <$UPLOADER_EMAIL>" bamliquidatorbatch_$VERSION.orig.tar.gz  
cp python-bamliquidatorbatch.preinst deb_dist/bamliquidatorbatch-$VERSION/debian/
cp python-bamliquidatorbatch.control deb_dist/bamliquidatorbatch-$VERSION/debian/control

cp -R deb_dist/bamliquidatorbatch-$VERSION deb_dist/bamliquidatorbatch-$VERSION-precise
pushd deb_dist/bamliquidatorbatch-$VERSION-precise
printf "$CHANGELOG" bamliquidatorbatch precise precise > debian/changelog
debuild $debuild_args
popd

cp -R deb_dist/bamliquidatorbatch-$VERSION deb_dist/bamliquidatorbatch-$VERSION-trusty
pushd deb_dist/bamliquidatorbatch-$VERSION-trusty
printf "$CHANGELOG" bamliquidatorbatch trusty trusty > debian/changelog
debuild $debuild_args
popd

mv bamliquidator-$VERSION.tar.gz bamliquidator_$VERSION.orig.tar.gz 
tar xf bamliquidator_$VERSION.orig.tar.gz
cp -R debian bamliquidator-$VERSION
cp -R bamliquidator-$VERSION bamliquidator-$VERSION-precise
mv bamliquidator-$VERSION bamliquidator-$VERSION-trusty
printf "$CHANGELOG" bamliquidator precise precise > bamliquidator-$VERSION-precise/debian/changelog 
printf "$CHANGELOG" bamliquidator trusty trusty > bamliquidator-$VERSION-trusty/debian/changelog 

pushd bamliquidator-$VERSION-precise/debian
debuild $debuild_args
popd

pushd bamliquidator-$VERSION-trusty/debian
debuild $debuild_args
popd
