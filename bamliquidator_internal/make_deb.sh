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

for ubuntu_version in "trusty" "xenial" "bionic"
do
  cp -R deb_dist/bamliquidatorbatch-$VERSION deb_dist/bamliquidatorbatch-$VERSION-$ubuntu_version
  pushd deb_dist/bamliquidatorbatch-$VERSION-$ubuntu_version
  printf "$CHANGELOG" bamliquidatorbatch $ubuntu_version $ubuntu_version > debian/changelog
  debuild $debuild_args
  popd
done

mv bamliquidator-$VERSION.tar.gz bamliquidator_$VERSION.orig.tar.gz 
tar xf bamliquidator_$VERSION.orig.tar.gz

cp -R debian bamliquidator-$VERSION
for ubuntu_version in "trusty" "xenial" "bionic"
do
  cp -R bamliquidator-$VERSION bamliquidator-$VERSION-$ubuntu_version
  printf "$CHANGELOG" bamliquidator $ubuntu_version $ubuntu_version > bamliquidator-$VERSION-$ubuntu_version/debian/changelog
  pushd bamliquidator-$VERSION-$ubuntu_version/debian
  debuild $debuild_args
  popd
done
rm -rf bamliquidator-$VERSION
