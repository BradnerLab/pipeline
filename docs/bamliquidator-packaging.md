# Bamliquidator Packaging

### Overview

* `make dput` needs to be run from a proper environment to update the bamliquidator ppa version
    * prerequisites in addition to [developer setup](https://github.com/BradnerLab/pipeline/wiki/bamliquidator#check-list) on Ubuntu 18.04: `sudo apt install devscripts python-all python-stdeb dput`
* the Bradner Lab main EC2 instance is configured for packaging with the user `packager`
    * you can use the ubuntu ssh key to as packager
    * if you need the gpg passphrase, contact jdimatteo@gmail.com 
* Ubuntu Launchpad account name is dfbradnerlab@gmail.com
    * Upload errors and acknowledgments are emailed to dfbradnerlab@gmail.com
* it takes a few minutes for new versions to be published
    * it takes a couple minutes for the upload acknowledgment to be emailed
    * next, it takes a couple more minutes to go from status "Pending" to "Published" at https://launchpad.net/~bradner-computation/+archive/ubuntu/pipeline/+packages

### Example Updating the Package Version

```bash
$ ssh packager@ec2
$ cd pipeline/bamliquidator_internal/
$ git pull
$ make dput
```

### Example Testing Package Before Uploading

```bash
$ ssh packager@ec2
$ cd pipeline/bamliquidator_internal/
$ git pull
$ make deb
$ scp deb_dist/python-bamliquidatorbatch_*trusty*.deb bamliquidator_*trusty*.deb ubuntu@some-other-ec2-machine:/home/ubuntu
```

```bash
$ ssh ubuntu@some-other-ec2-machine
$ sudo apt-get install gdebi-core
$ sudo gdebi python-bamliquidatorbatch*.deb
$ sudo gdebi bamliquidator*.deb
$ # do testing
$ sudo dpkg --remove bamliquidator python-bamliquidatorbatch
```

### Environment Setup Notes

1. Follow the steps on https://github.com/BradnerLab/pipeline/wiki/bamliquidator#check-list
2. `sudo apt-get install devscripts debhelper gnupg-agent pinentry-curses python-stdeb python-all-dev`
3. generate a gpg key, e.g. `gpg --gen-key`
4. register the gpg key, e.g. see https://help.ubuntu.com/community/GnuPrivacyGuardHowto#Uploading_the_key_to_Ubuntu_keyserver
5. import the gpg key at https://launchpad.net/~bradner-computation/+editpgpkeys
    * you will get an email "Launchpad: Confirm your OpenPGP Key"
    * copy from "-----BEGIN PGP MESSAGE-----" to "-----END PGP MESSAGE-----" (including the BEGIN/END PGP MESSAGE lines) to a file doc.encrypted on the system with the gpg key
    * run `gpg --decrypt doc.encrypted` and enter your passphrase when prompted
    * follow the instructions in the decrypted message
6. to reduce number of times entering the gpg passphrase, configure gpg-agent (e.g. see http://unix.stackexchange.com/questions/46960/how-to-configure-gpg-to-enter-passphrase-only-once-per-session)
