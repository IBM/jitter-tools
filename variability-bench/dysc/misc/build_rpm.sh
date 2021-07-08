#!/bin/bash
blddir=~/rpmbuild

rm -Ir ${blddir}

mkdir -p ${blddir}/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
mkdir -p ${blddir}/SOURCES/jitter-tools
cp -a ../../* ${blddir}/SOURCES/jitter-tools/.
cd ${blddir}/SOURCES
#git clone git@github.ibm.com:dcs-work/jitter-tools.git
tar czf jitter-tools.tar.gz jitter-tools
cd ../..

cd ${blddir}/SPECS
ln -sf ../SOURCES/jitter-tools/dysc/misc/dysc.spec
cd ../..

sed -i "s|%define _topdir.*|%define _topdir ${blddir}|g" ${blddir}/SPECS/dysc.spec

rpmbuild -v -ba --clean ${blddir}/SPECS/dysc.spec
