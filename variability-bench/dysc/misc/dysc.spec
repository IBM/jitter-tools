
%define _topdir     /home/psiegl/rpmbuild

Name:           dysc
Version:        0.1
Release:        0
Summary:        IBM DYnamic System Cgroup daemon

Group:          Development/Tools
License:        IBM
URL:            https://github.ibm.com/dcs-work/jitter-tools/tree/master/dysc
Source0:        jitter-tools.tar.gz

BuildRequires:  gcc-c++

%description
The IBM DYnamic System Cgroup daemon increases and decreases the system cgroup on demand.

%prep
%setup -q -n jitter-tools

cd %{_builddir}/jitter-tools/dysc && sh bld.sh

%install
mkdir -p %{buildroot}/%{_sbindir}
cp %{_builddir}/jitter-tools/dysc/dysc %{buildroot}/%{_sbindir}/dysc

mkdir -p %{buildroot}/%{_exec_prefix}/lib/systemd/system/
cp %{_builddir}/jitter-tools/dysc/misc/dysc.service %{buildroot}/%{_exec_prefix}/lib/systemd/system/dysc.service

mkdir -p %{buildroot}/%{_sysconfdir}/sysconfig
cp %{_builddir}/jitter-tools/dysc/misc/dysc.cfg %{buildroot}/%{_sysconfdir}/sysconfig/dysc

%clean
rm -rf %{buildroot}

%files
%{_sbindir}/dysc
%{_exec_prefix}/lib/systemd/system/dysc.service
%{_sysconfdir}/sysconfig/dysc
