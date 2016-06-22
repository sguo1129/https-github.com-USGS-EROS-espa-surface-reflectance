
# This spec file can be used to build an RPM package for installation.
# **NOTE**
#     Version, Release, and tagname information should be updated for the
#     particular release to build an RPM for.


%define project espa-surface-reflectance
%define algorithm ledaps-aux
%define build_timestamp %(date +"%%Y%%m%%d%%H%%M%%S")
# Specify the repository tag/branch to clone and build from
%define tagname dev_aug2016
# Specify the name of the directory to clone into
%define clonedname %{name}-%{tagname}
# Change the default rpm name format for the rpm built by this spec file
%define _build_name_fmt %%{NAME}.%%{VERSION}.%%{RELEASE}%{?dist}.%{ARCH}.rpm


# ----------------------------------------------------------------------------
Name:		%{project}-%{algorithm}
Version:	2.7.0
Release:	1.%{build_timestamp}
Summary:	ESPA Surface Reflectance Auxiliary Software - LEDAPS

Group:		ESPA
License:	NASA Open Source Agreement
URL:		https://github.com/USGS-EROS/espa-surface-reflectance.git

BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
BuildArch:	x86_64
Packager:	USGS EROS LSRD

BuildRequires:	espa-product-formatter >= 1.8.0

%description
Provides executables for retrieveing Ledaps specific auxiliary data.


# ----------------------------------------------------------------------------
%prep
# We don't need to perform anything here

%build
# Start with a clean clone of the repo
rm -rf %{clonedname}
git clone --depth 1 --branch %{tagname} %{url} %{clonedname}
# Build the applications
cd %{clonedname}
make all-ledaps-aux BUILD_STATIC=yes

%install
# Start with a clean installation location
rm -rf %{buildroot}
# Install the applications for a specific path
cd %{clonedname}
make install-ledaps-aux PREFIX=%{buildroot}/usr/local

%clean
# Cleanup our cloned repository
rm -rf %{clonedname}
# Cleanup our installation location
rm -rf %{buildroot}


# ----------------------------------------------------------------------------
%files
%defattr(-,root,root,-)
# All sub-directories are automatically included
/usr/local/bin/*
/usr/local/%{project}/%{algorithm}


# ----------------------------------------------------------------------------
%changelog
* Tue Jun 21 2016 Ronald D Dilley <rdilley@usgs.gov>
- Initial Version for August 2016 release
