
# This spec file can be used to build an RPM package for installation.
# **NOTE**
#     Version, Release, and tagname information should be updated for the
#     particular release to build an RPM for.

# ----------------------------------------------------------------------------
# Disable the creation of debug-info RPM package
# that contains stripped binary debug symbols
%global _enable_debug_package 0
%global debug_package %{nil}
# Disable stripped binary generation
%define __os_install_post %{nil}

# ----------------------------------------------------------------------------
# Change the default rpm name format for the rpm built by this spec file
%define _build_name_fmt %%{NAME}-aux.%%{VERSION}.%%{RELEASE}.rpm

Name:		espa-surface-reflectance
Version:	201512
Release:	1%{?dist}
Summary:	ESPA Surface Reflectance Auxiliary Software

Group:		ESPA
License:	Nasa Open Source Agreement
URL:		https://github.com/USGS-EROS/espa-surface-reflectance.git

BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-aux-%{version}-%{release}-XXXXXX)
BuildArch:	x86_64
Packager:	USGS EROS LSRD

BuildRequires:	espa-common
Requires:	espa-common >= 1.5.0

# ----------------------------------------------------------------------------
%description
Provides science application executables for generating surface reflectance products.  This is a C implementation which is statically built.


# ----------------------------------------------------------------------------
# Specify the repository tag/branch to clone and build from
#%define tagname dev_nov2015
%define tagname updated-makefiles
# Specify the name of the directory to clone into
%define clonedname %{name}-%{tagname}


# ----------------------------------------------------------------------------
%prep
# We don't need to perform anything here


# ----------------------------------------------------------------------------
%build

# Start with a clean clone of the repo
rm -rf %{clonedname}
git clone --depth 1 --branch %{tagname} %{url} %{clonedname}
# Build the applications
cd %{clonedname}
make all-aux BUILD_STATIC=yes


# ----------------------------------------------------------------------------
%install
# Start with a clean installation location
rm -rf %{buildroot}
# Install the applications for a specific path
cd %{clonedname}
make install-aux PREFIX=%{buildroot}/usr/local

# ----------------------------------------------------------------------------
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
/usr/local/%{name}/l8_sr/bin/combine_l8_aux_data
/usr/local/%{name}/l8_sr/bin/updatelads.py
/usr/local/%{name}/ledaps/bin/convert_ozone
/usr/local/%{name}/ledaps/bin/ncep_repackage
/usr/local/%{name}/ledaps/bin/updatencep.py
/usr/local/%{name}/ledaps/bin/updatetoms.py


# ----------------------------------------------------------------------------
%changelog
* Wed Nov 04 2015 Ronald D Dilley <rdilley@usgs.gov>
- Build for Dec 2015 release
- Initial implementation
