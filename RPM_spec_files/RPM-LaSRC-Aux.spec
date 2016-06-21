
# This spec file can be used to build an RPM package for installation.
# **NOTE**
#     Version, Release, and tagname information should be updated for the
#     particular release to build an RPM for.

# ----------------------------------------------------------------------------

%define project espa-surface-reflectance
%define algorithm lasrc-aux
%define build_timestamp %(date +"%%Y%%m%%d%%H%%M%%S")

# ----------------------------------------------------------------------------
# Change the default rpm name format for the rpm built by this spec file
%define _build_name_fmt %%{NAME}.%%{VERSION}.%%{RELEASE}%{?dist}.%{ARCH}.rpm

Name:		%{project}-%{algorithm}
Version:	0.6.1
Release:	1.%{build_timestamp}
Summary:	ESPA Landsat Surface Reflectance Code Auxiliary Software

Group:		ESPA
License:	Nasa Open Source Agreement
URL:		https://github.com/USGS-EROS/espa-surface-reflectance.git

BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
BuildArch:	x86_64
Packager:	USGS EROS LSRD

#BuildRequires:	espa-product-formatter
#Requires:	espa-product-formatter >= 1.6.0

# ----------------------------------------------------------------------------
%description
Provides application executables for building the LaSRC auxiliary data.  This is a C implementation which is statically built.


# ----------------------------------------------------------------------------
# Specify the repository tag/branch to clone and build from
%define tagname dev_hotfix_l8aux
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
make all-l8-sr-aux BUILD_STATIC=yes


# ----------------------------------------------------------------------------
%install
# Start with a clean installation location
rm -rf %{buildroot}
# Install the applications for a specific path
cd %{clonedname}
make install-l8-sr-aux PREFIX=%{buildroot}/usr/local

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
/usr/local/%{project}/l8_sr/bin/combine_l8_aux_data
/usr/local/%{project}/l8_sr/bin/updatelads.py


# ----------------------------------------------------------------------------
%changelog
* Tue Jun 07 2016 Ronald D Dilley <rdilley@usgs.gov>
- Build for 0.6.1 release
- Initial implementation
