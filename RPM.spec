
# This spec file can be used to build an RPM package for installation.
# **NOTE**
#     Version, Release, and tagname information should be updated for the
#     particular release to build an RPM for.

# ----------------------------------------------------------------------------
Name:		espa-surface-reflectance
Version:	201607
Release:	2%{?dist}
Summary:	ESPA Surface Reflectance Software

Group:		ESPA
License:	Nasa Open Source Agreement
URL:		https://github.com/USGS-EROS/espa-surface-reflectance.git

BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
BuildArch:	x86_64
Packager:	USGS EROS LSRD

BuildRequires:	espa-product-formatter
Requires:	espa-product-formatter >= 201605


# ----------------------------------------------------------------------------
%description
Provides science application executables for generating surface reflectance products.  These applications are implementated in C and are statically built.


# ----------------------------------------------------------------------------
# Specify the repository tag/branch to clone and build from
%define tagname dev_l8sr_v3.0
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
make BUILD_STATIC=yes


# ----------------------------------------------------------------------------
%install
# Start with a clean installation location
rm -rf %{buildroot}
# Install the applications for a specific path
cd %{clonedname}
make install PREFIX=%{buildroot}/usr/local

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
/usr/local/%{name}


# ----------------------------------------------------------------------------
%changelog
* Wed Jun 15 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated for July 2016 release

* Wed Apr 13 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated for May 2016 release
* Mon Mar 07 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated release number for a recompile against a support library
* Thu Mar 03 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated release number for a recompile against a support library
* Fri Feb 12 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated release number for a bug fix Mar 2016 release
* Mon Jan 25 2016 Ronald D Dilley <rdilley@usgs.gov>
- Updated for Mar 2016 release
* Wed Dec 02 2015 Ronald D Dilley <rdilley@usgs.gov>
- Changed release number for a recompile against the product formatter for Dec 2015 release
* Wed Nov 04 2015 Ronald D Dilley <rdilley@usgs.gov>
- Updated for Dec 2015 release
