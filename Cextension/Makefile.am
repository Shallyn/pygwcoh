# Makefile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
pkgpyexecdir = ../_core

pkgpyexec_LTLIBRARIES = PyGWCOH.la

PyGWCOH_la_SOURCES = \
py_interface.c cohDatatypes.h \
cohTime.h cohTime.c \
cohUtils.h cohUtils.c \
cohDetector.h cohDetector.c \
cohCore.h cohCore.c

PyGWCOH_la_LIBADD = -lgsl -lgslcblas
PyGWCOH_la_LDFLAGS = -module -avoid-version

all: $(pkgpyexec_LTLIBRARIES:.la=.so)

.la.so:
	ln -s .libs/$@ $@ || true
