# Makefile for minfit4; 
#
# toggle comments on library definitions below to compile 
SLATEC=../../slatec/libslatec.a
LAPACK=../../lapack/lapack-3.2.1/lapack_LINUX.a 
TLIB=../../tlib/tool.a
MFLIB=libminfit.a
MGOLIBS=-L/usr/X11R6/lib -lmongo -lX11 -lreadline -lhistory -lcurses
OPTIONS = -O
#	This is the list of object modules inside the library that we want
#	make to keep up-to-date. 
LOBJS=$(MFLIB)(getatomic.o)  $(MFLIB)(getlines.o)  $(MFLIB)(getparms.o)  \
      $(MFLIB)(getregions.o) $(MFLIB)(readdata.o)  $(MFLIB)(strappend.o) \
      $(MFLIB)(fcn.o)        $(MFLIB)(inform.o)    $(MFLIB)(getlist.o)   \
      $(MFLIB)(results.o)    $(MFLIB)(voigt.o)     $(MFLIB)(defcoeffs.o) \
      $(MFLIB)(errors.o)     $(MFLIB)(fileowc.o)   $(MFLIB)(mgosubs.o)   \
      $(MFLIB)(instrument.o) $(MFLIB)(convolve.o)  $(MFLIB)(fft.o)       \
      $(MFLIB)(fitregion.o)  $(MFLIB)(features.o)  $(MFLIB)(adjustfit.o) \
      $(MFLIB)(convpix.o)    $(MFLIB)(getlimits.o) $(MFLIB)(ftest.o)
FPP=f
FC=f95
.SUFFIXES:
.SUFFIXES: .out .o .c .e .r .F .f .y .l .s .p
.F.o:
	/lib/cpp $(CPPFLAGS) $*.F $*.$(FPP)
	$(FC) -c $(FFLAGS) $*.$(FPP)
	rm -f $*.$(FPP)
.F.a:
	/lib/cpp $(CPPFLAGS) $< $*.$(FPP)
	$(FC) -c $(FFLAGS) $*.$(FPP)
	$(AR) rv $@ $*.o
	rm -f $*.$(FPP) $*.o
.f.a:
	$(FC) -c $(OPTIONS) $(FFLAGS) $*.f
	$(AR) rv $@ $*.o
	rm -f $*.o
#########
#
default:: minfit
#
#       Don't let make delete the tools object library!
.PRECIOUS: $(MFLIB)
#
#	Tools library depends upon the object modules
#	to make tools library, use ar and ranlib
lib: $(LOBJS) minfit.h
	ranlib $(MFLIB)
minfit: libminfit.a minfit.f access.f
	$(FC) $(OPTIONS) minfit.f access.f -o minfit $(MFLIB) $(SLATEC) $(MGOLIBS) $(LAPACK)
#
#	remove all binaries and other typical garbage files
clean:
	-rm -f *.o *.a a.out *.trace \#* *.bak *~
tidy:
	-rm -f a.out *.trace \#* *.bak
