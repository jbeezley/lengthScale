F90=gfortran
NETCDF=/usr/local
UMFPACK=/usr/local
DEBUG=-O3
FFLAGS=$(DEBUG) -I$(NETCDF)/include
LFLAGS=$(DEBUG) `${NETCDF}/bin/nf-config --flibs`
CFLAGS=$(DEBUG) -I$(UMFPACK)/include
%.o: %.F90
	$(F90) -c $(FFLAGS) $< -o $@

%.o: %.c
	$(F90) -c $(CFLAGS) $< -o $@

%.o: %.f
	$(F90) -c $(FFLAGS) $< -o $@

OBJS= pointMorphing.o morphing.o pointPert.o lengthScale.o solver.o \
	  correlation.o derivative.o fileio.o sparse.o image_interp.o   \
	  common.o
#csparse.o coocsr.o
#umf4_f77wrapper.o coocsr.o
test: testPointPert.x testLengthScale.x
clean:
	rm -f *.o *.x *.mod *.a

umf4_f77wrapper.o: umf4_f77wrapper.c
common.o: common.F90
fileio.o: fileio.F90 sparse.o common.o
pointPert.o: pointPert.F90 common.o
image_interp.o: image_interp.F90
morphing.o: morphing.F90 common.o image_interp.o correlation.o derivative.o
derivative.o: derivative.F90 common.o
correlation.o: correlation.F90 common.o
sparse.o: sparse.F90 common.o
lengthScale.o: lengthScale.F90 common.o pointMorphing.o correlation.o sparse.o
solver.o: solver.F90 common.o sparse.o
csparse.o: csparse.c
pointMorphing.o: pointMorphing.F90 pointPert.o morphing.o common.o
testPointPert.o: testPointPert.F90 pointPert.o fileio.o common.o
testPointPert.x: testPointPert.o pointPert.o fileio.o common.o
	$(F90) $^ $(LFLAGS) -o $@
testLengthScale.o: testLengthScale.F90 $(OBJS)
testLengthScale.x: testLengthScale.o $(OBJS)
	$(F90) $< $(OBJS) $(LFLAGS) -o $@ 

libpspline.a: 
	(cd pspline ; FC=$(F90) FC90=$(F90) NETCDF=$(NETCDF) make)
	ln -f pspline/libpspline.a .
	ln -f pspline/*.mod .

