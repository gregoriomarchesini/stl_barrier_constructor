MSKMAJORVER=11
MSKMINORVER=0
CXX=c++
LD=c++
AR=ar
STD=c++11

CFLAGS_DEBUG=-g
CFLAGS_RELEASE=-O2

CFLAGS=-O2 -fPIC -std=$(STD)
IPATHS=-I../../h -I.
LFLAGS=-std=$(STD)

libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER): fusion.os BaseModel.os Debug.os mosektask.os IntMap.os SolverInfo.os StringBuffer.os
	$(LD) -shared -L../../bin -Wl,-rpath-link,../../bin '-Wl,-rpath,$$ORIGIN' -o libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER) $(LFLAGS) fusion.os mosektask.os BaseModel.os Debug.os IntMap.os SolverInfo.os StringBuffer.os -lmosek64

libfusion64.a:      fusion.os BaseModel.os Debug.os mosektask.os IntMap.os SolverInfo.os StringBuffer.os
	$(AR) r libfusion64.a fusion.os BaseModel.os Debug.os mosektask.os IntMap.os SolverInfo.os StringBuffer.os

%.os: %.cc
	$(CXX) $(IPATHS) $(CFLAGS) -c -o $@ $<

.PHONY: clean install

install: ../../bin/libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER) ../../bin/libfusion64.so

../../bin/libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER): libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER) libfusion64.a
	install libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER) ../../bin
	install libfusion64.a      ../../bin
../../bin/libfusion64.so: ../../bin/libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER)
	rm -f ../../bin/libfusion64.so && ln -s ./libfusion64.so.$(MSKMAJORVER).$(MSKMINORVER) ../../bin/libfusion64.so

clean:
	rm -rf *.os
