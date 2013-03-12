# Default target of this Makefile
all::

CCC = /usr/bin/g++
STATIC = 
INCLUDES = -I /usr/include
LDFLAGS =  -lboost_program_options -lboost_system -lboost_filesystem
CCFLAGS = -O3 -Wall -Wno-deprecated

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
datarootdir = ${prefix}/share

INSTALLDIR = ${bindir}

EXEC_SRC = ukb_wsd.cc ukb_ppv.cc compile_kb.cc convert2.0.cc
EXEC = $(notdir $(basename $(EXEC_SRC)))

SOURCES= common.cc globalVars.cc configFile.cc fileElem.cc kbGraph.cc kbGraph_common.cc kbGraph_v16.cc disambGraph.cc csentence.cc csentenceSSI.cc wdict.cc

MEMBERS=$(SOURCES:.cc=.o)

.PHONY: all install clean

all:: $(EXEC)

version.h: FORCE
	@$(SHELL_PATH) ./UKB-VERSION-GEN 2.0

globalVars.o: version.h

FORCE:

%.o : %.cc
	$(CCC) -c $(CCFLAGS) -o $@ $(INCLUDES) $<

$(EXEC): $(MEMBERS) $(SRCS) $(EXEC_SRC)
	$(CCC) $(CCFLAGS) -o $@ $@.cc $(MEMBERS) $(INCLUDES) $(STATIC) $(LDFLAGS)

install: $(EXEC)
	install -d -m 755 $(INSTALLDIR)
	@for exec in $(EXEC); do \
	cp -f $$exec $(INSTALLDIR) > /dev/null 2>&1; \
	done

clean:
	find . -type f -name '*.o' | xargs rm -f
	rm -f $(EXEC)
	rm -f version.h
