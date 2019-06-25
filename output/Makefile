# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

RC     := root-config
ifeq ($(findstring $(MAKECMDTARGET),clean),)
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
RCONFIG := $(wildcard $(shell $(RC) --incdir)/RConfigure.h)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
include $(ROOTSYS)/etc/Makefile.arch
endif
endif

#------------------------------------------------------------------------------

PRODUCTIONO       = production.$(ObjSuf) productionDict.$(ObjSuf)
PRODUCTIONS       = production.$(SrcSuf) productionDict.$(SrcSuf)
PRODUCTIONSO      = lib_production.$(DllSuf)

OBJS         = $(PRODUCTIONO) 
PROGRAMS     = $(PRODUCTIONSO) 

ifeq ($(ARCH),aix5)
MAKESHARED    = /usr/vacpp/bin/makeC++SharedLib
endif

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


ifeq ($(ARCH),aix5)
		$(MAKESHARED) $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif

$(PRODUCTIONSO):    $(PRODUCTIONO)
ifeq ($(ARCH),aix5)
		$(MAKESHARED) $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif


clean:
		@rm -f $(OBJS) $(TRACKMATHSRC) core *Dict.* cernstaff.root sg*


.SUFFIXES: .$(SrcSuf)

###
production.$(ObjSuf): production.h 
productionDict.$(SrcSuf): production.h 
	@echo "Generating dictionary $@..."
	$(ROOTCLING) -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

ifeq ($(GCC_MAJOR),4)
ifeq ($(GCC_MINOR),1)
productionDict.o: CXXFLAGS += -Wno-strict-aliasing
endif
endif
