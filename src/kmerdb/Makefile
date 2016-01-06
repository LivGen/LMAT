include ../../Makefile.inc

CXXHDR = Utils.hpp KmerFileMetaData.hpp KmerNode.hpp \
         KmerFileMetaData.hpp TaxTree.hpp TaxNode.hpp  \
         TaxTable.hpp SortedDb.hpp

CXXSRCS = Utils.cpp KmerFileMetaData.cpp TaxTree.cpp TaxNode.cpp TaxTable.cpp SortedDb.hpp

CXXOBJS = Utils.o KmerFileMetaData.o TaxTree.o TaxNode.o TaxTable.o SortedDb.o

lib: libmetag.a


libmetag.a: $(CXXHDR) $(CXXSRCS) $(CXXOBJS)
	ar -cvq libmetag.a $(CXXOBJS)
	cp libmetag.a lib

 #clean $(CXXHDR) $(CXXSRCS) $(CXXOBJS)
clean : 
	rm -f *.o tmp *.mc a.out *TVD* *.a
	rm -f lib/libmetag.a
