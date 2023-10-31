C = gcc
CFLAGS = -lstdc++

all: kmc_query kmc_api

KMC_API_DIR = kmc_api
KMC_API_OBJS = \
	mmer.o \
	kmc_file.o \
	kmer_api.o
KMC_API_SRC_FILES = $(KMC_API_DIR)/*.cpp

kmc_api: $(KMC_API_SRC_FILES) 
	$(C) $(CFLAGS) -o $(KMC_API_OBJS) 

kmc_query: kmc_query.cpp 
	$(C) $(CFLAGS) -o kmc_query.out kmc_query.cpp kmc_api/*.cpp

clean:
	-rm -f kmc_query.out
	-rm -f $(KMC_API_DIR)/*.o
