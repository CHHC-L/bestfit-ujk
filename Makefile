# Best Fit under U{j,k}: build, verify, test.
#   make            build the tools
#   make verify     verify every certificate (two exact C++ verifiers) and regenerate the result tables
#   make verify-py  the same, plus the Python verifier on certificates with <= 300k patterns (slow)
#   make test       quick self-tests (chain moves vs explicit Best Fit, a few certificates end to end)
#   make simcheck   simulate every certified pair and compare with the certificate (sanity check)
#   make figures    figures/classification.png, figures/boundary.png (needs matplotlib)
CXXFLAGS ?= -O3 -std=c++17 -Wall
CFLAGS ?= -O3 -Wall
PYTHON ?= python3
JOBS ?= 8

BINS = bfcert bfvi verify2 bfsim bfstats

all: $(BINS)

bfcert: src/bfcert.cpp src/bfchain.hpp
	$(CXX) $(CXXFLAGS) -o $@ src/bfcert.cpp
bfvi: src/bfvi.cpp src/bfchain.hpp
	$(CXX) $(CXXFLAGS) -pthread -o $@ src/bfvi.cpp
verify2: src/verify2.cpp
	$(CXX) $(CXXFLAGS) -o $@ src/verify2.cpp
bfsim: src/bfsim.c
	$(CC) $(CFLAGS) -o $@ src/bfsim.c
bfstats: src/bfstats.c
	$(CC) $(CFLAGS) -o $@ src/bfstats.c -lm

verify: bfcert verify2
	$(PYTHON) scripts/verify_all.py --jobs $(JOBS)
	$(PYTHON) scripts/make_tables.py

verify-py: bfcert verify2
	$(PYTHON) scripts/verify_all.py --jobs $(JOBS) --python 300000
	$(PYTHON) scripts/make_tables.py

test: bfcert verify2
	./bfcert test 8 11 2000000
	./bfcert test 12 24 2000000
	./bfcert test 5 25 2000000 smalld0
	./bfcert test 21 25 2000000
	$(PYTHON) scripts/verify_all.py --only 8_11_ --python 300000
	$(PYTHON) scripts/e2e_check.py certs/8_11_linear_vi_M1-3-7-9-10_R.cert.gz \
	    certs/9_13_linear_vi_M1-2-1-2-2-2_R.cert.gz certs/8_12_bounded_vi_M12-2-1-2-2_R.cert.gz \
	    certs/9_13_bounded2_vi_M2-2-1-2-2-2_R.cert.gz certs/9_13_ylo_vi_M0-6-6-6-6-6_R.cert.gz \
	    certs/9_13_yhi_vi_M0-6-6-6-6-6_R.cert.gz

simcheck: bfsim
	$(PYTHON) scripts/sim_check.py --jobs $(JOBS)

figures:
	$(PYTHON) scripts/make_figures.py

clean:
	rm -f $(BINS)
	rm -rf lp

.PHONY: all verify verify-py test simcheck figures clean
