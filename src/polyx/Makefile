all: polyx polyx_threads

polyx:
	nim c --threads:on -p:release --opt:speed  -o:polyx polyX_example.nim

polyx_threads:
	nim c --threads:on -p:release --opt:speed -o:polyx_threads polyX_threads.nim

test.fq.gz:
	python3 makePolyX.py -n 100000 -o test.fq.gz --max-polytail 1000 --max-length 1000

clean:
	if [ -f polyx ]; then rm polyx; fi
	if [ -f polyx_threads ]; then rm polyx_threads; fi
	if [ -f test.fq.gz ]; then rm test.fq.gz; fi
	if [ -f test_trimmed.fq ]; then rm test_trimmed.fq; fi
	if [ -f test_trimmed_threads.fq ]; then rm test_trimmed_threads.fq; fi

test: polyx polyx_threads test.fq.gz
	./polyx test.fq.gz > test_trimmed.fq
	./polyx_threads test.fq.gz > test_trimmed_threads.fq
	# if hyperfine is available, use it to benchmark
	# first check if hyperfine is installed
	if command -v hyperfine &> /dev/null; then \
		hyperfine --warmup 1 --max-runs 8 --export-markdown benchmark.md './polyx test.fq.gz' './polyx_threads test.fq.gz'; \
	else \
		echo "hyperfine not found, skipping benchmark"; \
	fi
	md5sum test*fq