.PHONY: all clean verify

PRECISION ?= dd

all:
	@mkdir -p cpp/build && cd cpp/build && cmake -DBUILD_SHARED_LIBS=ON .. && make -j

clean:
	rm -rf cpp/build

verify: all
	./cpp/build/verify qaoa_data/ --precision $(PRECISION)
