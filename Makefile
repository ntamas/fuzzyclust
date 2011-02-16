IGRAPH_CFLAGS=$(shell pkg-config igraph --cflags)
IGRAPH_LDFLAGS=$(shell pkg-config igraph --libs)
STATIC_IGRAPH_LDFLAGS=$(shell pkg-config igraph --libs --static)

CFLAGS=$(IGRAPH_CFLAGS) -g -Wall -O3
LDFLAGS=$(IGRAPH_LDFLAGS)
STATIC_LDFLAGS=$(STATIC_IGRAPH_LDFLAGS)

BINARY=fuzzyclust
STATIC_BINARY=$(BINARY)-static
SRC=src/fuzzyclust.c
TARGETS=build/fuzzyclust.o
DIST=$(SRC) doc/COPYING doc/fuzzyclust.1 data/ Makefile

.PHONY: prepare clean release

all: prepare $(BINARY)
static: prepare $(STATIC_BINARY)

prepare:
	mkdir -p build

$(BINARY): $(TARGETS)
	$(CC) -o $@ $(LDFLAGS) $<
$(STATIC_BINARY): $(TARGETS)
	$(CC) -static -static-libgcc -o $@ $< $(STATIC_LDFLAGS)

dist:
	tar -cvvzf $(BINARY).tar.gz $(DIST)
dist-static: $(STATIC_BINARY)
	strip $(STATIC_BINARY)
	tar -cvvzf $(BINARY)-static.tar.gz $(DIST) $(STATIC_BINARY)

clean:
	-rm -rf $(BINARY) $(STATIC_BINARY) $(TARGETS)

build/%.o: src/%.c
	$(CC) -c -o $@ $< $(CFLAGS)

