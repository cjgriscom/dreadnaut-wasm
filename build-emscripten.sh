#!/bin/bash
set -e

# Emscripten build script for dreadnaut to produce a Web Worker

emconfigure ./configure
emmake make dreadnaut

# 3. Final linking for Web Worker
# -s ASYNCIFY: allows synchronous-looking I/O (like reading from stdin)
# -s MODULARIZE: wraps the output in a function for better control
# -s EXPORT_ES6: uses ES6 modules
# -s EXPORTED_RUNTIME_METHODS: allows us to call FS methods to handle I/O
# -s ALLOW_MEMORY_GROWTH: dreadnaut can use a lot of memory for large graphs

echo "Finalizing dreadnaut.js for Web Worker..."

emcc -O3 \
    -s ASYNCIFY=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORTED_RUNTIME_METHODS="['FS', 'callMain', 'TTY']" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s ENVIRONMENT=node,worker \
    -o dreadnaut.js \
    dreadnaut.c \
    naututil.o \
    nautinv.o \
    gtools.o \
    traces.o \
    nauty.o \
    nautil.o \
    nausparse.o \
    naugraph.o \
    schreier.o \
    naurng.o

echo "Build complete: dreadnaut.js and dreadnaut.wasm"
