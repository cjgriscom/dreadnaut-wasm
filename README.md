# Nauty and Traces v2.9.3

An unofficial mirror of nauty and Traces v2.9.3 - see README

# dreadnaut-wasm

A WebAssembly port of the `dreadnaut` program (from the [nauty and Traces](https://pallini.di.uniroma1.it/) graph isomorphism suite) designed to run as a Web Worker in the browser or via `worker_threads` in Node.js.

## Overview

This project modifies the original `dreadnaut.c` to compile with Emscripten, utilizing `ASYNCIFY` to perfectly replicate the interactive CLI experience asynchronously. Rather than waiting for the entire program to execute, the Web Worker mimics a true terminal:
- Automatically yields execution when it awaits input.
- Accepts standard dreadnaut commands (like `n 5 g 0:1 1:2 . x`) seamlessly via message passing.
- Flushes standard output (`stdout`) directly back to the parent script.

## Prerequisites

To build this project, you need the [Emscripten SDK (emsdk)](https://emscripten.org/docs/getting_started/downloads.html) installed and activated in your environment, providing the `emcc` compiler.

## Building

To build the project into JavaScript and WebAssembly (`dreadnaut.js` and `dreadnaut.wasm`), simply run the included build script:

```bash
./build-emscripten.sh
```

## Testing

A Node.js test script `test-emscripten.js` is included to verify the build successfully processes commands asynchronously and returns the correct calculations. Run it using:

```bash
node test-emscripten.js
```

## Usage

Integrating `dreadnaut-wasm` into your JavaScript application is straightforward through the `dreadnaut-worker.js` wrapper.

### 1. Instantiate the Worker

**In a Browser Web App:**
```javascript
const worker = new Worker('dreadnaut-worker.js');
```

**In Node.js:**
```javascript
const { Worker } = require('worker_threads');
const worker = new Worker('./dreadnaut-worker.js');
```

### 2. Communicate with Dreadnaut

Send standard `dreadnaut` interactive commands to the worker just as you would type them in the CLI. For instance, to calculate the automorphism group of a 5-cycle graph:
```javascript
// Wait for Dreadnaut to initialize
worker.onmessage = (e) => {
    if (e.data.type === 'ready') {
        // Send our command string
        worker.postMessage({ type: 'command', data: 'n 5 g 0:1 1:2 2:3 3:4 4:0 . x' });
    }
    else if (e.data.type === 'output') {
        console.log("Dreadnaut output:\n" + e.data.data);
    }
    else if (e.data.type === 'error') {
        console.error("Dreadnaut error:\n" + e.data.data);
    }
};
```

## Technical Details

- **Emscripten Asyncify:** `dreadnaut.wasm` is compiled with Emscripten's `ASYNCIFY` flag, mapping the synchronous character readings in C via `getc()` down to yielded JavaScript `Promises`.
- **Node vs Browser Setup:** The `dreadnaut-worker.js` gracefully determines whether it is executing inside a web browser or Node.js environment, falling back to ES6/CommonJS worker patterns seamlessly.

## License

This project adheres to the same licensing terms as the original `nauty and Traces` package. Check the local `COPYRIGHT` or `LICENSE-2.0.txt` files provided by the original authors.
