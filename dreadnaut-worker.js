import MyModule from './dreadnaut.js';

let isNode = typeof process !== 'undefined' && process.versions && process.versions.node;

async function initWorkerEnv() {
	if (isNode) {
		const wt = await import('worker_threads');
		globalThis.postWorkerMessage = (msg) => wt.parentPort.postMessage(msg);
		wt.parentPort.on('message', (data) => {
			handleCommand({ data });
		});
	} else {
		globalThis.postWorkerMessage = (msg) => postMessage(msg);
		globalThis.onmessage = (e) => handleCommand(e);
	}
}

function handleCommand(e) {
	const { type, data } = e.data;
	if (type === 'command') {
		if (!globalThis.inputBuffer) globalThis.inputBuffer = [];
		for (let i = 0; i < data.length; i++) {
			globalThis.inputBuffer.push(data.charCodeAt(i));
		}
		if (!data.endsWith('\n')) {
			globalThis.inputBuffer.push(10); // newline
		}
		if (globalThis.resolveInput) {
			let wakeUp = globalThis.resolveInput;
			globalThis.resolveInput = null;
			wakeUp(globalThis.inputBuffer.shift());
		}
	}
}

const setupWorker = async () => {
	await initWorkerEnv();
	await MyModule({
		print: (text) => {
			if (!globalThis.accumulatedOut) globalThis.accumulatedOut = "";
			globalThis.accumulatedOut += text + '\n';
		},
		printErr: (text) => {
			if (!globalThis.accumulatedOut) globalThis.accumulatedOut = "";
			globalThis.accumulatedOut += text + '\n';
		},
	});
};

setupWorker().then(() => {
	globalThis.postWorkerMessage({ type: 'ready' });
}).catch(err => {
	if (globalThis.postWorkerMessage) {
		globalThis.postWorkerMessage({ type: 'error', data: err.message });
	}
});
