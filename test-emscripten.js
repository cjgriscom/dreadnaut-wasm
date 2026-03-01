import { Worker } from 'worker_threads';

async function testWorker() {
	console.log("Starting worker...");
	const worker = new Worker('./dreadnaut-worker.js');

	worker.on('message', (msg) => {
		if (msg.type === 'ready') {
			console.log("Worker ready. Sending command...");
			worker.postMessage({ type: 'command', data: 'n 5 g 0:1 1:2 2:3 3:4 4:0 . x' });
		} else if (msg.type === 'output') {
			console.log("OUTPUT RECEIVED:\n" + msg.data);
			if (msg.data.includes("1 orbit; grpsize=10; 2 gens; 6 nodes; maxlev=3")) {
				console.log("Test passed!");
				process.exit(0);
			}
		} else if (msg.type === 'error') {
			console.error("Worker error message:", msg.data);
			process.exit(1);
		}
	});

	worker.on('error', (err) => {
		console.error("Worker error:", err);
		process.exit(1);
	});

	worker.on('exit', (code) => {
		console.log(`Worker stopped with exit code ${code}`);
	});
}

testWorker();
