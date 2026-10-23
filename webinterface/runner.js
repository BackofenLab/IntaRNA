// The main thread stays available even while synchronous C++ is calculating.
export class PredictionRunner {
  constructor(workerFactory = () => new Worker(new URL('./intarna-worker.js', import.meta.url), { type: 'module' })) {
    this.workerFactory = workerFactory;
    this.worker = null;
    this.cancel = null;
  }

  run(request, onStatus = () => {}) {
    if (this.worker) return Promise.reject(new Error('A prediction is already running.'));
    return new Promise((resolve, reject) => {
      const worker = this.workerFactory();
      this.worker = worker;
      const finish = (error, result) => {
        if (this.worker !== worker) return; // Ignore queued messages after abort.
        worker.terminate();
        this.worker = null;
        this.cancel = null;
        if (error) reject(error);
        else resolve(result);
      };
      this.cancel = () => finish(new DOMException('Computation aborted.', 'AbortError'));
      worker.onmessage = ({ data }) => {
        if (this.worker !== worker) return;
        if (data.type === 'status') onStatus(data.message);
        else if (data.type === 'result') finish(null, data);
        else if (data.type === 'error') {
          const error = new Error(data.message);
          error.diagnostics = data.diagnostics;
          finish(error);
        }
      };
      worker.onerror = (event) => {
        event.preventDefault();
        finish(new Error('Could not load or run the prediction engine. Reload the page and check that intarna.js and intarna.wasm are served together.'));
      };
      worker.onmessageerror = () => finish(new Error('Could not read the prediction engine response.'));
      try {
        worker.postMessage(request);
      } catch (error) {
        finish(error);
      }
    });
  }

  abort() {
    this.cancel?.();
  }
}
