import { test } from 'node:test';
import assert from 'node:assert/strict';
import { PredictionRunner } from '../runner.js';
import { DEFAULTS } from '../protocol.js';

function fakeWorker() {
  return { terminated: false, terminate() { this.terminated = true; }, postMessage() {}, send(data) { this.onmessage({ data }); } };
}

test('abort terminates computation, rejects and permits a clean restart', async () => {
  const workers = [];
  const runner = new PredictionRunner(() => { const worker = fakeWorker(); workers.push(worker); return worker; });
  const first = runner.run(DEFAULTS);
  await assert.rejects(runner.run(DEFAULTS), /already running/);
  runner.abort();
  await assert.rejects(first, { name: 'AbortError' });
  assert.equal(workers[0].terminated, true);
  const second = runner.run(DEFAULTS);
  workers[0].send({ type: 'result', csv: 'stale' });
  workers[1].send({ type: 'result', csv: 'fresh' });
  assert.equal((await second).csv, 'fresh');
  assert.equal(workers[1].terminated, true);
});

test('runtime, load and message errors release the worker', async () => {
  for (const kind of ['error', 'load', 'message', 'post']) {
    const worker = fakeWorker();
    if (kind === 'post') worker.postMessage = () => { throw new Error('post failed'); };
    const runner = new PredictionRunner(() => worker);
    const result = runner.run(DEFAULTS);
    if (kind === 'error') worker.send({ type: 'error', message: 'allocation failed', diagnostics: 'details' });
    if (kind === 'load') worker.onerror({ preventDefault() {} });
    if (kind === 'message') worker.onmessageerror();
    await assert.rejects(result);
    assert.equal(worker.terminated, true);
    assert.equal(runner.worker, null);
  }
});
