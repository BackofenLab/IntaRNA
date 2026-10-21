// Record the versions of the compiled executable, alongside its exact assets.
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
const [repo, site, emscripten, commit, changes] = process.argv.slice(2);
const { default: createIntaRNA } = await import(pathToFileURL(resolve(site, 'wasm/intarna.js')));
const lines = [];
const module = await createIntaRNA({ noInitialRun: true, locateFile: () => pathToFileURL(resolve(site, 'wasm/intarna.wasm')).href, print: (line) => lines.push(line), printErr: (line) => lines.push(line) });
if (module.callMain(['--version']) !== 0) throw new Error('Could not read compiled engine versions');
const version = lines.join('\n');
const match = /IntaRNA\s+(\S+)[\s\S]*Vienna RNA package\s+(\S+)\s+and boost\s+(\S+)/i.exec(version);
if (!match) throw new Error(`Unrecognized engine version output: ${version}`);
const metadata = {
  source: 'https://github.com/BackofenLab/IntaRNA', commit,
  dirty: !!changes,
  intarna: match[1], viennarna: match[2], boost: match[3], emscripten, zlib: '1.3.1',
  threading: 'single worker per task',
  sha256: Object.fromEntries(['intarna.js', 'intarna.wasm'].map((name) => [name, createHash('sha256').update(readFileSync(resolve(site, 'wasm', name))).digest('hex')])),
};
writeFileSync(resolve(site, 'build-info.json'), JSON.stringify(metadata, null, 2) + '\n');
writeFileSync(resolve(site, 'build-info.txt'), `${version}\nIntaRNA repository: ${metadata.source}\nIntaRNA revision: ${metadata.commit}${metadata.dirty ? ' (modified sources)' : ''}\nEmscripten: ${emscripten}\nzlib: ${metadata.zlib}\n${Object.entries(metadata.sha256).map(([name, hash]) => `SHA-256 ${name}: ${hash}`).join('\n')}\n`);
console.log(`Recorded compiled versions: IntaRNA ${match[1]}, ViennaRNA ${match[2]}, Boost ${match[3]}`);
