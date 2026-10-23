import { makeArguments } from './protocol.js';
import { projectCsv } from './results.js';

const encoder = new TextEncoder();
const crcTable = Uint32Array.from({ length: 256 }, (_, value) => {
  for (let bit = 0; bit < 8; bit++) value = (value >>> 1) ^ (value & 1 ? 0xedb88320 : 0);
  return value >>> 0;
});
function crc32(bytes) {
  let crc = 0xffffffff;
  for (const byte of bytes) crc = crcTable[(crc ^ byte) & 255] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

// A standard, uncompressed ZIP: local generation needs no CDN or dependency.
// Fixed timestamps and UTF-8 names make the same completed run reproducible.
export function zipFiles(files) {
  const local = [], central = [], names = new Set();
  let offset = 0, directorySize = 0;
  if (files.length > 65535) throw new Error('Too many files for a ZIP archive.');
  for (const file of files) {
    if (!file.name || /[/\\\x00]/.test(file.name) || ['.', '..'].includes(file.name) || names.has(file.name)) throw new Error('Invalid or duplicate archive filename.');
    names.add(file.name);
    const name = encoder.encode(file.name), data = encoder.encode(file.text), crc = crc32(data);
    if (name.length > 65535 || offset + name.length + data.length + 30 > 0xffffffff) throw new Error('ZIP archive is too large.');
    const header = new Uint8Array(30), view = new DataView(header.buffer);
    view.setUint32(0, 0x04034b50, true); view.setUint16(4, 20, true);
    view.setUint16(6, 0x800, true); view.setUint16(12, 33, true); // 1980-01-01
    view.setUint32(14, crc, true); view.setUint32(18, data.length, true); view.setUint32(22, data.length, true);
    view.setUint16(26, name.length, true);
    local.push(header, name, data);
    const entry = new Uint8Array(46), record = new DataView(entry.buffer);
    record.setUint32(0, 0x02014b50, true); record.setUint16(4, 20, true); record.setUint16(6, 20, true);
    record.setUint16(8, 0x800, true); record.setUint16(14, 33, true);
    record.setUint32(16, crc, true); record.setUint32(20, data.length, true); record.setUint32(24, data.length, true);
    record.setUint16(28, name.length, true); record.setUint32(42, offset, true);
    central.push(entry, name); directorySize += entry.length + name.length;
    offset += header.length + name.length + data.length;
  }
  const end = new Uint8Array(22), view = new DataView(end.buffer);
  view.setUint32(0, 0x06054b50, true);
  view.setUint16(8, files.length, true); view.setUint16(10, files.length, true);
  view.setUint32(12, directorySize, true); view.setUint32(16, offset, true);
  const zip = new Uint8Array(offset + directorySize + end.length);
  let position = 0;
  for (const bytes of [...local, ...central, end]) { zip.set(bytes, position); position += bytes.length; }
  return zip;
}

export function archiveFiles(request, result, provenance) {
  if (!provenance?.trim()) throw new Error('Compiled version information is unavailable');
  const versionComments = provenance.trimEnd().split(/\r?\n/).map((line) => `# ${line}`).join('\n');
  const reports = [['normal', 'N'], ['detailed', 'D'], ['ensemble', 'E']].filter(([key]) => request.outputs.includes(key));
  const config = (format, filename) => {
    // Reuse the worker's effective arguments, replacing only MEMFS paths and
    // internal visualization columns. Optional data are shipped as text files.
    const args = makeArguments(request, format).filter((arg) => !['--qIdxPos0=1', '--tIdxPos0=1'].includes(arg));
    const settings = args.map((arg) => {
      if (arg.startsWith('--outCsvCols=')) return `outCsvCols = ${request.csvColumns.join(',')}`;
      return arg.slice(2).replace('=', ' = ').replace(/\/(query\.fa|target\.fa|[qt](?:Shape|AccFile)\.txt)/g, '$1')
        .replace('/results.csv', 'intarna-interactions.csv').replace('/report.txt', `intarna-${filename}.txt`).replace('/output/', '');
    });
    return `${versionComments}\n# Run from the extracted archive directory with these compiled library versions.\n# Reproduce: IntaRNA --parameterFile=${filename}.cfg\n# Only non-default analysis settings are listed; input/output routing is explicit.\n${settings.join('\n')}\n`;
  };
  const parameters = config('C', 'parameters') + reports.map(([key]) => `# Additional report: IntaRNA --parameterFile=${key}.cfg\n`).join('');
  const files = [
    { name: 'query.fa', text: request.query.fasta }, { name: 'target.fa', text: request.target.fasta },
    { name: 'parameters.cfg', text: parameters }, { name: 'build-info.txt', text: provenance },
    { name: 'intarna-interactions.csv', text: projectCsv(result.csv, request.csvColumns) },
    ...result.files,
    ...reports.map(([key, format]) => ({ name: `${key}.cfg`, text: config(format, key) })),
  ];
  for (const [key, value] of Object.entries(request.parameters)) {
    if (/^(q|t)(Shape|AccFile)$/.test(key) && value) files.push({ name: `${key}.txt`, text: value + '\n' });
  }
  return files;
}

export function createArchive(request, result, provenance) {
  return zipFiles(archiveFiles(request, result, provenance));
}
