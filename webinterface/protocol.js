// Browser limits bound memory use; the native CLI remains unrestricted.
export const LIMITS = Object.freeze({ minLength: 7, maxLength: 2000, maxPairs: 1000000 });
export const CSV_COLUMNS = 'id1,start1,end1,id2,start2,end2,E';
export const DEFAULTS = Object.freeze({
  query: '>query\nGGGAGGAGGG',
  target: '>target\nCCCUCCUCCC',
  count: 10,
});

export function parseSequence(text, label) {
  if (typeof text !== 'string' || text.length > 10000) {
    throw new Error(`${label}: enter one sequence of at most ${LIMITS.maxLength} nucleotides.`);
  }
  const lines = text.trim().split(/\r?\n/).map((line) => line.trim());
  let id = label.toLowerCase();
  if (lines[0].startsWith('>')) {
    id = lines.shift().slice(1).trim();
    if (!id || id.length > 200 || /[\x00-\x1f\x7f]/.test(id)) {
      throw new Error(`${label}: the FASTA header needs a name (at most 200 characters).`);
    }
  }
  if (lines.some((line) => line.includes('>'))) {
    throw new Error(`${label}: enter exactly one FASTA record, with its header on the first line.`);
  }
  const sequence = lines.join('').replace(/\s/g, '').toUpperCase();
  if (!sequence || /[^ACGUTN]/.test(sequence)) {
    throw new Error(`${label}: use only A, C, G, U, T or N in the sequence.`);
  }
  if (sequence.length < LIMITS.minLength || sequence.length > LIMITS.maxLength) {
    throw new Error(`${label}: use ${LIMITS.minLength}–${LIMITS.maxLength} nucleotides (the default seed needs 7).`);
  }
  const rna = sequence.replace(/T/g, 'U');
  return { id, sequence: rna, fasta: `>${id}\n${rna}\n` };
}

export function prepareRequest(request) {
  const query = parseSequence(request?.query, 'Query');
  const target = parseSequence(request?.target, 'Target');
  if (query.sequence.length * target.sequence.length > LIMITS.maxPairs) {
    throw new Error('For this browser version, query length × target length must be at most 1,000,000. Use shorter sequences or the command line.');
  }
  const count = Number(request.count);
  if (!Number.isInteger(count) || count < 1 || count > 100) {
    throw new Error('Choose between 1 and 100 interactions.');
  }
  return { query, target, count };
}

export function makeArguments(count) {
  return [
    '--query=/query.fa', '--target=/target.fa', '--outMode=C',
    `--outNumber=${count}`, '--outOverlap=B',
    `--outCsvCols=${CSV_COLUMNS}`, '--out=/results.csv',
  ];
}
