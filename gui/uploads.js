import { parseFasta, LIMITS } from './protocol.js';

/** Read and validate a local text file before replacing a FASTA input. */
export async function readFastaFile(file, label, task) {
  if (!file || file.size > LIMITS.maxFastaBytes) throw new Error(`${label}: choose a FASTA text file of at most 1 MB.`);
  const text = await file.text();
  parseFasta(text, label, task !== 'evaluation');
  return text;
}
