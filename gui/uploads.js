import { parseFasta } from './protocol.js';

/** Read and validate a local text file before replacing a FASTA input. */
export async function readFastaFile(file, label, task) {
  if (!file || file.size > 60000) throw new Error(`${label}: choose a FASTA text file of at most 60 KB.`);
  const text = await file.text();
  parseFasta(text, label, task !== 'evaluation');
  return text;
}
