// Column IDs and descriptions from OutputHandlerCsv.h; seq1 = target, seq2 = query.
export const CSV_CATALOG = [
  {
    "key": "id1",
    "label": "Target",
    "description": "id of first sequence",
    "numeric": false
  },
  {
    "key": "id2",
    "label": "Query",
    "description": "id of second sequence",
    "numeric": false
  },
  {
    "key": "seq1",
    "label": "seq1",
    "description": "full first sequence",
    "numeric": false
  },
  {
    "key": "seq2",
    "label": "seq2",
    "description": "full second sequence",
    "numeric": false
  },
  {
    "key": "subseq1",
    "label": "subseq1",
    "description": "interacting subsequence of first sequence",
    "numeric": false
  },
  {
    "key": "subseq2",
    "label": "subseq2",
    "description": "interacting subsequence of second sequence",
    "numeric": false
  },
  {
    "key": "subseqDP",
    "label": "subseqDP",
    "description": "hybrid subsequences compatible with hybridDP",
    "numeric": false
  },
  {
    "key": "subseqDB",
    "label": "subseqDB",
    "description": "hybrid subsequences compatible with hybridDB",
    "numeric": false
  },
  {
    "key": "start1",
    "label": "Target start",
    "description": "start index of hybrid in seq1",
    "numeric": true
  },
  {
    "key": "end1",
    "label": "Target end",
    "description": "end index of hybrid in seq1",
    "numeric": true
  },
  {
    "key": "start2",
    "label": "Query start",
    "description": "start index of hybrid in seq2",
    "numeric": true
  },
  {
    "key": "end2",
    "label": "Query end",
    "description": "end index of hybrid in seq2",
    "numeric": true
  },
  {
    "key": "hybridDP",
    "label": "hybridDP",
    "description": "hybrid in VRNA dot-bracket notation",
    "numeric": false
  },
  {
    "key": "hybridDB",
    "label": "hybridDB",
    "description": "hybrid in dot-bar notation",
    "numeric": false
  },
  {
    "key": "hybridDPfull",
    "label": "hybridDPfull",
    "description": "hybrid in VRNA dot-bracket notation for full sequence lengths",
    "numeric": false
  },
  {
    "key": "hybridDBfull",
    "label": "hybridDBfull",
    "description": "hybrid in dot-bar notation for full sequence lengths",
    "numeric": false
  },
  {
    "key": "bpList",
    "label": "bpList",
    "description": "list of hybrid base pairs like (1,3),(4,2),...",
    "numeric": false
  },
  {
    "key": "E",
    "label": "Energy (kcal/mol)",
    "description": "overall interaction energy",
    "numeric": true
  },
  {
    "key": "Etotal",
    "label": "Etotal",
    "description": "overall energy of the interaction including the intra-molecular ensemble energies (E+Eall1+Eall2)",
    "numeric": true
  },
  {
    "key": "ED1",
    "label": "ED1",
    "description": "ED value of seq1",
    "numeric": true
  },
  {
    "key": "ED2",
    "label": "ED2",
    "description": "ED value of seq2",
    "numeric": true
  },
  {
    "key": "Pu1",
    "label": "Pu1",
    "description": "probability to be accessible for seq1",
    "numeric": true
  },
  {
    "key": "Pu2",
    "label": "Pu2",
    "description": "probability to be accessible for seq2",
    "numeric": true
  },
  {
    "key": "E_init",
    "label": "E_init",
    "description": "initiation energy",
    "numeric": true
  },
  {
    "key": "E_loops",
    "label": "E_loops",
    "description": "sum of loop energies (excluding E_init)",
    "numeric": true
  },
  {
    "key": "E_dangleL",
    "label": "E_dangleL",
    "description": "dangling end contribution of base pair (start1,end2)",
    "numeric": true
  },
  {
    "key": "E_dangleR",
    "label": "E_dangleR",
    "description": "dangling end contribution of base pair (end1,start2)",
    "numeric": true
  },
  {
    "key": "E_endL",
    "label": "E_endL",
    "description": "penalty of closing base pair (start1,end2)",
    "numeric": true
  },
  {
    "key": "E_endR",
    "label": "E_endR",
    "description": "penalty of closing base pair (end1,start2)",
    "numeric": true
  },
  {
    "key": "E_hybrid",
    "label": "E_hybrid",
    "description": "energy of hybridization only = E - ED1 - ED2",
    "numeric": true
  },
  {
    "key": "E_norm",
    "label": "E_norm",
    "description": "length normalized energy = E/ln(length(seq1)*length(seq2))",
    "numeric": true
  },
  {
    "key": "E_hybridNorm",
    "label": "E_hybridNorm",
    "description": "length normalized energy of hybridization only = E_hybrid / ln(length(seq1)*length(seq2))",
    "numeric": true
  },
  {
    "key": "E_add",
    "label": "E_add",
    "description": "user provided energy shift",
    "numeric": true
  },
  {
    "key": "seedStart1",
    "label": "seedStart1",
    "description": "start index of the seed in seq1",
    "numeric": true
  },
  {
    "key": "seedEnd1",
    "label": "seedEnd1",
    "description": "end index of the seed in seq1",
    "numeric": true
  },
  {
    "key": "seedStart2",
    "label": "seedStart2",
    "description": "start index of the seed in seq2",
    "numeric": true
  },
  {
    "key": "seedEnd2",
    "label": "seedEnd2",
    "description": "end index of the seed in seq2",
    "numeric": true
  },
  {
    "key": "seedE",
    "label": "seedE",
    "description": "overall energy of the seed only (including accessibility etc)",
    "numeric": true
  },
  {
    "key": "seedED1",
    "label": "seedED1",
    "description": "ED value of seq1 of the seed only (excluding rest)",
    "numeric": true
  },
  {
    "key": "seedED2",
    "label": "seedED2",
    "description": "ED value of seq2 of the seed only (excluding rest)",
    "numeric": true
  },
  {
    "key": "seedPu1",
    "label": "seedPu1",
    "description": "probability of seed region to be accessible for seq1",
    "numeric": true
  },
  {
    "key": "seedPu2",
    "label": "seedPu2",
    "description": "probability of seed region to be accessible for seq2",
    "numeric": true
  },
  {
    "key": "w",
    "label": "w",
    "description": "Boltzmann weight of the interaction energy",
    "numeric": true
  },
  {
    "key": "Eall",
    "label": "Eall",
    "description": "ensemble energy of all interactions (outConstraint.needZall)",
    "numeric": true
  },
  {
    "key": "Eall1",
    "label": "Eall1",
    "description": "ensemble energy of all intra-molecular structures of seq1",
    "numeric": true
  },
  {
    "key": "Eall2",
    "label": "Eall2",
    "description": "ensemble energy of all intra-molecular structures of seq2",
    "numeric": true
  },
  {
    "key": "Zall",
    "label": "Zall",
    "description": "partition function of all interactions (outConstraint.needZall)",
    "numeric": true
  },
  {
    "key": "Zall1",
    "label": "Zall1",
    "description": "partition function of all intra-molecular structures of seq1",
    "numeric": true
  },
  {
    "key": "Zall2",
    "label": "Zall2",
    "description": "partition function of all intra-molecular structures of seq2",
    "numeric": true
  },
  {
    "key": "EallTotal",
    "label": "EallTotal",
    "description": "total ensemble energy (Eall+Eall1+Eall2) of all interactions including the intra-molecular ensemble energies (outConstraint.needZall)",
    "numeric": true
  },
  {
    "key": "P_E",
    "label": "P_E",
    "description": "probability of mfe within interaction ensemble (outConstraint.needZall)",
    "numeric": true
  },
  {
    "key": "RT",
    "label": "RT",
    "description": "the scaled temperature used for Boltzmann weight computation",
    "numeric": true
  }
];
