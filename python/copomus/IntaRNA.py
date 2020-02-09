# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

import sys
from typing import Dict, Union, List
from subprocess import Popen, PIPE, run


class IntaRNA:
    """IntaRNA python3 wrapper"""
    def __init__(self):
        pass

    @staticmethod
    def check_binary(required_version: str):
        """
        Checks weather or not the IntaRNA binary can be found in PATH and exits with returncode 1 if not found.
        :param required_version: The minimum required version for IntaRNA
        """
        p = run('IntaRNA --version', shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if p.returncode:
            sys.stderr.write('Please add the IntaRNA binary to your PATH.\n')
            sys.exit(1)
        if p.stdout.split('\n')[0].split()[1] < required_version:
            sys.stderr.write(f'This tool requires IntaRNA version {required_version}, please upgrade\n')
            sys.exit(1)

    def run(self, query: str, target: str, qidxpos0: int, tidxpos0: int, outcsvcols: str, threads: int, n: int = 1,
            param_file: str = '', extra_params: list = '', raw_stdout: bool = False) -> \
            Union[Dict[str, any], List[Dict[str, any]], str]:
        """
        Runs IntaRNA on a pair of sequences and returns a dictionary with the output
        :param query: The query sequence
        :param target: The target sequence
        :param qidxpos0: The starting index for the query sequence
        :param tidxpos0: The starting index for the target sequence
        :param outcsvcols: The csvcols returned by IntaRNA
        :param threads: Thread count used by IntaRNA
        :param n: The number of suboptimal interactions returned
        :param param_file: Optional path to a file with additional parameters
        :param extra_params: Optional list of additional parameters
        :param raw_stdout: Weather or not stdout should be returned in raw format, not as a dict
        :return: Dictionary with column name as key and payload as value
        """
        if not extra_params:
            extra_params = []
        param_file = f'--parameterFile={param_file}' if param_file else ''
        p = Popen(['IntaRNA', '-q', query, '-t', target, '--outMode=C', f'--outcsvcols={outcsvcols}',
                   f'--qIdxPos0={qidxpos0}', f'--tIdxPos0={tidxpos0}', f'--outNumber={n}', f'--threads={threads}', param_file]
                  + extra_params, stdout=PIPE, stderr=PIPE, universal_newlines=True)

        stdout, stderr = p.communicate()
        if p.returncode:
            sys.stderr.write(f'IntaRNA Error: {stdout} {stderr}\n')
            sys.exit(p.returncode)
        if not raw_stdout:
            if n == 1:
                d = self.output_to_dict(stdout)
                if type(d) != bool:
                    return d[0]
                else:
                    return d
            return self.output_to_dict(stdout)
        else:
            return stdout

    @staticmethod
    def output_to_dict(data: str) -> Union[List[Dict[str, any]], bool]:
        """
        Transforms IntaRNA output into a dictionary
        :param data: Raw IntaRNA output
        :return: dictionary with colname as key and payload as value (casted into fitting type)
        >>> d = IntaRNA.output_to_dict('start1;end1;id1;id2;E\\n34;42;target;query;-4.7\\n')
        >>> d
        [{'start1': 34, 'end1': 42, 'id1': 'target', 'id2': 'query', 'E': -4.7}]
        >>> [type(x) for x in d[0].values()]
        [<class 'int'>, <class 'int'>, <class 'str'>, <class 'str'>, <class 'float'>]
        """
        subopts = []

        lines = data.rstrip().split('\n')
        if len(lines) == 2:
            headers, entries = lines
            entries = [entries]
        elif len(lines) > 2:
            headers = lines[0]
            entries = lines[1:]
        else:
            headers = lines[0]
            entries = ['']
        headers = headers.split(';')

        for values in entries:
            d = {}
            values = values.split(';')
            if len(headers) != len(values):  # no interaction between query and target
                return False

            for i, h in enumerate(headers):
                try:
                    value = int(values[i])  # try to cast to int
                except ValueError:
                    try:
                        value = float(values[i])  # try to cast to float
                    except ValueError:
                        value = values[i]  # leave as string
                d[h] = value
            subopts.append(d)

        for subopt in subopts:
            for k, v in subopt.items():
                if not v:
                    subopt[k] = None

        return subopts

    @staticmethod
    def to_fasta(sequences: Dict[str, str]) -> str:
        """
        Transforms sequences into a FASTA string
        :param sequences: A dictionary with sequence name as key and raw sequence as tuple
        :return: A FASTA string containing all sequences
        """
        fasta_str = ''
        for name, seq in sequences.items():
            fasta_str += f'>{name}\n{seq}\n'
        return fasta_str
