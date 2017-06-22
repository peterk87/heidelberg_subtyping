import os
from subprocess import Popen, PIPE
import logging

import re


def run_command(cmdlist):
    p = Popen(cmdlist,
              stdout=PIPE,
              stderr=PIPE)
    exit_code = p.wait()
    stdout, stderr = p.communicate()
    if isinstance(stdout, bytes):
        stdout = stdout.decode()
    if isinstance(stderr, bytes):
        stderr = stderr.decode()
    return exit_code, stdout, stderr


def exc_exists(exc_name):
    """Check if an executable exists

    Args:
        exc_name (str): Executable name or path (e.g. "blastn")

    Returns:
        bool: Does the executable exists in the user's $PATH?
    """
    cmd = ['which', exc_name]
    exit_code, stdout, stderr = run_command(cmd)
    if exit_code == 0:
        return True
    else:
        logging.warning('which exited with non-zero code {} with command "{}"'.format(exit_code, ' '.join(cmd)))
        logging.warning(stderr)
        return False


def genome_name_from_fasta_path(fasta_path):
    """Extract genome name from fasta filename

    Get the filename without directory and remove the file extension.

    Example:
        With fasta file path ``/path/to/genome_1.fasta``::

            fasta_path = '/path/to/genome_1.fasta'
            genome_name = genome_name_from_fasta_path(fasta_path)
            print(genome_name)
            # => "genome_1"

    Args:
        fasta_path (str): fasta file path

    Returns:
        str: genome name
    """
    filename = os.path.basename(fasta_path)
    return re.sub(r'(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)|(\.\w{1,}$)', '', filename)


def compare_subtypes(a, b):
    for x, y in zip(a, b):
        if x != y:
            return False
    return True


def find_inconsistent_subtypes(subtypes):
    from collections import Counter
    incon = []
    for i in range(len(subtypes) - 1):
        a = subtypes[i]
        for j in range(i + 1, len(subtypes)):
            b = subtypes[j]
            is_consistent = compare_subtypes(a, b)
            if not is_consistent:
                incon.append((a, b))
    l = []
    for a, b in incon:
        astr = '.'.join([str(x) for x in a])
        bstr = '.'.join([str(x) for x in b])
        l += [astr, bstr]
    c = Counter(l)
    incon_subtypes = []
    for subtype, freq in c.most_common():
        if freq > 1:
            incon_subtypes.append(subtype)
        else:
            break
    return incon_subtypes
