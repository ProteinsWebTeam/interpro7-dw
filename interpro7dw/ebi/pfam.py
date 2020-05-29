# -*- coding: utf-8 -*-

import json
import math
import re
from io import StringIO

import MySQLdb
import MySQLdb.cursors

from interpro7dw.utils import url2dict

p7H_NTRANSITIONS = 7
eslUNKNOWN = 0
eslRNA = 1
eslDNA = 2
eslAMINO = 3
LOG2E = math.log2(math.e)


class Alphabet(object):
    def __init__(self, s):
        s = s.lower()

        if s == 'rna':
            self._type = eslRNA
            self.symbols = 'ACGU'
        elif s == 'dna':
            self._type = eslDNA
            self.symbols = 'ACGT'
        elif s == 'amino':
            self._type = eslAMINO
            self.symbols = 'ACDEFGHIKLMNPQRSTVWY'
        else:
            self._type = eslUNKNOWN
            self.symbols = ''

        self.size = len(self.symbols)

    @property
    def type(self):
        return ['unk', 'rna', 'dna', 'aa'][self._type]

    def get_bg(self):
        # null (background) model
        if self._type == eslAMINO:
            return [
                .0787945,  # A
                .0151600,  # C
                .0535222,  # D
                .0668298,  # E
                .0397062,  # F
                .0695071,  # G
                .0229198,  # H
                .0590092,  # I
                .0594422,  # K
                .0963728,  # L
                .0237718,  # M
                .0414386,  # N
                .0482904,  # P
                .0395639,  # Q
                .0540978,  # R
                .0683364,  # S
                .0540687,  # T
                .0673417,  # V
                .0114135,  # W
                .0304133  # Y
            ]
        else:
            return [1 / self.size for _ in range(self.size)]


class HMMFile(object):
    def __init__(self, stream):
        self._fh = stream
        self._reader = self._peek()

        self.length = None

        self.abc = None

        self.map = None
        self.ins = None
        self.t = None

        self.map = None
        self.consensus = None
        self.rf = None
        self.mm = None
        self.cs = None

        # Flags
        self._flag_map = False
        self._flag_cons = False
        self._flag_rf = False
        self._flag_mm = False
        self._flag_cs = False

        self.read()

    @property
    def ali_map(self):
        return [self.map[i] for i in range(1, self.length)]

    def _peek(self):
        line = next(self._fh)

        if re.match(r'HMMER3/[abcdef]', line) is not None:
            return self._parse_hmmer3
        elif re.match(r'HMMER2.0/[abcdef]', line) is not None:
            return None
        else:
            return None

    @staticmethod
    def _parse_header_line(line):
        try:
            tag, rest = line.strip().split(' ', 1)
        except ValueError:
            raise
        else:
            return tag.strip(), rest.strip()

    def _create_body(self):
        m = self.length
        k = self.abc.size

        self.mat = [[1.0 for _ in range(k)] for _ in range(m + 1)]
        self.ins = [[1.0 for _ in range(k)] for _ in range(m + 1)]
        self.t = [[1.0 for _ in range(k)] for _ in range(m + 1)]

        if self._flag_map:
            self.map = [0 for _ in range(m + 1)]

        if self._flag_cons:
            self.consensus = [None for _ in range(m + 2)]

        if self._flag_rf:
            self.rf = [None for _ in range(m + 2)]

        if self._flag_mm:
            self.mm = [None for _ in range(m + 2)]

        if self._flag_cs:
            self.cs = [None for _ in range(m + 2)]

    def _parse_hmmer3(self):
        # Parse header
        while True:
            tag, val = self._parse_header_line(next(self._fh))

            if tag == 'LENG':
                self.length = int(val)
            elif tag == 'ALPH':
                self.abc = Alphabet(val)
            elif tag == 'MAP':
                self._flag_map = val.lower() == 'yes'
            elif tag == 'CONS':
                self._flag_cons = val.lower() == 'yes'
            elif tag == 'RF':
                self._flag_rf = val.lower() == 'yes'
            elif tag == 'MM':
                self._flag_mm = val.lower() == 'yes'
            elif tag == 'CS':
                self._flag_cs = val.lower() == 'yes'

            if tag == 'HMM':
                break

        # Skip line after HMM, and allocate body of HMM
        next(self._fh)
        self._create_body()

        line = next(self._fh)
        tag, val = self._parse_header_line(line)

        if tag == 'COMPO':
            """
            Skip model composition (optional): 
            select next line (1st line of main model)
            """
            line = next(self._fh)

        """
        # First two lines are node 0: 
        insert emissions, then transitions from node 0
        """
        arr = line.strip().split()
        for x in range(0, self.abc.size):
            if arr[x] == '*':
                self.ins[0][x] = 0
            else:
                self.ins[0][x] = math.exp(-float(arr[x]))
        arr = next(self._fh).strip().split()
        for x in range(p7H_NTRANSITIONS):
            if arr[x] == '*':
                self.t[0][x] = 0
            else:
                self.t[0][x] = math.exp(-float(arr[x]))

        # Main model
        for k in range(1, self.length + 1):
            # Match emission
            arr = next(self._fh).strip().split()

            if k != int(arr[0]):  # verifies the node number
                raise RuntimeError('')

            x = 0
            for x in range(self.abc.size):
                if arr[x + 1] == '*':
                    self.mat[k][x] = 0
                else:
                    self.mat[k][x] = math.exp(-float(arr[x + 1]))

            if self._flag_map:
                self.map[k] = int(arr[x + 2])

            if self._flag_cons:
                self.consensus[k] = arr[x + 3]

            if self._flag_rf:
                self.rf[k] = arr[x + 4]

            if self._flag_mm:
                self.mm[k] = arr[x + 5]

            if self._flag_cs:
                self.cs[k] = arr[x + 6]

                # Insert emission
            arr = next(self._fh).strip().split()
            for x in range(self.abc.size):
                if arr[x] == '*':
                    self.ins[k][x] = 0
                else:
                    self.ins[k][x] = math.exp(-float(arr[x]))

            # State transition
            arr = next(self._fh).strip().split()
            for x in range(p7H_NTRANSITIONS):
                if arr[x] == '*':
                    self.t[k][x] = 0
                else:
                    self.t[k][x] = math.exp(-float(arr[x]))

        line = next(self._fh)  # todo: check the line is //

    def read(self):
        if self._reader is None:
            raise RuntimeError('invalid format')

        self._reader()

    def logo_max_height(self):
        bg = self.abc.get_bg()

        min_p = 1
        for i in range(self.abc.size):
            min_p = min(min_p, bg[i])

        return LOG2E * math.log(1 / min_p)

    def relative_entropy_all(self):
        bg = self.abc.get_bg()
        m = self.length + 1
        k = self.abc.size
        logodds = 0

        probs = [[1.0 for _ in range(k)] for _ in range(m)]
        heights = [[1.0 for _ in range(k)] for _ in range(m)]
        rel_ents = [0 for _ in range(m)]

        for i in range(1, m):
            rel_ents[i] = 0
            for j in range(k):
                probs[i][j] = self.mat[i][j]

                if probs[i][j] > 0:
                    logodds = LOG2E * math.log(probs[i][j] / bg[j])
                    rel_ents[i] += probs[i][j] * logodds

            for j in range(k):
                heights[i][j] = rel_ents[i] * probs[i][j]

        return heights[1:], probs[1:]

    def indel_value(self):
        m = self.length + 1
        insert_p = [0 for _ in range(m)]
        insert_expl = [0 for _ in range(m)]
        occupancy = [0 for _ in range(m)]

        insert_p[1] = self.t[1][1]
        insert_expl[1] = 1 / (1 - self.t[1][4])
        occupancy[1] = self.t[0][0] + self.t[0][1]

        for i in range(2, m):
            insert_p[i] = self.t[i][1]
            insert_expl[i] = 1 / (1 - self.t[i][4])
            occupancy[i] = occupancy[i - 1] * (
                    self.t[i - 1][0] + self.t[i - 1][1]
            ) + (1 - occupancy[i - 1]) * self.t[i - 1][5]

        insert_p[self.length] = 0
        insert_expl[self.length] = 0

        return insert_p[1:], insert_expl[1:], occupancy[1:]

    def get_mm(self):
        return [
            1 if self.mm is not None and self.mm[i] == 'm' else 0
            for i in range(1, self.length + 1)
        ]

    def get_alignment_map(self):
        return [
            self.map[i] for i in range(1, self.length + 1)
        ]


def hmm_to_logo(stream, method='info_content_all', processing='observed'):
    hmm = HMMFile(stream)
    max_height_theoretical = 0
    max_height_observed = 0
    min_height_observed = 0

    if method == 'info_content_all':
        max_height_theoretical = hmm.logo_max_height()
        heights, probs = hmm.relative_entropy_all()
    elif method == 'info_content_above':
        raise NotImplementedError()
    elif method == 'score':
        raise NotImplementedError()
    else:
        raise RuntimeError()

    for i, row in enumerate(heights):
        char_heights = {}
        height_sum = 0
        neg_height_sum = 0

        for h, symbol in zip(row, hmm.abc.symbols):
            if h > 0:
                height_sum += h
            else:
                neg_height_sum += h

            char_heights[symbol] = h

        if height_sum > max_height_observed:
            max_height_observed = height_sum

        if neg_height_sum < min_height_observed:
            min_height_observed = neg_height_sum

        heights[i] = [
            f"{k}:{char_heights[k]:.3f}"
            for k in sorted(char_heights, key=lambda k: char_heights[k])
        ]

    for i, row in enumerate(probs):
        char_probs = {}

        for p, symbol in zip(row, hmm.abc.symbols):
            char_probs[symbol] = p

        probs[i] = [
            f"{k}:{char_probs[k]:.3f}"
            for k in sorted(char_probs, key=lambda k: char_probs[k])
        ]

    insert_p, insert_len, delete_p = hmm.indel_value()

    return {
        'alphabet': hmm.abc.type,
        'max_height_theory': max_height_theoretical,
        'max_height_obs': max_height_observed,
        'min_height_obs': min_height_observed,
        'height_arr': heights,
        'insert_probs': insert_p,
        #'insert_probs': list(map('{:.2f}'.format, insert_p)),
        'insert_lengths': insert_len,
        # 'insert_lengths': [round(v, 1) for v in insert_len],
        'delete_probs': delete_p,
        #'delete_probs': list(map('{:.2f}'.format, delete_p)),
        'mmline': hmm.get_mm(),
        'ali_map': hmm.get_alignment_map(),
        'height_calc': method,
        'processing': processing,
        'probs_arr': probs
    }


def get_annotations(url: str):
    con = MySQLdb.connect(**url2dict(url))
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT pfamA_acc, hmm
        FROM pfamA_HMM
        WHERE hmm IS NOT NULL
        """
    )
    for accession, hmm in cur:
        yield (
            accession,
            "hmm",
            hmm,
            "text/plain",
            None  # number of sequences
        )

        # Generate logo from HMM
        with StringIO(hmm.decode()) as stream:
            logo = hmm_to_logo(stream,
                               method="info_content_all",
                               processing="hmm")

        yield (
            accession,
            "logo",
            json.dumps(logo),
            "application/json",
            None  # number of sequences
        )

    cur.execute(
        """
        SELECT pfamA_acc, num_seed, num_full, number_rp15, number_rp35, 
               number_rp55, number_rp75, number_uniprot, number_ncbi, 
               number_meta
        FROM pfamA
        """
    )
    counts = {}
    for row in cur:
        counts[row[0]] = {
            "seed": row[1],
            "full": row[2],
            "rp15": row[3],
            "rp35": row[4],
            "rp55": row[5],
            "rp75": row[6],
            "uniprot": row[7],
            "ncbi": row[8],
            "meta": row[9]
        }

    cur.execute(
        """
        SELECT pfamA_acc, type, alignment
        FROM alignment_and_tree
        WHERE alignment IS NOT NULL
        """
    )
    for accession, aln_type, aln_bytes in cur:
        yield (
            accession,
            aln_type,
            aln_bytes,  # gzip-compressed steam
            "application/gzip",
            counts[accession][aln_type]
        )
    cur.close()
    con.close()
