import math
import re


HMMFILE_FMT = {
    "2": 0,
    "3a": 1,
    "3b": 2,
    "3c": 3,
    "3d": 4,
    "3e": 5,
    "3f": 6,
}
TRANSITIONS = {
    "MM": 0,
    "MI": 1,
    "MD": 2,
    "IM": 3,
    "II": 4,
    "DM": 5,
    "DD": 6
}
NTRANSITIONS = len(TRANSITIONS)
UNKNOWN = 0
RNA = 1
DNA = 2
AMINO = 3
LOG2E = math.log2(math.e)


class Alphabet:
    def __init__(self, s):
        s = s.lower()

        if s == "rna":
            self.type = RNA
            self.symbols = "ACGU"
        elif s == "dna":
            self.type = DNA
            self.symbols = "ACGT"
        elif s == "amino":
            self.type = AMINO
            self.symbols = "ACDEFGHIKLMNPQRSTVWY"
        else:
            self.type = UNKNOWN
            self.symbols = ''

        self.K = len(self.symbols)

    def get_type(self):
        return ['unk', 'rna', 'dna', 'aa'][self.type]

    def get_bg(self):
        # See p7_bg_Create <p7_bg.c>

        if self.type == AMINO:
            # See p7_AminoFrequencies <hmmer.c>
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
            return [1 / self.K] * self.K


class HMMFile:
    def __init__(self, stream):
        # See typedef struct p7_hmm_s <hmmer.h>
        self.fh = stream
        self.format = None
        self.reader = self.peek()

        self.acc = None  # Accession number
        self.name = None  # Model name

        self.M = None  # length of model
        self.t = None  # transition prob's
        self.mat = None  # match emissions
        self.ins = None  # insert emissions
        self.abc = None  # Alphabet object
        self.rf = None  # reference line from alignment
        self.mm = None  # model mask line from alignment
        self.consensus = None  # consensus residue line
        self.cs = None  # consensus structure line
        self.map = None  # map of alignment cols onto model

        # Flags
        self.flag_rf = False  # RF annotation available
        self.flag_mmask = False  # MM annotation available
        self.flag_cons = False  # consensus residue line available
        self.flag_cs = False  # CS annotation available
        self.flag_map = False  # alignment map is available

        self.read()

    def peek(self):
        line = next(self.fh)

        m = re.match(r'HMMER3/([abcdef])', line)
        if m:
            self.format = HMMFILE_FMT["3" + m.group(1)]
            return self.parse_hmmer3
        elif re.match(r'HMMER2.0', line) is not None:
            self.format = HMMFILE_FMT["2"]
            return None
        else:
            return None

    @staticmethod
    def parse_header_line(line):
        try:
            tag, rest = line.strip().split(' ', 1)
        except ValueError:
            raise
        else:
            return tag.strip(), rest.strip()

    def create_body(self):
        # See p7_hmm_CreateBody <p7_hmm.c>
        self.t = [[0] * NTRANSITIONS for _ in range(self.M + 1)]
        self.mat = [[0] * self.abc.K for _ in range(self.M + 1)]
        self.ins = [[0] * self.abc.K for _ in range(self.M + 1)]

        # Allocation for flag-related variables
        self.rf = [None] * (self.M + 2)
        self.mm = [None] * (self.M + 2)
        self.consensus = [None] * (self.M + 2)
        self.cs = [None] * (self.M + 2)
        # self.ca = [None] * (self.M + 2)  # not initialized
        self.map = [None] * (self.M + 1)

    def parse_hmmer3(self):
        # See read_asc30hmm <p7_hmmfile.c>

        null = re.compile(r"(\*|0\.0+)$")

        # Parse header
        while True:
            tag, val = self.parse_header_line(next(self.fh))

            if tag == "ACC":
                self.acc = val
            elif tag == "NAME":
                self.name = val
            elif tag == 'LENG':
                self.M = int(val)
            elif tag == 'ALPH':
                self.abc = Alphabet(val)
            elif tag == 'RF':
                self.flag_rf = val.lower() == 'yes'
            elif tag == 'MM':
                self.flag_mmask = val.lower() == 'yes'
            elif tag == 'CONS':
                self.flag_cons = val.lower() == 'yes'
            elif tag == 'CS':
                self.flag_cs = val.lower() == 'yes'
            elif tag == 'MAP':
                self.flag_map = val.lower() == 'yes'

            if tag == 'HMM':
                break

        # Skip line after HMM, and allocate body of HMM
        next(self.fh)
        self.create_body()

        line = next(self.fh)
        tag, val = self.parse_header_line(line)

        if tag == 'COMPO':
            """
            Skip model composition (optional):
            select next line (1st line of main model)
            """
            line = next(self.fh)

        """
        First two lines are node 0:
        insert emissions, then transitions from node 0
        """
        values = line.strip().split()
        for x in range(self.abc.K):
            val = values[x]
            # self.ins[0][x] = 0 if val == '*' else math.exp(-float(val))
            self.ins[0][x] = 0 if null.match(val) else math.exp(-float(val))

        values = next(self.fh).strip().split()
        for x in range(NTRANSITIONS):
            val = values[x]
            # self.t[0][x] = 0 if val == '*' else math.exp(-float(val))
            self.t[0][x] = 0 if null.match(val) else math.exp(-float(val))

        # Main model
        for k in range(1, self.M + 1):
            # Match emission line
            values = iter(next(self.fh).strip().split())

            # Node number check
            val = next(values)
            if int(val) != k:
                raise ValueError(
                    f"Expected match line to start with {k} (of {self.M}); saw {val}")

            # Match emissions
            for x in range(self.abc.K):
                val = next(values)
                # self.mat[k][x] = 0 if val == '*' else math.exp(-float(val))
                self.mat[k][x] = 0 if null.match(val) else math.exp(
                    -float(val))

            # MAP annotation
            val = next(values)
            if self.flag_map:
                self.map[k] = int(val)

            if self.format >= HMMFILE_FMT["3e"]:
                # CONS consensus residue
                val = next(values)
                if self.flag_cons:
                    self.consensus[k] = val

            # RF annotation
            val = next(values)
            if self.flag_rf:
                self.rf[k] = val

            if self.format >= HMMFILE_FMT["3f"]:
                # MM mask value
                val = next(values)
                if self.flag_mmask:
                    self.mm[k] = val

            # CS annotation
            val = next(values)
            if self.flag_cs:
                self.cs[k] = val

            # Insert emission line
            values = iter(next(self.fh).strip().split())
            for x in range(self.abc.K):
                val = next(values)
                # self.ins[k][x] = 0 if val == '*' else math.exp(-float(val))
                self.ins[k][x] = 0 if null.match(val) else math.exp(
                    -float(val))

            # State transition line
            values = iter(next(self.fh).strip().split())
            for x in range(NTRANSITIONS):
                val = next(values)
                # self.t[k][x] = 0 if val == '*' else math.exp(-float(val))
                self.t[k][x] = 0 if null.match(val) else math.exp(-float(val))

        val = next(self.fh).strip()
        if val != "//":
            raise ValueError(f"Expected closing //; found {val} instead")

    def read(self):
        if self.reader is None:
            raise RuntimeError('invalid format')

        self.reader()

    """
    Below are functions for logos
    """

    def logo_max_height(self):
        # See hmmlogo_maxHeight <hmmlogo.c>
        bg = self.abc.get_bg()

        min_p = 1
        for i in range(self.abc.K):
            min_p = min(min_p, bg[i])

        return LOG2E * math.log(1 / min_p)

    def relative_entropy_all(self):
        # hmmlogo_RelativeEntropy_all <hmmlogo.c>
        bg = self.abc.get_bg()

        rel_ents = [0] * (self.M + 1)
        heights = [[0] * self.abc.K for _ in range(self.M + 1)]
        probs = [[0] * self.abc.K for _ in range(self.M + 1)]

        logodds = 0
        for i in range(1, self.M + 1):
            # height of column, to be split among the residues
            rel_ents[i] = 0
            for j in range(self.abc.K):
                probs[i][j] = self.mat[i][j]

                if probs[i][j] > 0:
                    logodds = LOG2E * math.log(probs[i][j] / bg[j])
                    rel_ents[i] += probs[i][j] * logodds

            # height of residues
            for j in range(self.abc.K):
                heights[i][j] = rel_ents[i] * probs[i][j]

        return heights[1:], probs[1:]

    def indel_value(self):
        # See hmmlogo_IndelValues <hmmlogo.c>

        # probability of inserting after this match
        insert_p = [0] * (self.M + 1)
        # expected length of the insert, if it happens
        insert_expl = [0] * (self.M + 1)

        mm = TRANSITIONS["MM"]
        mi = TRANSITIONS["MI"]
        ii = TRANSITIONS["II"]
        dm = TRANSITIONS["DM"]

        for i in range(1, self.M):
            insert_p[i] = self.t[i][mi]
            insert_expl[i] = 1 / (1 - self.t[i][ii])

        occupancy = [0] * (self.M + 1)
        occupancy[1] = self.t[0][mi] + self.t[0][mm]
        for k in range(2, self.M + 1):
            occupancy[k] = occupancy[k - 1] * (
                        self.t[k - 1][mm] + self.t[k - 1][mi]) + (
                                       1.0 - occupancy[k - 1]) * self.t[k - 1][
                               dm]

        return insert_p[1:], insert_expl[1:], occupancy[1:]

    def get_mm(self):
        mm = []
        for i in range(1, self.M + 1):
            if self.mm[i] == 'm':
                mm.append(1)
            else:
                mm.append(0)

        return mm

    def get_alignment_map(self):
        return [self.map[i] for i in range(1, self.M + 1)]

    def logo(self, method="info_content_all", processing="observed"):
        max_height_theoretical = 0
        max_height_observed = 0
        min_height_observed = 0
        if method == 'info_content_all':
            max_height_theoretical = self.logo_max_height()
            heights, probs = self.relative_entropy_all()
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

            for h, symbol in zip(row, self.abc.symbols):
                if h > 0:
                    height_sum += h
                else:
                    neg_height_sum += h

                char_heights[symbol] = h

            if height_sum > max_height_observed:
                max_height_observed = height_sum

            if neg_height_sum < min_height_observed:
                min_height_observed = neg_height_sum

            # Sort by height
            heights[i] = [
                f"{k}:{char_heights[k]:.3f}"
                for k in sorted(char_heights, key=lambda k: char_heights[k])
            ]

        for i, row in enumerate(probs):
            char_probs = {}

            for p, symbol in zip(row, self.abc.symbols):
                char_probs[symbol] = p

            probs[i] = [
                f"{k}:{char_probs[k]:.3f}"
                for k in sorted(char_probs, key=lambda k: char_probs[k])
            ]

        insert_p, insert_len, delete_p = self.indel_value()

        return {
            'alphabet': self.abc.get_type(),
            'max_height_theory': max_height_theoretical,
            'max_height_obs': max_height_observed,
            'min_height_obs': min_height_observed,
            'height_arr': heights,
            'insert_probs': insert_p,
            # 'insert_probs': list(map('{:.2f}'.format, insert_p)),
            'insert_lengths': insert_len,
            # 'insert_lengths': [round(v, 1) for v in insert_len],
            'delete_probs': delete_p,
            # 'delete_probs': list(map('{:.2f}'.format, delete_p)),
            'mmline': self.get_mm(),
            'ali_map': self.get_alignment_map(),
            'height_calc': method,
            'processing': processing,
            'probs_arr': probs
        }
