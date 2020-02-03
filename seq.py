import os
import re
import random

try:
    import numpy as np
    import pandas as pd
except:
    raise ModuleNotFoundError("Pandas or numpy module not found ")


class Protein:
    """
    protein analyzes
    mutation detection and time
    """

    __peptide = []
    __a_acids = {
        "Fenilalanin": "F",
        "Alanin": "A",
        "Arginin": "R",
        "Asparagin": "N",
        "Aspartik asit": "D",
        "Sistein": "C",
        "Glutamik asit": "E",
        "Glutamin": "Q",
        "Glisin": "G",
        "Histidin": "H",
        "Izolösin": "I",
        "Lösin": "L",
        "Lizin": "K",
        "Metiyonin": "M",
        "Prolin": "P",
        "Serin": "S",
        "Treonin": "T",
        "Triptofan": "W",
        "Trozin": "Y",
        "Valin": "V"
    }
    __particular = {
        "F": (165.19, 1.36, 283),
        "A": (89.1, 1.42, 295),
        "R": (174.2, 0.7, 238),
        "N": (132.12, "NA", 234),
        "D": (133.1, "NA", 270),
        "C": (121.16, "NA", 220),
        "E": (169.111, 1.54, 205),
        "Q": (146.15, "NA", 185),
        "G": (75.07, 1.607, 232),
        "I": (131.17, "NA", 284),
        "L": (131.18, 1.29, 293),
        "K": (146.19, "NA", 224),
        "M": (149.21, 1.34, 280),
        "P": (115.13, 1.36, 220),
        "S": (115.09, "NA", 215),
        "T": (119.12, "NA", 255),
        "W": (204.23, "NA", 290),
        "Y": (181.19, 1.46, 342),
        "V": (117.15, 1.316, 295),
        "info": ("m", "d", "en")
    }

    @classmethod
    def __control(cls, seq):
        snc = []
        kontrol = "AUGCT"
        for bs in seq:
            if bs in kontrol:
                snc.append(True)
            else:
                snc.append(False)
        return all(snc)

    @classmethod
    def __rRNA(cls, seq):
        """
        first python exp...
        :param seq: sequence
        :return:
        """
        cls.__codons = []
        rstx = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        mismatch = []
        for base in rstx:
            base_kontrol = []
            for nkltt in base:
                base_kontrol.append(cls.__control(nkltt))
            if all(base_kontrol) == True:
                if (base == "UUU" or base == "UUC"):
                    cls.__codons.append(cls.__a_acids.get("Fenilalanin"))
                elif \
                        (
                                base == "UUA" or base == "UUG" or base == "CUU" or base == "CUC" or base == "CUA" or base == "CUG"):
                    cls.__codons.append(cls.__a_acids.get("Lösin"))
                elif (base == "AUU" or base == "AUC" or base == "AUA"):
                    cls.__codons.append(cls.__a_acids.get("Izolösin"))
                elif (base == "AUG"):
                    cls.__codons.append(cls.__a_acids.get("Metiyonin"))
                elif (base == "GUU" or base == "GUC" or base == "GUA" or base == "GUG"):
                    cls.__codons.append(cls.__a_acids.get("Valin"))
                elif (base == "UCU" or base == "UCC" or base == "UCA" or base == "UCG"):
                    cls.__codons.append(cls.__a_acids.get("Serin"))
                elif (base == "CCU" or base == "CCC" or base == "CCA" or base == "CCG"):
                    cls.__codons.append(cls.__a_acids.get("Prolin"))
                elif (base == "ACU" or base == "ACC" or base == "ACA" or base == "ACG"):
                    cls.__codons.append(cls.__a_acids.get("Treonin"))
                elif (base == "GCU" or base == "GCC" or base == "GCA" or base == "GCG"):
                    cls.__codons.append(cls.__a_acids.get("Alanin"))
                elif (base == "UAU" or base == "UAC"):
                    cls.__codons.append(cls.__a_acids.get("Trozin"))
                elif (base == "UAA" or base == "UAG" or base == "UGA"):
                    if (base == "UAA"):
                        cls.__codons.append("*")
                    elif (base == "UGA"):
                        cls.__codons.append("*")
                    else:
                        cls.__codons.append("*")
                elif (base == "CAU" or base == "CAC"):
                    cls.__codons.append(cls.__a_acids.get("Histidin"))
                elif (base == "CAA" or base == "CAG"):
                    cls.__codons.append(cls.__a_acids.get("Glutamin"))
                elif (base == "AAU" or base == "AAC"):
                    cls.__codons.append(cls.__a_acids.get("Asparagin"))
                elif (base == "AAA" or base == "AAG"):
                    cls.__codons.append(cls.__a_acids.get("Lizin"))
                elif (base == "GAU" or base == "GAC"):
                    cls.__codons.append(cls.__a_acids.get("Aspartik asit"))
                elif (base == "GAA" or base == "GAG"):
                    cls.__codons.append(cls.__a_acids.get("Glutamik asit"))
                elif (base == "UGU" or base == "UGC"):
                    cls.__codons.append(cls.__a_acids.get("Sistein"))
                elif (base == "UGG"):
                    cls.__codons.append(cls.__a_acids.get("Triptofan"))
                elif (
                        base == "CGU" or base == "CGC" or base == "CGA" or base == "CGG" or base == "AGA" or base == "AGG"):
                    cls.__codons.append(cls.__a_acids.get("Arginin"))
                elif (base == "AGU" or base == "AGC"):
                    cls.__codons.append(cls.__a_acids.get("Serin"))
                elif (base == "GGU" or base == "GGC" or base == "GGA" or base == "GGG"):
                    cls.__codons.append(cls.__a_acids.get("Glisin"))
                else:
                    mismatch.append(base)
                    cls.__codons.append("-")
        return cls.__codons

    @classmethod
    def cntrl_dogma(cls, seq="ATGC", arg=""):
        """
        ------chain flip or swap----
        ------standard operations for central dogma-----
        :param seq: target sequence
        :param arg: tips(tc, tl or cm)
        :return:
        """
        if arg == "tc":
            return cls.__mRNA(seq)
        if arg == "tl":
            return cls.__orf(seq)
        if arg == "cm":
            return cls.__cmplmntr(seq)
        else:
            raise ValueError("arg in = tc : Transcription, tl : Translation, cm : Complement")

    @classmethod
    def __orf(cls, seq="AUGC"):
        rslt_str = ""
        for base in cls.__rRNA(seq):
            rslt_str += base
        return rslt_str

    @classmethod
    def __mRNA(cls, seq="ATGC"):
        mrna = ""
        for base in seq:
            if base == "T":
                mrna += "U"
            else:
                mrna += base
        return mrna

    @classmethod
    def __cmplmntr(cls, seq="ATGC"):
        cmpl = ""
        mismatch = ""
        for base in seq:
            if base == "A":
                cmpl += "T"
            elif base == "T":
                cmpl += "A"
            elif base == "G":
                cmpl += "C"
            elif base == "C":
                cmpl += "G"
            else:
                mismatch += base
        return cmpl

    @classmethod
    def __revrse_a_acids(cls, a):
        reverse_acids = {}
        for acid in cls.__a_acids.items():
            reverse_acids[acid[1]] = acid[0]
        return reverse_acids.get(a, "*")

    @classmethod
    def mutant_detection(cls, seq1="ATGCU", seq2="ATGCU", show="pep"):
        """
        -----homolog sequence for mutation detection----
        :param seq1: homolog sequence first
        :param seq2: homolog sequence second
        :param show: output kind (pep,peptide or None)
        :return: list
        """
        n = 0
        rslt = []
        seqx = cls.__mRNA(seq1)
        seqy = cls.__mRNA(seq2)
        _rstx = [seqx[i:i + 3] for i in range(0, len(seqx), 3)]
        _rstx_2 = [seqy[i:i + 3] for i in range(0, len(seqy), 3)]
        rstx = [i for i in _rstx if len(i) == 3]
        rstx_2 = [i for i in _rstx_2 if len(i) == 3]
        acid = cls.cntrl_dogma(seqx, "tl")
        acid_2 = cls.cntrl_dogma(seqy, "tl")
        while n < min(len(rstx), len(rstx_2)):
            if show == "pep":
                out_p = acid[n], acid_2[n]
            elif show == "peptide":
                out_p = cls.__revrse_a_acids(acid[n]), cls.__revrse_a_acids(acid_2[n])
            else:
                out_p = rstx[n], rstx_2[n]
            if cls.__control(rstx[n]) == True and cls.__control(rstx_2[n]) == True:
                if rstx[n] != rstx_2[n] and acid[n] == acid_2:
                    rslt.append(out_p[0] + " -> " + out_p[1] + " : s")
                if rstx[n] != rstx_2[n] and acid[n] != acid_2[n] and acid_2[n] == "*":
                    rslt.append(out_p[0] + " -> " + out_p[1] + " : n")
                elif rstx[n] != rstx_2[n] and acid[n] != acid_2[n]:
                    rslt.append(out_p[0] + " -> " + out_p[1] + " : m")
            n += 1
        if show == "info":
            return "s : Silent Mutation", "n : Nonesense Mutation", "m : Missense Mutation"
        return rslt

    @classmethod
    def mutant_maker(cls, seq="ATGC", rate=20):
        """
        --------Creates sequence mutation--------
        :param seq: target sequence
        :param rate: mutation rate
        :return: mutant sequence
        """
        s_seqx = ""
        seqx = list(str(seq))
        liste = list(range(0, len(seqx)))
        n = 0
        x = int(len(seqx) * rate / 100)
        while n < 3:
            random.shuffle(liste)
            n += 1
        for ndx in random.sample(liste, x):
            m = ["A", "T", "G", "C"]
            base = seqx[ndx]
            m.remove(base)
            m_base = random.choice(m)
            seqx.pop(ndx)
            seqx.insert(ndx, m_base)
        for bs in seqx:
            s_seqx += bs
        return s_seqx


class N_tears:
    """
    test class
    """
    def __init__(self):
        pass

    @classmethod
    def _komplex_alignment(cls, liste, liste_2):
        """
        not recommended

        :param liste:list
        :param liste_2: list
        :return:
        """
        n = 0
        m = ""
        match_list_1 = []
        match_seperator = []
        match_liste_2 = []
        for _base in liste:
            if re.match("[0-9]", _base[0]):
                object_1 = cls.backmelting(_base)
                _base_1 = object_1[0]
            if re.match("[ATGC]", _base[0]):
                _base_1 = _base
            for base_ in liste_2:
                if re.match("[0-9]", base_[0]):
                    object_2 = cls.backmelting(base_)
                    _base_2 = object_2[0]
                if re.match("[ATGC]", base_[0]):
                    _base_2 = base_
                m_1 = ""
                y = 0
                while y < len(_base_2):
                    value_1 = _base_1[y]
                    value_2 = _base_2[y]
                    if value_1 == value_2:
                        if value_1 != "-" and value_2 != "-":
                            if (cls._nucleotide_Control(value_1) == True) \
                                    and \
                                    (cls._nucleotide_Control(value_2) == True):
                                m += "|"
                                m_1 += "|"
                    if value_1 != value_2:
                        if (cls._nucleotide_Control(value_1) == True) \
                                and \
                                (cls._nucleotide_Control(value_2) == True):
                            m_1 += "*"
                    y += 1
                match_seperator.append(m_1)
                match_list_1.append(_base)
                match_liste_2.append(base_)
                m += " "

            n += 1
        return match_seperator, match_list_1, match_liste_2

    @classmethod
    def backmelting(cls, base):
        if isinstance(base, str):
            marker = base.index(" <")
            marker_2 = base.index("<= ")
            position_1 = base[:marker]
            position_2 = base[marker_2 + 2:]
            out_base = (base[marker + 2:marker_2])
            return out_base, marker + 2, marker_2, (out_base, position_1, position_2)
        if isinstance(base, list):
            result = list()
            for _base_ in base:
                marker = _base_.index(" <")
                marker_2 = _base_.index("<= ")
                position_1 = _base_[:marker]
                position_2 = _base_[marker_2 + 2:]
                out_base = (_base_[marker + 2:marker_2])
                result.append((out_base, position_1, position_2))
            return result

    @staticmethod
    def base_rate(seq):
        info = {}
        adenin = seq.count("A")
        guanin = seq.count("G")
        sitozin = seq.count("C")
        timin = seq.count("T")
        info["A"] = adenin
        info["G"] = guanin
        info["C"] = sitozin
        info["T"] = timin
        info["A/T"] = f'{adenin / timin:.3}'
        info["C/G"] = f'{sitozin / guanin:.3}'
        return info

    @staticmethod
    def _Anti_Parallel(seq1, seq2):
        out_put = ""
        r_eve_r = reversed(seq2)
        for s in r_eve_r:
            out_put += s
        if seq1 == out_put:
            return True
        else:
            return False

    @staticmethod
    def _nucleotide_Control(base, type="ATGC"):
        _control = type
        return base in _control

    @staticmethod
    def restriction(seq="MustafaUyar", splite="5", specifity=False):
        """
        to splite the sequence
        :param seq: array to be plotted
        :param splite: form of parceling
        :param specifity: attribute
        :return:list
        """
        splite = int(splite)
        non_id_rstx_list = []
        id_rstx_list = []
        upper = 0
        _Treu_seq = ""
        for base in seq:
            if base != "\n":
                _Treu_seq += base
        while upper < splite:
            restriction = [_Treu_seq[i:i + splite] for i in range(upper, len(seq), splite)]
            position = 0
            for base in restriction:
                if not specifity:
                    if len(base) == splite:
                        non_id_rstx_list.append(base)
                if specifity:
                    if len(base) == splite:
                        id_rstx_list.append((str(upper + (splite * position)) +
                                             " <"
                                             + base +
                                             "<= "
                                             + str(((splite * (position + 1)) - (splite - len(base))) + upper)))
                position += 1
            upper += 1
        if specifity:
            return id_rstx_list
        if not specifity:
            return non_id_rstx_list


class test():
    @staticmethod
    def polindromic_array(seq1, seq2, minsize=4, show="|both", specifty=True):
        """
        ------------find two sequential polyindomic strings--------
        :param seq1: first sequence to be scanned
        :param seq2: second sequence to be scanned
        :param minsize: lowest sequence length
        :param show: display preference ( '|' = preferred lock )
        :param specifty: identification
        :return:
        """

        def output_regulator(s, x):
            if show == "seq1" or show == "|" + "seq1":
                result.append(s)
            if show == "seq2" or show == "|" + "seq2":
                result.append(x)
            if show == "both" or show == "|" + "both":
                ekle = s, x
                result.append(ekle)

        result = list()
        n = 0
        mediator = min(len(seq1), len(seq2))
        while n < mediator:
            rstx_1 = Sequence.rstx_smash(seq1, mediator - n, specifty)
            rstx_2 = Sequence.rstx_smash(seq2, mediator - n, specifty)
            for base_1 in rstx_1:
                for base_2 in rstx_2:
                    _base_1 = base_1
                    _base_2 = base_2
                    if specifty == True:
                        _base_1 = N_tears().backmelting(base_1)[3][0]
                        _base_2 = N_tears().backmelting(base_2)[3][0]
                    query = N_tears()._Anti_Parallel(_base_1, _base_2)
                    if query == True and _base_1 != _base_2:
                        if "|" == show[0]:
                            if len(_base_1) == minsize and len(_base_2) == minsize:
                                output_regulator(base_1, base_2)
                        else:
                            if len(_base_1) >= minsize and len(_base_2) >= minsize:
                                output_regulator(base_1, base_2)
            n += 1
        return result


class Sequence(test, N_tears):
    _dir_name = "random_sequence_id%s.txt"
    __random_seq_folder_id = 0
    rsfid = _random_seq_folder_id_dict_in_the_memory = dict()

    def __init__(self):
        super().__init__()

    @classmethod
    def show(cls):
        return cls._random_seq_folder_id_dict_in_the_memory

    @staticmethod
    def revers(seq):
        return "".join([i for i in reversed(seq)])

    @staticmethod
    def repeat_array(sequence, min_unit=2, maxperlength=14, specifity=False):
        """
        -------------find repeating motifs----------
        :param sequence: sequence to be scanned
        :param min_unit: smallest repeat unit
        :param maxperlength: the highest array size
        :param specifity: identification
        :return:
        """
        result = []
        for base in Sequence.rstx_smash(sequence, maxperlength, specifity):
            base_1 = base
            if specifity:
                base_1 = N_tears().backmelting(base)[3][0]
            n = min_unit
            while n < len(base_1):
                match = set([base_1[x:x + n] for x in range(0, len(base_1), n)])
                if len(base_1) == base_1.count("A") or \
                        len(base_1) == base_1.count("T") or \
                        len(base_1) == base_1.count("G") or \
                        len(base_1) == base_1.count("C"):
                    pass
                else:
                    if len(match) == 1:
                        add = base, n
                        result.append(add)
                n += 1
        return result

    @staticmethod
    def snp_array_(seq1="main", seq2="side", difference="seq1"):
        """
        ----------the only nucleotide or array change between two sequences is found---------
        not : offers convenient use for homolog sequence
        :param seq1: first main sequence to be scanned
        :param seq2: second side sequence to be scanned
        :param difference: difference preference
        :return:
        """
        PUBG = "SNP"
        o_u_t__p_u_t = {}

        def diff(difference):
            if difference == "both":
                o_u_t__p_u_t["both"] = list(_rstx_.symmetric_difference(_rstx_2))
            if difference == "seq1":
                o_u_t__p_u_t["seq1"] = list(_rstx_ - _rstx_2)
            if difference == "seq2":
                o_u_t__p_u_t["seq2"] = list(_rstx_2 - _rstx_)
            return o_u_t__p_u_t

        if PUBG == "SNP":
            size = 1
            specifıty = True
            _rstx_ = set(Sequence.rstx_smash(seq1, size, specifıty))
            _rstx_2 = set(Sequence.rstx_smash(seq2, size, specifıty))
            return diff(difference)
        return False

    @staticmethod
    def common_array_(seq1, seq2, PUBG="solo", size=1, sensitive=False):
        """
        -------------find common motifs between the two series------------
        not : offers convenient use for heterolog sequence
        :param seq1:  first sequence to be scanned
        :param seq2: second sequence to be scanned
        :param PUBG: ( solo, duo, squad ) scanning preference
        :param size: PUBG="solo" for range
        :param sensitive: scanning sensitivity
        :return:
        """
        if PUBG == "solo":
            _rstx_ = set(Sequence.rstx_smash(seq1, size, sensitive))
            _rstx_2 = set(Sequence.rstx_smash(seq2, size, sensitive))
            if not _rstx_.isdisjoint(_rstx_2):
                return _rstx_ & _rstx_2
            if _rstx_.isdisjoint(_rstx_2):
                return None
        if PUBG == "duo" or PUBG == "squad":
            result = []
            n = 0
            mediator = min(len(seq1), len(seq2))
            while n < mediator:
                _rstx_ = set(Sequence.rstx_smash(seq1, mediator - n, sensitive))
                _rstx_2 = set(Sequence.rstx_smash(seq2, mediator - n, sensitive))
                if not _rstx_.isdisjoint(_rstx_2):
                    result.append((_rstx_ & _rstx_2))
                    if PUBG == "duo":
                        return _rstx_ & _rstx_2
                    if PUBG == "duo":
                        break
                n += 1
            if result.__len__() == 0:
                return None
            return result
        return False

    @classmethod
    def random_seq(cls, len, cache="", redo=1):
        sekans = random.sample("ATGC" * len, len)
        random.shuffle(sekans)
        sequence = "".join(sekans)
        cls.__random_seq_folder_id += 1
        cls._dir = cls._dir_name % (str(cls.__random_seq_folder_id))
        if cache == "::":
            cls._random_seq_folder_id_dict_in_the_memory[str(cls.__random_seq_folder_id)] = str(sequence)
        if cache == "":
            try:
                os.remove(cls._dir)
            except:
                pass
            return sequence
        if redo == 1 and cache == ":":
            f = open(cls._dir, "w")
            print("{}".format(str(sequence)), file=f, flush=True)
        if cache == ":":
            if not os.path.exists(cls._dir):
                f = open(cls._dir, "w")
                print("{}".format(str(sequence)), file=f, flush=True)
            else:
                pass
        if not cache == "::":
            dosya = open(cls._dir, "r")
            oku = dosya.read()
            outer = ""
            for i in oku:
                if i != "\n":
                    outer += i
            return outer

    @classmethod
    def clear(cls, ignore="", stress=None):
        if stress is None:
            select = cls.__random_seq_folder_id
        else:
            select = int(stress)
        for i in range(select):
            if str(i + 1) not in list(ignore):
                try:
                    die = cls._dir_name % (str(i + 1))
                    os.remove(die)
                except:
                    pass

    @classmethod
    def open_random_seq(cls, ignore=""):
        top_seq = []
        for i in range(cls.__random_seq_folder_id):
            if str(i + 1) not in list(ignore):
                outer = ""
                direc_ = cls._dir_name % (str(i + 1))
                seq_open = open(direc_, "r")
                for base in seq_open.read():
                    if base != "\n":
                        outer += base
                top_seq.append(outer)
        return top_seq

    @staticmethod
    def rstx_smash(sequence, splite, specifity=False):
        """
        ----------- assigns an ID by truncating the sequence---------------
        :param sequence: sequence to be identified
        :param splite:  amount of identification
        :param specifity: identity preference
        :return:
        """
        return N_tears.restriction(sequence, splite, specifity)

    @staticmethod
    def simple_align(list=list("ATGC"), list2=list("ATGC"), kind="dflt", rate=75):
        """
        ---------- !! CALL NOT RECOMMENDED !!----------
        ----------development process in progress--------
        ---------- is used by some functions in the shell --------
        :param list:
        :param list2:
        :param kind:
        :param rate:
        :return:
        """

        def puanla_motif():
            if kind == "list":
                n = 0
                output = []
                for s in a[0]:
                    rate_1 = a[0][n].count("|")
                    rate_2 = a[0][n].count("*")
                    anlamsal_oran = rate_1 / (rate_1 + rate_2 + (1 * 10 ** (-49))) * 100
                    if anlamsal_oran >= int(rate) and a[1][n] != "-" * len(a[1][n]) and a[2][n] != "-" * len(a[2][n]):
                        output.append((a[2][n], int(rate_1 / (rate_1 + rate_2 + 1 * 10 ** (-49)) * 100)))
                    n += 1
                return output

        a = N_tears()._komplex_alignment(list, list2)
        if kind == "list":
            return puanla_motif()
        if kind == "dflt":
            return zip(a[1], a[0], a[2])

    @staticmethod
    def find_motif(sequence, motif, rate=99, specifity=False):
        """
        -----------find motifs on the layout--------
        :param motif: motif
        :param sequence: sequence to be scanned
        :param rate: purity percentage
        :param specifity: identification
        :return:
        """
        motif = motif.__str__()
        motif_alfa = []
        motif_alfa.append(motif)
        result = Sequence.simple_align(motif_alfa, Sequence.rstx_smash(sequence, len(motif), specifity), rate=rate,
                                       kind="list")
        # if kind == "df":
        #    df = pd.DataFrame(np.array(list(result)), index=list(len(result)*["find_motif"]), columns=["sequence", "rate"])
        #    return df
        return result

    @staticmethod
    def imperfect_seq(sequence, motif1="TATAA", motif2="TATAA", i_volume=1, i_seq="|GC", specifity=True):
        """
        ------------it finds leaking lines-----------
        :param sequence:   sequence to be scanned
        :param motif1: first motif
        :param motif2: second motif
        :param i_volume: infiltration amount
        :param i_seq: Determining the desired leaky together knee ( '|' = preferred lock )
        :param specifity: identification
        :return:
        """
        nfrtt = N_tears()
        o_u_t__p_u_t = []
        _rstx = Sequence.rstx_smash(sequence + " " * i_volume, splite=len(motif1) + len(motif2) + i_volume,
                                    specifity=specifity)
        if specifity:
            _melting = [nfrtt.backmelting(i)[3] for i in _rstx]
            rstx_ = _melting
            for base in rstx_:
                if i_seq[0] == "|":
                    objct_ = re.match("(%s)[%s]{%s}(%s)" % (motif1, i_seq[0:], i_volume, motif2), base[0])
                else:
                    objct_ = re.match("(%s)[%s]+(%s)" % (motif1, i_seq, motif2), base[0])
                if objct_:
                    add_ = (base[1], objct_.group(), base[2])
                    o_u_t__p_u_t.append(add_)
        if not specifity:
            for base in _rstx:
                if i_seq[0] == "|":
                    _objct = re.match("(%s)[%s]{%s}(%s)" % (motif1, i_seq[0:], i_volume, motif2), base)
                else:
                    _objct = re.match("(%s)[%s]+(%s)" % (motif1, i_seq, motif2), base)
                if _objct:
                    _add = _objct.group()
                    o_u_t__p_u_t.append(_add)
        return o_u_t__p_u_t

    @staticmethod
    def distance_seq(sequence, motif1="TATAA", motif2="TATAA", kind="list", defect=None):
        """
        -----------measure distance between two motifs-----------
        :param sequence: sequence to be scanned
        :param motif1: first motif
        :param motif2: second motif
        :param kind: converter (list,dict, df = dataframe
        :param defect: defect
        :return:
        """
        kind = str(kind.lower())
        results = []
        for _base in [N_tears().backmelting(i)[3] for i in Sequence.rstx_smash(sequence, len(motif1), True)]:
            if _base[0] == motif1:
                meaningful_value = [(str(_base[1]) + " <" + _base[0] + "<=" + str(_base[2]))]
                for base_ in [N_tears().backmelting(i)[3] for i in Sequence.rstx_smash(sequence, len(motif2), True)]:
                    if base_[0] == motif2:
                        if defect is None:
                            if _base[0][-1] == base_[0][0] or _base[0][0] == base_[0][-1]:
                                defect = False
                            else:
                                defect = True
                        matrix_1 = int(_base[2]) - int(base_[1])
                        matrix_2 = int(_base[1]) - int(base_[2])
                        if not defect:
                            dfct_1 = int(_base[1]) - int(base_[1])
                            dfct_2 = int(_base[2]) - int(base_[2])
                            dfct_R = max(abs(dfct_1), abs(dfct_2))
                            dfct_B = max(len(_base[0]), len(base_[0]))
                            if dfct_R < dfct_B:
                                if _base[0] == base_[0] and _base[1] == base_[1] and _base[2] == base_[2]:
                                    values_ = (base_[0], "match")
                                    meaningful_value.append(values_)
                                else:
                                    values_ = (base_[0], "inside")
                                    meaningful_value.append(values_)
                            else:
                                if matrix_1 < 0 and matrix_2 < 0:
                                    negative_matrix_main = max(matrix_1, matrix_2)
                                    values_ = (base_[0], negative_matrix_main)
                                    meaningful_value.append(values_)
                                else:
                                    positive_matrix_main = min(abs(matrix_1), abs(matrix_2))
                                    values_ = (base_[0], positive_matrix_main)
                                    meaningful_value.append(values_)
                        if defect:
                            if matrix_1 < 0 and matrix_2 < 0:
                                negative_matrix_main = max(matrix_1, matrix_2)
                                values_ = (base_[0], negative_matrix_main)
                                meaningful_value.append(values_)
                            else:
                                positive_matrix_main = min(abs(matrix_1), abs(matrix_2))
                                values_ = (base_[0], positive_matrix_main)
                                meaningful_value.append(values_)
                results.append(meaningful_value)

        def convert(type="dict"):
            if type == "dict" or type == "df":
                result = {}
            if type == "list":
                result = []
            return result

        def tool(results, srch_seq, sp1=" <", sp2="<= "):
            sonuç = list()
            for i in results:
                seq, start, end = Sequence.backmelting(i[0])[3]
                m = 1
                ara_sonuç = list()
                ara_sonuç.append(i[0])
                while m < len(i):
                    variable = i[m][1]
                    if variable < 0:
                        konum = srch_seq[int(end) - int(variable): int(end) - int(variable) + len(seq)]
                        cıktı = (str(int(end) - int(variable)) + sp1 + str(konum) + sp2 + str(
                            int(end) - int(variable) + len(seq)))
                        ara_sonuç.append(cıktı)
                    elif variable > 0:
                        konum = srch_seq[int(start) - int(variable) - len(seq):int(start) - int(variable)]
                        cıktı = (str(int(start) - int(variable) - len(seq)) + sp1 + konum + sp2 + str(
                            int(start) - int(variable)))
                        ara_sonuç.append(cıktı)
                    m += 1
                sonuç.append(ara_sonuç)
            return sonuç

        def pretty():
            main = list()
            for i in results:
                ara = list()
                m = 1
                ara.append(i[0])
                while m < len(i):
                    ara.append(i[m][1])
                    m += 1
                main.append(ara)
            return main

        if kind == "sample":
            return tool(results, sequence)
        if kind == list:
            return convert(type="list")
        if kind == "prettify":
            return pretty()
        if kind == "":
            return ("use : prettify", "sample", "list")
        else:
            return results


class Series(N_tears):
    """
    purporse : Creating DNA sequences
    testing class
    """
    def __init__(self, seq, smash=5, expand=0):
        self.__expand = expand
        self.__smash = smash
        if isinstance(seq, str) == False:
            raise ValueError("only str format")
        self.__sequence = seq
        if expand == 1:
            self.__sırala = sorted(Sequence.rstx_smash(seq, smash, True),
                                   key=lambda x: int(N_tears().backmelting(x)[3][1]))
        else:
            self.__sırala = []
            upper = 0
            for base in [seq[i:i + smash] for i in range(0, len(seq), smash)]:
                self.__sırala.append(str(upper) + " <" + base + "<= " + str(upper + len(base)))
                upper += smash

    @classmethod
    def idseqsort(cls, seq):
        def control(seq):
            result = [i for i in seq
                      if i.__str__().isalnum() == False \
                      and i.__str__().isalpha() == False \
                      and i.__str__().isdigit() == False
                      ]
            return result

        if isinstance(seq[0], tuple):
            return sorted(control(seq), key=lambda x: int(N_tears.backmelting(x[0])[3][1]))
        elif isinstance(seq, list):
            return sorted(control(seq), key=lambda x: int(N_tears.backmelting(x)[3][1]))
        else:
            raise TypeError("isproperty alone list[i] = ' integer <string <= integer '\nsamples: [2 <ATGC<= 6]\n" +
                            " " * "samples".__len__() + "  [int <str<= int]")

    def __repr__(self):
        return self.__sırala

    def __str__(self):
        return "  ".join(self.__sırala)

    def __add__(self, other):
        return self.__sırala + other.__sırala

    def prettify(self, squad="3 Group", kind="df", show=1):
        col = ["seq", "start", "end"]
        if squad == None:
            ndx = ["rstx"] * self.__sırala.__len__()
        else:
            size, name = squad.split(" ")
            n = 0
            ndx = []
            for base in self.__sırala:
                if len(ndx) % int(size) == 0:
                    n += 1
                ndx.append("{}_{}".format(str(name), n))
        if kind == "info":
            return list(set(ndx))
        if kind == "df":
            df = pd.DataFrame([N_tears().backmelting(i)[3] for i in self.__sırala],
                              columns=col,
                              index=ndx)
            n = 0
            _base = ["A", "T", "G", "C"]
            if show == 1:
                while n < 4:
                    son = [i.count(_base[n]) for i in list(df["seq"])]
                    df[_base[n]] = son
                    n += 1
            return df
        if kind == "list":
            return [N_tears().backmelting(i)[3] for i in self.__sırala]
        if kind == "str":
            if self.__expand == 1:
                return " <---> ".join(self.__sırala)
            else:
                return " ---> ".join(self.__sırala)

    def rates(self):
        info = {}
        adenin = self.__sequence.count("A")
        guanin = self.__sequence.count("G")
        sitozin = self.__sequence.count("C")
        timin = self.__sequence.count("T")
        # info["A/T"] = "{:.2f}".format(adenin/timin)
        # info["G/C"] = "{:.2f}".format(guanin/sitozin)
        info["A"] = adenin
        info["G"] = guanin
        info["C"] = sitozin
        info["T"] = timin
        return info


if __name__ == "__main__":
    print("test")
    print(Series(Sequence.random_seq(100), 10).prettify("1 seq"))
