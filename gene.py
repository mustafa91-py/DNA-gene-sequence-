from seq import Sequence, Series, Protein, N_tears


class Distance(object):
    """
    Class Distance :
    purporse : Measure and map all possible distances between two motifs within the DNA sequence
    """
    def __init__(self, seq, motif1=None, motif2=None, defect=True):
        """

        :param seq: sequence
        :param motif1:first motif
        :param motif2: second motif
        :param defect: condition of defects that may occur
        """
        self.seq = seq
        self.__mtf1, self.__mtf2 = motif1, motif2
        self.__dfct = defect
        if self.__mtf1 is not None and self.__mtf2 is not None:
            self.__data__ = Sequence.distance_seq(self.seq, self.__mtf1, self.__mtf2, "prettify", self.__dfct)

    def srange(self, mini: int = -10, maxi: int = 10) -> list:
        """

        :param mini:filter data min range
        :param maxi: filter data max range
        :return: list
        """
        result = list()
        for i in self.__data__:
            name = [i[0]]
            range_ = filter(lambda x: mini <= x <= maxi, i[1:])
            name.extend(list(range_))
            if len(name) != 1:
                result.append(name)
        return result

    def nearest_distance(self, exp=0):
        """
        purporse: Returns the closest distance between filtered or unfiltered data
        :param exp:Description
        :return:list
        """
        main_motif, side_motif = self.__mtf1, self.__mtf2
        results = []
        out = {}
        put = {}
        DATA = self.__data__
        try:
            if len(DATA[0]) < 2:
                return None
        except IndexError:
            return "NOT FOUND"
        mirror = False
        for base in DATA:
            try:
                minimum = min(map(abs, base[1:]))
                out[int(minimum)] = base[0]
                put[base[0]] = int(minimum)
            except TypeError:
                raise TypeError(f'bad parameters : {main_motif} or {side_motif}')
        min_ = min(out.keys())
        for base in DATA:
            if min_ == put.get(base[0]):
                negative = []
                positive = []
                query = None
                _x = None
                _y = None
                for x in base[1:]:
                    if x < 0:
                        negative.append(x)
                    if x >= 0:
                        positive.append(x)
                if len(negative) != 0:
                    _x = max(negative)
                    query = _x
                if len(positive) != 0:
                    _y = min(positive)
                    query = _y
                if len(negative) != 0 and len(positive) != 0:
                    if abs(_x) < abs(_y):
                        query = _x
                    else:
                        query = _y
                if len(positive) != 0 and len(negative) != 0:
                    if abs(max(negative)) == abs(min(positive)):
                        query = 0
                        mirror = True
                if query < 0:
                    d = f'explanation : {min_} units of "{side_motif}" to the right of "{base[0]}"'
                    if exp:
                        v = (base[0], "right", side_motif, min_, d)
                        results.append(v)
                    if not exp:
                        v = (base[0], "right", side_motif, min_)
                        results.append(v)
                if query > 0:
                    d = f'explanation : {min_} units of "{side_motif}" to the left of "{base[0]}"'
                    if exp:
                        v = (base[0], "left", side_motif, min_, d)
                        results.append(v)
                    if not exp:
                        v = (base[0], "left", side_motif, min_)
                        results.append(v)
                if query == 0:
                    if mirror:
                        disc = "mirror"
                    else:
                        disc = "adjoining"
                    d = f'explanation : {min_} units of "{side_motif}" to the {disc} of "{base[0]}"'
                    if exp:
                        v = (base[0], disc, side_motif, min_, d)
                        results.append(v)
                    if not exp:
                        v = (base[0], disc, side_motif, min_)
                        results.append(v)
        return results

    def get(self, iterable, specifity=False):
        """

        :param iterable:list(-nearest_distance-)
        :param specifity: id
        :return: list
        """
        put = []
        if isinstance(iterable, list):
            for base in iterable:
                main_seq, x1, x2 = Sequence.backmelting(base=base[0])[3]
                x1, x2 = int(x1), int(x2)
                sea_lane = base[1]
                side_seq = base[2]
                step_lie = base[3]
                if sea_lane == "right":
                    terminal_1 = x1
                    terminal_2 = x2 + step_lie + len(side_seq)
                    definition = self.seq[terminal_1:terminal_2]
                    if not specifity:
                        put.append(definition)
                    if specifity:
                        put.append(f'{terminal_1} <' + definition + f'<= {terminal_2}')
                if sea_lane == "left":
                    terminal_1 = x1 - step_lie - len(side_seq)
                    terminal_2 = x2
                    definition = self.seq[terminal_1:terminal_2]
                    if x1 == 0:
                        definition = self.seq[0:terminal_2]
                    if not specifity:
                        put.append(definition)
                    if specifity:
                        put.append(f'{terminal_1} <' + definition + f'<= {terminal_2}')
                if sea_lane == "adjoining":
                    terminal_1 = x1
                    terminal_2 = x2
                    definition_1 = self.seq[terminal_1 - len(side_seq):terminal_2]
                    definition_2 = self.seq[terminal_1:terminal_2 + len(side_seq)]
                    if definition_1 == main_seq + side_seq or definition_1 == side_seq + main_seq:
                        if not specifity:
                            put.append(definition_1)
                        if specifity:
                            put.append(f'{terminal_1 - len(side_seq)} <' + definition_1 + f'<= {terminal_2}')
                    if definition_2 == main_seq + side_seq or definition_2 == side_seq + main_seq:
                        if not specifity:
                            put.append(definition_2)
                        if specifity:
                            put.append(f'{terminal_1} <' + definition_2 + f'<= {terminal_2 + len(side_seq)}')
                if sea_lane == "mirror":
                    terminal_1 = x1 - len(side_seq) - step_lie
                    terminal_2 = x2 + len(side_seq) + step_lie
                    definition = self.seq[terminal_1:terminal_2]
                    if not specifity:
                        put.append(definition)
                    if specifity:
                        put.append(f'{terminal_1} <' + definition + f'<= {terminal_2}')
        return put


class Test():
    def __init__(self, seq):
        self.seq = seq
        self.dfct = True
        self.__control = 0
        self.__y = []
        self.__x = []
        self.__valuer_turn = 1

    @property
    def defect(self):
        self.dfct = None

    @defect.setter
    def defect(self, bool):
        self.dfct = bool

    def promotor(self, motif):
        out = list()
        for base in list(Sequence.find_motif(sequence=self.seq, motif=motif, rate=99, specifity=True)):
            out.append(Series.backmelting(base=base[0])[3])
        return out

    def fibonacci_seq(self, motif1, motif2):
        """
        testing...
        purporse : creates the fibonacci algorithm between two sequences
        :param motif1:first motif
        :param motif2:second motif
        :return:
        """
        a = Distance(seq=self.seq, motif1=motif1, motif2=motif2).nearest_distance()
        b = Distance(seq=self.seq, motif1=motif1, motif2=motif2).get(a)
        self.__y.append(a)
        self.__control += 1
        if self.__control > len(self.seq) / 50:
            return None
        for base in b:
            return self.fibonacci_seq(motif1=base, motif2=motif1)

    @property
    def fiboncci(self):
        for i in self.__y:
            self.__x.extend(Distance(seq=self.seq).get(i, False))
        ddd = sorted(list(set(self.__x)), key=lambda x: len(x), reverse=True)
        self.__x = []
        self.__y = []
        if self.__valuer_turn == 0:
            try:
                return ddd[0].__str__()
            except IndexError:
                return None
        else:
            return ddd[:self.__valuer_turn]

    @fiboncci.setter
    def value(self, value=1):
        self.__valuer_turn = value

    def ekok(self, motif1, motif2):
        """
        testing...
        :param motif1:
        :param motif2:
        :return:
        """
        self.value = 0
        self.fibonacci_seq(motif1=motif1, motif2=motif2)
        a = self.fiboncci
        self.fibonacci_seq(motif1=motif2, motif2=motif1)
        b = self.fiboncci
        self.value = 1
        return Sequence.common_array_(a, b, "duo")


class Gen_control(Test):
    """
    testing...
    """
    value: int
    seq = str

    def __init__(self, seq):
        super(Gen_control, self).__init__(seq)


if __name__ == "__main__":
    gen = Gen_control(Sequence.random_seq(100))
    gen.value = 20
    fib = gen.fibonacci_seq("TA", "CC")
    print(gen.fiboncci)
