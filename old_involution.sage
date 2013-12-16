from collections import deque


def get_alphabet(num_letters):
    """
    Returns an alphabet of specified length.

    If the length of the alphabet is at most 26, then
    letters of the English alphabet are used, otherwise
    numbers, starting from 1.

    INPUT:

    - ``num_letters`` - a positive integer

    OUTPUT:

    list -- the list containing the letters of the alphabet

    TESTS::

        sage: get_alphabet(3)
        ['a', 'b', 'c']
        sage: get_alphabet(50) == range(1,51)
        True
    """
    if num_letters <= 0:
        raise ValueError('The length of the alphabet must be '
                    'at least 1')
    if num_letters > 26:
        return range(1, num_letters + 1)
    else:
        return map(chr, range(97, 97 + num_letters))

class Involution(object):
    """
    Fixed point-free involution of a finite set.

    In general, the syntax for using this class is similar to
    the combinatorial permutation classes that are already
    part of Sage.

    INPUT:

    - ``pairing`` - a string or list of letters separated by spaces

    - ``flips`` - list of letters corresponding to flipped pairs,
      empty by default
      
    EXAMPLES::

        sage: Involution('a a b b c c')
        a a b b c c
        sage: Involution(['a','a','b','b','c','c'])
        a a b b c c
        sage: Involution('a a b b c c', flips = ['a', 'b'])
        -a -a -b -b c c

    Letters must come in pairs::

        sage: Involution('a a a b b a')
        Traceback (most recent call last):
        ...
        ValueError: Letters must come in pairs
        
    The flip list must contain valid letters::

        sage: Involution('a b c a b c', flips = ['a', 'd'])
        Traceback (most recent call last):
        ...
        TypeError: The flip list is not valid
        
    If a letter occurs in the flips more than once,
    it is still counted only once::

        sage: Involution('a e c c e a', flips = ['a', 'a'])
        -a e c c e -a
    """
    def __init__(self, pairing, flips = []):
        if isinstance(pairing, str):
            l = pairing.split()
        elif isinstance(pairing, list):
            l = pairing[:]
        letters = set(l)
        for letter in letters:
            if l.count(letter) != 2: 
                raise ValueError('Letters must come in pairs')

        flips = set(flips)
        
        for letter in flips:
            if letter not in letters:
                raise TypeError('The flip list is not valid')
        self._pairing = l
        self._flips = flips
        self._index = {}

        d = {}
        count = 0
        perm_list = [0] * len(l)
        for i in range(len(l)):
            if l[i] in d:
                perm_list[i] = d[l[i]] + 1
                perm_list[d[l[i]]] = i + 1
            else:
                d[l[i]] = i
                self._index[l[i]] = count
                count += 1
        self._permutation = Permutation(perm_list)

    def __repr__(self):
        def negate(x):
            if x in self._flips:
                return '-' + x
            else:
                return x
        l = map(negate, self._pairing)
        return ' '.join(l)

    def __getitem__(self, index):
        """
        Returns the letter at position 'index'.

        TESTS::

            sage: i = Involution('a a b b')
            sage: i[0]
            'a'
            sage: i[3]
            'b'
        """
        return self._pairing[index]

    def __eq__(self, other):
        """
        Tests equality up to renaming the letters.

        TESTS::
            
            sage: a = Involution('a a b b c c')
            sage: b = Involution('b b x x h h')
            sage: c = Involution('a b u u b a')
            sage: d = Involution('a a b b c c', flips = ['a'])
            sage: a == b
            True
            sage: a == c
            False
            sage: a == d
            False
        """
        if self._permutation != other._permutation:
            return False
        return map(self.index, self.flips()) == map(other.index,
                other.flips())

    @classmethod
    def _random(cls, num_letters, without_flips = True):
        """
        Constructs a random Involution.

        INPUT:

        - ``num_letters`` - The number of letters used for the
          Involution

        - ``without_flips`` - If True, the contructed Involution
          will not have flips, otherwise it is guaranteed to have
          flips. By default is is set to True.

        OUTPUT:

        - Involution - A random Involution

        TESTS::

            sage: i = Involution._random(3)
            sage: len(i.flips()) == 0
            True
            sage: i.alphabet() == {'a', 'b', 'c'}
            True
            sage: j = Involution._random(4, without_flips = False)
            sage: len(j.flips()) > 0
            True
        """
        import random
        l = 2 * get_alphabet(num_letters)
        random.shuffle(l)
        if without_flips:
            flips = []
        else:
            k = random.randint(1, num_letters)
            flips = random.sample(set(l), k)
        return Involution(l, flips)

    @classmethod
    def random_for_pseudo_anosov(cls, num_letters,
            without_flips = True):
        """
        Constructs a random involution arising from a pseudo-anosov.

        There are various reasons an involution does not arise from
        a pseudo-anosov on a closed surface:

        1. There is a one-pronged singularity

        2. It the involution yields a two-pronged singularity, it
          means that the same interval exchange can be specified by
          an involution on fewer intervals.

        3. There is a cyclic rotation of the letters where each 
          letter appears in the first and second half as well.
          This results in a saddle connection independently of the
          length parameters since we are now on the boundary of 
          a Mobius band so the starting point and the midpoint of
          the interval are connected.

        INPUT:

        - ``num_letters`` - The number of letters used for the
          Involution

        - ``without_flips`` - If True, the contructed Involution
          will not have flips, otherwise it is guaranteed to have
          flips. By default is is set to True.

        OUTPUT:

        - Involution - A random Involution arising form a PA

        TESTS:

        There is only one good involution without flips on 3 letters::

            sage: Involution.random_for_pseudo_anosov(3)
            a a b b c c

        A good involution with flips must be on at least 4 letters::
         
            sage: Involution.random_for_pseudo_anosov(3, 
            ....: without_flips = False)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        if num_letters <= 2:
            raise ValueError('A foliation without saddle connnection'
                    ' needs an involution with an alphabet of length'
                    ' at least three')
        if num_letters == 3:
            if without_flips:
                return Involution('a a b b c c')
            else:
                raise ValueError('A non-orientable pseudo-anosov '
                        'foliation must come from an involution '
                        'with alphabet of length at least four.')
        while True:
            inv = Involution._random(num_letters, without_flips)
            s = inv.singularity_type()
            if 1 not in s and 2 not in s and\
                    not inv.has_saddle_connection(): 
                return inv 

    @classmethod
    def arnoux_yoccoz_involution(cls, num_letters, 
            without_flips=True):
        """
        Returns the involution of the Arnoux-Yoccoz foliations.

        INPUT:

        - ``num_letters`` - number of letters in the alphabet

        - ``without_flips`` - if True, the the output is without
          flips, if False, all letters are flipped. The default is
          True

        OUTPUT:

        Involution - an Arnoux-Yoccoz involution

        EXAMPLES::

            sage: Involution.arnoux_yoccoz_involution(3)
            a a b b c c
            sage: Involution.arnoux_yoccoz_involution(4,
            ....: without_flips = False)
            -a -a -b -b -c -c -d -d
        """
        a = get_alphabet(num_letters)
        if without_flips:
            flips = []
        else:
            flips = a
        return Involution(sorted(2 * a), flips)

    def alphabet(self):
        """
        Returns the alphabet used for the encoding of the involution.

        OUTPUT:

        - set - the set of letters used for defining the involution

        EXAMPLES::

            sage: i = Involution('a a b c c b')
            sage: i.alphabet() == {'a', 'b', 'c'}
            True
            sage: j = Involution('b b c c a a', flips = ['a'])
            sage: j.alphabet() == {'a', 'b', 'c'}
            True
        """
        return set(self._pairing)

    def flips(self):
        """
        Returns the set of flips.

        OUTPUT:

        - set - the set of flips

        EXAMPLES::

            sage: i = Involution('a b a c c b', flips = ['b', 'c'])
            sage: i.flips() == set(['b', 'c'])
            True
        """
        return set(self._flips)

    def is_flipped(self, position):
        """
        Decides if the letter at specified position is flipped.

        INPUT:

        - ``position`` - a non-negative integer, less than twice
          the number of letters

        OUTPUT:

        - boolean - True if the letter at ``position`` is flipped, 
          otherwise False

        EXAMPLES::

            sage: i = Involution('a b c c b a', flips = ['a'])
            sage: [i.is_flipped(pos) for pos in range(6)]
            [True, False, False, False, False, True]
        """
        return self[position] in self._flips 

    def index(self, letter):
        """
        Returns the index of a letter in the alphabet.

        If there are n letters, there is unique bijection between the
        letters and range(n) - called the index - such that a letter
        with smaller index appears earlier it the encoding of the
        involution.

        INPUT:

        - ``letter`` - a letter of the alphabet of the involution

        OUTPUT:

        - non-negative integer - the index of the letter

        EXAMPLES::

            sage: i = Involution('a c c b b a')
            sage: i.index('a')
            0
            sage: i.index('b')
            2
            sage: i.index('c')
            1
        """
        return self._index[letter]

    def to_permutation(self):
        """
        Creates a Permutation from the Involution.

        OUTPUT:

        - Permutation - the permutation corresponding to the
          Involution

        EXAMPLES::

            sage: i = Involution('a a b b c c')
            sage: i.to_permutation()
            [2, 1, 4, 3, 6, 5]

        ::

            sage: j = Involution('a a b b c c', flips = ['a'])
            sage: i.to_permutation()
            [2, 1, 4, 3, 6, 5]

        """
        return self._permutation

    def list(self):
        """
        Returns the encoding of the involution as a list.

        OUTPUT:

        - list - the encoding list of the involution, ignoring the
          flips

        EXAMPLES::

            sage: i = Involution('a a b b c c', flips = ['a'])
            sage: i.list()
            ['a', 'a', 'b', 'b', 'c', 'c']
        """
        return self._pairing[:]

    def rotated(self, k):
        """
        Returns the Involution after cyclically rotating by k.

        INPUT:

        - ``k`` - any integer

        OUTPUT:

        - Involution - the resulting Involution after rotation

        EXAMPLES::

            sage: i = Involution('a a b b c c', flips = ['b']); i
            a a -b -b c c
            sage: i.rotated(4)
            -b -b c c a a
            sage: i.rotated(-1)
            a -b -b c c a
        """
        new_list = deque(self._pairing)
        new_list.rotate(k)
        return Involution(list(new_list), self._flips)

    def reversed(self):
        """
        Returns the Involution after reversing it.

        OUTPUT:

        - Involution - the reversed Involution

        EXAMPLES::

            sage: i = Involution('a a b b c c', flips = ['b']); i
            a a -b -b c c
            sage: i.reversed()
            c c -b -b a a
        """
        new_list = deque(self._pairing)
        new_list.reverse()
        return Involution(list(new_list), self._flips)

    def singularity_partition(self):
        """
        Returns the partitioning of separatrices in suspensions of
        the involution according to singularities.

        OUTPUT:

        - list of lists - partitioning of separatrices. The ordering
          of the lists is undefined, in particular they are not 
          necessarily sorted.

        EXAMPLES:

        The following involution identifies all separatrices::

            sage: i = Involution('a a b b c c')
            sage: set(i.singularity_partition()[0]) == {0,1,2,3,4,5}
            True

        The following involution yields foliations with 4
        singularities, one 3-pronged and three 1-pronged::

            sage: i = Involution('a a b b c c', flips = ['a','b','c'])
            sage: i.singularity_partition()     #not tested
            [[0, 2, 4], [1], [3], [5]]
        """
        done = set()
        partition = []
        n = len(self._pairing)
        for i in range(n):
            if i in done:
                continue
            j = i
            partition.append([])
            to_left = True
            while True:
                if to_left:
                    j = self._permutation[(j - 1) % n] - 1
                    if not self.is_flipped(j):
                        j = (j + 1) % n
                        to_left = False
                else:
                    j = self._permutation[j] - 1
                    if not self.is_flipped(j):
                        to_left = True
                    else:
                        j = (j + 1) % n

                partition[-1].append(j)
                done.add(j)
                if j == i:
                    break
        return partition

    def singularity_type(self):
        """
        Returns the singularity type of suspensions of this involution

        OUTPUT:

        - tuple - with number of prongs of singularities in
          decreasing order

        EXAMPLES::

            sage: i = Involution('a a b b c c')
            sage: i.singularity_type()
            (6,)

        ::

            sage: j = Involution('a a b b c c', flips = ['a','b','c'])
            sage: j.singularity_type()
            (3, 1, 1, 1)
        """
        return tuple(sorted(map(len, self.singularity_partition()), 
            reverse = True))

    def has_saddle_connection(self):
        """
        Returns True if any choice of length parameters results in a
        foliation with saddle connection.

        This happend if there is a cyclic rotation of the letters 
        where each 
        letter appears in the first and second half as well.
        This results in a saddle connection independently of the
        length parameters since we are now on the boundary of 
        a Mobius band so the starting point and the midpoint of
        the interval are connected.

        OUTPUT:

        - boolean - True if any choice of length parameters results
          in a saddle connection, otherwise False

        EXAMPLES::

            sage: i = Involution('a a b b c c')
            sage: i.has_saddle_connection()
            False

        ::

            sage: i = Involution('a a b b')
            sage: i.has_saddle_connection()
            True

        ::

            sage: i = Involution('a b c a d b c d', flips = ['d'])
            sage: i.has_saddle_connection()
            True
        """
        n = len(self._pairing) / 2
        for i in range(n):
            if len(set(self._pairing[i:i+n])) == n:
                return True
        return False









class TwoSidedInvolution(SageObject):
    def __init__(self, top_letters, bottom_letters, flips = []):
        self._gen_perm = iet.GeneralizedPermutation(\
                top_letters, bottom_letters, flips = flips)

    def __repr__(self):
        return repr(self._gen_perm)

    def __getitem__(self, index):
        return self._gen_perm[index]

    def __eq__(self, other):
        return self._gen_perm == other._gen_perm

    @classmethod
    def _random(cls, top_length, bottom_length, without_flips = True):
#        if top_length <= 0 or bottom_length <= 0:
#            raise ValueError('The number of intervals on top and'
#                    ' bottom must be positive')
        if (top_length + bottom_length) % 2 != 0:
            raise ValueError('There should be even number of'
                    ' intervals in total')
        i = Involution.random((top_length + bottom_length) / 2, 
                without_flips)
        l = i.list()
        return TwoSidedInvolution(l[:top_length], l[top_length:],
                i.flips())


