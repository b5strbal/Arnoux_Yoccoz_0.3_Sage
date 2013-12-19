from collections import deque


class Involution(SageObject):
    """
    A wrapper for different kinds of Generalized Permutations.

    Foliation objects are based on Involutions which are 
    basically the same as GeneralizedPermutations already 
    in sage. Involution objects have a few extra methods.
   
    INPUT:

    - ``top_letters`` - a string where interval names are
      separated by spaces, or a list is interval names
    - ``bottom_letters`` - a string where interval names are
      separated by spaces, or a list is interval names

    - ``flips`` - a list or set of flipped interval names, or
      a string containing flipped letters in case the names 
      are single characters

    EXAMPLES:

    The top and bottom intervals can be specified as strings,
    the letters (or numbers or words) separated by spaces. The
    flips can be a string in containing the letters of flipped
    intervals::

        sage: Involution('a b c', 'c b a', flips = 'ab')
        -a -b  c
         c -b -a

    If the names of intervals are not single characters, this
    flips notation doesn't work::

        sage: Involution('1 2 33', '33 2 1', flips = '33')
        Traceback (most recent call last):
        ...
        TypeError: The flip list is not valid

    In this case the flipped intervals must be presented in
    a list or set so that element testing works::

        sage: Involution('1 2 33', '33 2 1', flips = ['33'])
         1  2 -33
        -33  2  1
        
    The top and bottom intervals can also be listed in a list::

        sage: Involution(['a','b','c'],['c','b','a'], 'ab')
        -a -b  c
         c -b -a

    If the second argument is omitted, the bottom side of the
    curve is considered as a Moebius band without punctures::

        sage: Involution('a a b b c c', flips = 'abc')
        -a -a -b -b -c -c
        Moebius band
        sage: _.singularity_type()
        (3, 1, 1, 1)

    It is only be omitting the bottom letters that one gets
    a Moebius band without punctures. The following results in
    a once punctured Moebius band on the bottom::

        sage: Involution('a a b b c c', 'd d', flips = 'abc')
        -a -a -b -b -c -c
         d  d
        sage: _.singularity_type()
        (3, 2, 1, 1, 1)

    If both arguments are omitted, then the trivial permutation
    is constructed whose suspension is a torus without
    punctures::

        sage: Involution()
        Torus
        sage: _.singularity_type()
        ()

    Again, this is the only way of specifying an unpunctured
    torus. The following constructs a once-punctured torus::

        sage: Involution('a', 'a')
        a
        a
        sage: _.singularity_type()
        (2,)

    And here is a twice-punctured torus::

        sage: Involution('a b', 'a b')
        a b 
        a b
        sage: _.singularity_type()
        (2, 2)

    """
    def __init__(self, top_letters = None, 
            bottom_letters = None, 
            flips = []):
        if top_letters == None: # torus
            top_letters = bottom_letters = 'JOKER'
        if bottom_letters == None: # bottom side is Moebius
            bottom_letters = 'JOKER JOKER'
        self._gen_perm = iet.GeneralizedPermutation(\
                top_letters, bottom_letters, flips = flips)

        # initializing self._pair
        self._pair = {}
        for i in range(2):
            for j in range(len(self[i])):
                for a in range(2):
                    for b in range(len(self[a])):
                        if (i != a or j != b) and \
                                self[i][j] == self[a][b]:
                                    self._pair[(i,j)] = (a,b)

        # initializing self._index
        self._index = {}
        count = 0
        done = set()
        for (i, j) in sorted(self._pair.keys()):
            letter = self[i][j]
            if letter in done:
                continue
            done.add(letter)
            self._index[letter] = count
            count += 1

        # initializing self._singularity_partition
        done = set()
        partition = []
        for (i, j) in self._pair:
            if (i, j) in done:
                continue
            (a, b) = (i, j)
            partition.append([])
            direction = 'left'
            while True:
                if direction == 'left':
                    (a, b) = self._pair[(a, (b - 1) % 
                        len(self[a]))]
                    if not self.is_flipped((a, b)):
                        b = (b + 1) % len(self[a])
                        direction = 'away'
                else:
                    (a, b) = self._pair[(a,b)]
                    if not self.is_flipped((a, b)):
                        direction = 'left'
                    else:
                        b = (b + 1) % len(self[a])
                partition[-1].append((a, b))
                done.add((a, b))
                if (a, b) == (i, j):
                    break
        self._singularity_partition = partition


    def __repr__(self):
        """
        Returns a representation of self.

        EXAMPLES:

        Usually it is the same as the representation of a
        GeneralizedPermutation::

            sage: Involution('a a b b', 'c c', flips = 'ab')
            -a -a -b -b
             c  c

        It's different when the bottom side is a Moebius band::

            sage: Involution('a a b b', flips ='ab')
            -a -a -b -b
            Moebius band

        Or when the suspension is a torus::

            sage: Involution()
            Torus

        """
        if self.is_torus():
            return 'Torus'
        if self.is_bottom_side_moebius():
            return repr(self._gen_perm).split('\n')[0] + \
                    '\nMoebius band'
        return repr(self._gen_perm)

    def __getitem__(self, index):
        """
        Returns the list of top of bottom letters.

        INPUT:

        - ``index`` - 0 for the top letters, 1 for the bottom

        EXAMPLES::

            sage: i = Involution('a a b b','c c',flips = 'b');i
            a a -b -b
            c c
            sage: i[0]
            ['a', 'a', 'b', 'b']
            sage: i[1]
            ['c', 'c']

        """
        return self._gen_perm.list()[index]

    def __eq__(self, other):
        """
        Decides if two Involutions are the same up to renaming
        letters.

        EXAMPLES::
            
            sage: i = Involution('a a b b', 'c c')
            sage: j = Involution('1 1 2 2', '3 3')
            sage: k = Involution('a a b b', 'c c', flips = 'c')
            sage: l = Involution('a a b b')
            sage: i == j
            True
            sage: i == k
            False
            sage: i == l
            False
        """ 
        if self.is_torus():
            return other.is_torus()
        if self.is_bottom_side_moebius() !=\
                other.is_bottom_side_moebius():
                    return False
        return self._gen_perm == other._gen_perm

    def alphabet(self):
        """
        Returns the set of interval names.

        OUTPUT:

        set -- the set of interval names

        EXAMPLES::

            sage: i = Involution('a a b b','c c',flips='c')
            sage: i.alphabet() == {'a', 'b', 'c'}
            True

            sage: i = Involution()
            sage: i.alphabet()
            set([])
        """
        s = set(self._gen_perm.alphabet())
        s.discard('JOKER')
        return s

    def flips(self):
        """
        Returns the list of flips.

        OUTPUT:

        - list -- the list of flipped interval names

        EXMAPLES::

            sage: i = Involution('a a b b','c c',flips = 'ab')
            sage: i.flips()
            ['a', 'b']

            sage: i = Involution('a a', flips = 'a')
            sage: i.flips()
            ['a']

            sage: i = Involution()
            sage: i.flips()
            []

        """
        if isinstance(self._gen_perm[0][0], tuple):
            # self._gen_perm is a FlippedLabelledPermutationLI
            # that has a flips() method
            return self._gen_perm.flips()
        else: #self._gen_perm is a LabelledPermutationIET 
            # hence there are no flips
            return []

    def is_flipped(self, pos):
        """
        Decides if the interval at a certain position is 
        flipped.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.

        OUTPUT:

        - boolean -- True is the interval is flipped, False
          is not

        EXAMPLES::

            sage: i = Involution('a a b b','c c', flips='bc');i
            a a -b -b
            -c -c
            sage: i.is_flipped((0,0))
            False
            sage: i.is_flipped((0,1))
            False
            sage: i.is_flipped((0,2))
            True
            sage: i.is_flipped((0,3))
            True
            sage: i.is_flipped((1,0))
            True
            sage: i.is_flipped((1,1))
            True

        """
        x = self._gen_perm[pos[0]][pos[1]]
        if isinstance(x, tuple):
            # self._gen_perm is a FlippedLabelledPermutationLI
            return x[1] == -1
        #self._gen_perm is a LabelledPermutationIET 
        return False

    def index(self, letter):
        """
        Returns the index of an letter.

        If n letters are used to notate the involution, they
        are indexed from 0 to n-1 according to their first 
        occurrence when read from to to bottom, from left to
        right.

        INPUT:

        - ``letter`` - string

        OUTPUT:

        - integer - the index of the letter

        EXAMPLES::

            sage: i = Involution('a b a c','d c d b')
            sage: i.index('a')
            0
            sage: i.index('b')
            1
            sage: i.index('c')
            2
            sage: i.index('d')
            3

        """
        return self._index[letter]

    def pair(self, pos):
        """
        Returns the position of the pair of the interval at
        a specified position.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.

        OUTPUT:

        - tuple -- the position of the pair

        EXAMPLES::

            sage: i = Involution('a b a b','c c', flips = 'ab')
            sage: i.pair((0,0))
            (0, 2)
            sage: i.pair((0,2))
            (0, 0)
            sage: i.pair((1,1))
            (1, 0)

        """
        return self._pair[pos]

    @classmethod
    def orientable_arnoux_yoccoz(self, genus):
        """
        Returns the Involution of the Arnoux-Yoccoz foliations
        on orientable surfaces.

        INPUT:

        - ``genus`` - the genus of the surface

        OUTPUT:

        - Involution

        EXAMPLES::

            sage: Involution.orientable_arnoux_yoccoz(3)
            1 2 3 4 5 6
            2 1 4 3 6 5
            sage: Involution.orientable_arnoux_yoccoz(4)
            1 2 3 4 5 6 7 8
            2 1 4 3 6 5 8 7

        """
        if genus < 3:
            raise ValueError('The genus of an orientable'
                    'Arnoux-Yoccoz surface is at least 3')
        n = 2 * genus
        top = range(1, n + 1)
        def switch(k):
            if k % 2 == 0:
                return k + 1
            return k - 1
        bottom = [top[switch(i)] for i in range(n)]
        return Involution(top, bottom)

    @classmethod
    def nonorientable_arnoux_yoccoz(self, genus):
        """
        Returns the Involution of the Arnoux-Yoccoz foliations
        on non-orientable surfaces.

        Take care with the the genus here: 
        The orienting double cover of a closed genus 
        $g$ non-orientable surface is the closed genus $g-1$
        orientable surface.
        
        INPUT:

        - ``genus`` - the non-orientable genus of the surface

        OUTPUT:

        - Involution -

        EXAMPLES::

            sage: Involution.nonorientable_arnoux_yoccoz(4)
            1 1 2 2 3 3
            Moebius band
            sage: Involution.nonorientable_arnoux_yoccoz(5)
            1 1 2 2 3 3 4 4
            Moebius band

        """
        if genus < 4:
            raise ValueError('The genus of a non-orientable '
                    'Arnoux-Yoccoz surface is at least 4')
        top = sorted(2 * range(1, genus))
        return Involution(top)

    @classmethod
    def RP2_arnoux_yoccoz(self):
        """
        Returns the Involution of the Arnoux-Yoccoz foliation
        on the projective plane.

        OUTPUT:

        - Involution -

        EXAMPLES::

            sage: Involution.RP2_arnoux_yoccoz()
            -a -a -b -b -c -c
            Moebius band

        """
        return Involution('a a b b c c', flips = 'abc')

    def _top_deque(self):
        """
        Returns the list of top letters as a deque.

        OUTPUT:

        - deque -- the deque of top letters

        TESTS::

            sage: i = Involution('a a b b','c c', flips='ab')
            sage: i._top_deque()
            deque(['a', 'a', 'b', 'b'])

            sage: i = Involution('a a b b','c c')
            sage: i._top_deque()
            deque(['a', 'a', 'b', 'b'])
        """
        return deque(self._gen_perm.list()[0])

    def _bottom_deque(self):
        """
        Returns the list of bottom letters as a deque.

        OUTPUT:

        - deque -- the deque of bottom letters

        TESTS::

            sage: i = Involution('a a b b','c c', flips='ab')
            sage: i._bottom_deque()
            deque(['c', 'c'])

            sage: i = Involution('a a b b','c c')
            sage: i._bottom_deque()
            deque(['c', 'c'])

        """
        return deque(self._gen_perm.list()[1])

    def rotated(self, top, bottom):
        """
        Returns an involution where the top and bottom rows
        are rotated cyclically.

        INPUT:

        - ``top`` - an integer, shift the top letters
          cyclically by this amount

        - ``bottom`` - an integer, shift the bottom letters
          cyclically by this amount

        OUTPUT:

        - Involution - the rotated Involution

        EXAMPLES::

            sage: i = Involution('a b c b','c a d d', \
                    flips='bc');i
            a -b -c -b
            -c a d d
            sage: i.rotated(1, 1)
            -b a -b -c
            d -c a d
            sage: i.rotated(-6, 2)
            -c -b a -b
            d d -c a
            
        """
        top_list = self._top_deque() 
        bottom_list = self._bottom_deque() 
        top_list.rotate(top)
        bottom_list.rotate(bottom)
        return Involution(list(top_list), list(bottom_list),
                self.flips())


    def reversed(self):
        """
        Returns an involution where the top and bottom rows
        are reversed.

        OUTPUT:

        - Involution - the reversed Involution

        EXAMPLES::

            sage: i = Involution('a b c b','c a d d', \
                    flips='bc');i
            a -b -c -b
            -c a d d
            sage: i.reversed()
            -b -c -b a
            d d a -c

        """
        top_list = self._top_deque() 
        bottom_list = self._bottom_deque() 
        top_list.reverse()
        bottom_list.reverse()
        return Involution(list(top_list), list(bottom_list),
                self.flips())
        
    def singularity_type(self):
        """
        Returns the singularity type of self.

        The suspension of the Involution yields a foliation.
        The singularity type of that foliation is the tuple of
        the number of prongs at singularities.

        OUTPUT:

        - tuple -

        EXAMPLES::

        sage: Involution('a a b b c c', flips = 'abc')
        -a -a -b -b -c -c
        Moebius band
        sage: _.singularity_type()
        (3, 1, 1, 1)

        sage: Involution('a a b b c c', 'd d', flips = 'abc')
        -a -a -b -b -c -c
         d  d
        sage: _.singularity_type()
        (3, 2, 1, 1, 1)

        sage: Involution()
        Torus
        sage: _.singularity_type()
        ()

        sage: Involution('a', 'a')
        a
        a
        sage: _.singularity_type()
        (2,)

        sage: Involution('a b', 'a b')
        a b 
        a b
        sage: _.singularity_type()
        (2, 2)

        """
        if self.is_torus():
            return ()
        t = sorted([x for x in map(len, 
            self._singularity_partition)], reverse = True)
        if self.is_bottom_side_moebius():
            t.remove(2)
        return tuple(t)

    def is_torus(self):
        """
        Decides if the suspension of self is a torus without
        punctures.

        OUTPUT:

        - boolean

        EXAMPLES:

        This is the only way to get a torus::

            sage: i = Involution()
            sage: i.is_torus()
            True

        Even tori with punctures return False::

            sage: i = Involution('a', 'a')
            sage: i.is_torus()
            False

        """
        return self[0] == ['JOKER']

    def is_bottom_side_moebius(self):
        """
        Decides if the bottom side is Moebius band without
        punctures.

        This happends exactly when the second argument in
        the constructor is omitted.

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: i = Involution('a a b b')
            sage: i.is_bottom_side_moebius()
            True

        Here the bottom component is a once punctured Moebius
        band, so it returns False::

            sage: i = Involution('a a b b', 'c c')
            sage: i.is_torus()
            False

        """
        return self[1] == ['JOKER', 'JOKER']

    def is_surface_orientable(self):
        """
        Decides if the suspension surface is orientable.

        OUTPUT:

        - boolean -- 

        EXAMPLES:

        If the bottom side is a Moebius band, it is always
        non-orientable::

            sage: i = Involution('a a b b')
            sage: i.is_surface_orientable()
            False

        Or if there is flipped pair on different sides::

            sage: i = Involution('a b c','c b a', flips = 'a')
            sage: i.is_surface_orientable()
            False

        Or if there is not flipped pair on the same side::

            sage: i = Involution('a a b b', 'c c', flips='ac')
            sage: i.is_surface_orientable()
            False

        Otherwise it is orientable::

            sage: i = Involution('a b c', 'b c a')
            sage: i.is_surface_orientable()
            True

        """ 
        for x in self._pair:
            a = (x[0] == self._pair[x][0])
            b = self.is_flipped(x)
            if a != b:
                return False
        return True

    def is_foliation_orientable(self):
        """
        Decides if the suspension foliations are orientable.

        OUTPUT:

        - boolean --

        EXAMPLES:

        If there are flipped intervals, the foliation is 
        non-orientable::

            sage: i = Involution('a a b b', 'c c', flips ='a')
            sage: i.is_foliation_orientable()
            False

        If there are no flips, the foliation is orientable::

            sage: i = Involution('a a b b', 'c c')
            sage: i.is_foliation_orientable()
            True

        """
        return len(self.flips()) == 0
















#from involution import Involution
from bisect import bisect_left, bisect
from collections import namedtuple

epsilon = 1e-10

class PointWithCoefficients(SageObject):
    """
    Represents a real number together with its coefficients as
    a linear combination of the a certain set of real numbers.

    In practice this set is the set of interval lengths and 
    possibly the twist of a Foliation.

    INPUT:

    - ``value`` - a real number (can be floating point, 
      or an exact algebraic number)

    - ``coefficients`` - a list or tuple or vector of 
      the coefficients

    EXAMPLES:

    The coefficients can be specified as a list, tuple or 
    vector::

        sage: a = PointWithCoefficients(1.2, (3, -1, 0))
        sage: b = PointWithCoefficients(0.8, [2, 1, 5])
        sage: c = PointWithCoefficients(3.4, vector((2, 3)))

    One can add or subtract two objects as long as they are 
    of the same type:

        sage: a + b
        (2.00000000000000, (5, 0, 5))
        sage: a - b
        (0.400000000000000, (1, -2, -5))
        sage: a - c
        Traceback (most recent call last):
        ...
        TypeError: Cannot subtract two PointWithCoefficients of different types
    
    """
    def __init__(self, value, coefficients):
        self.value = value
        self.coefficients = vector(coefficients)

    def __repr__(self):
        """
        Returns the representation of self.

        TESTS::

            sage: PointWithCoefficients(3, (4, 3, 2))
            (3, (4, 3, 2))

        """
        return repr((self.value, self.coefficients))

    def __add__(self, other):
        """
        Adds the numbers and their coefficient vectors.

        OUTPUT:

        - PointWithCoefficients - the sum

        TESTS::

        sage: a = PointWithCoefficients(1.2, (3, -1, 0))
        sage: b = PointWithCoefficients(0.8, (2, 1, 5))
        sage: c = PointWithCoefficients(3.4, (2, 3))

        sage: a + b
        (2.00000000000000, (5, 0, 5))

        sage: a + c
        Traceback (most recent call last):
        ...
        TypeError: Cannot add two PointWithCoefficients of different types

        """
        try:
            return PointWithCoefficients(self.value + 
                    other.value,
                    self.coefficients + other.coefficients)
        except TypeError:
            raise TypeError("Cannot add two "
            "PointWithCoefficients of different types")
        
    def __sub__(self, other):
        """
        Subtracts the numbers and their coefficient vectors.

        OUTPUT:

        - PointWithCoefficients - the difference

        TESTS::

        sage: a = PointWithCoefficients(1.2, (3, -1, 0))
        sage: b = PointWithCoefficients(0.8, [2, 1, 5])
        sage: c = PointWithCoefficients(3.4, (2, 3))

        sage: a - b
        (0.400000000000000, (1, -2, -5))
        sage: a - c
        Traceback (most recent call last):
        ...
        TypeError: Cannot subtract two PointWithCoefficients of different types

        """
        try:
            return PointWithCoefficients(self.value - 
                    other.value,
                    self.coefficients - other.coefficients)
        except TypeError:
            raise TypeError("Cannot subtract two "
                    "PointWithCoefficients of different types")

    def mod_one(self):
        """
        Returns the PointWithCoefficients corresponding to the
        real number of self modulo 1.

        The sum of the numbers in the generating set used for
        linear combinations is 1 hence all but the twist
        coefficient is decreased by the floor of the real
        number of self.

        OUTPUT:

        - PointWithCoefficients --

        EXAMPLES::

            sage: p = PointWithCoefficients(3.2, (4, 5, 3))
            sage: p.mod_one()
            (0.200000000000000, (1, 2, 3))

            sage: q = PointWithCoefficients(-1.3, (0, 2, 4,-2))
            sage: q.mod_one()
            (0.700000000000000, (2, 4, 6, -2))

        """
        if self.value == 0: #need to check beacuse of a bug
            # with algebraic numbers
            n = 0
        else:
            n = floor(self.value)
        return PointWithCoefficients(self.value - n, [x - n 
            for x in self.coefficients[:-1]] + \
                    [self.coefficients[-1]])

def arnoux_yoccoz_factor(genus, field = RDF):
    """
    Returns the Perron root of the polynomials 
    $x^g - x^{g-1} - ... - 1$.

    INPUT:

    - ``genus`` - integer, the genus of the orientable
      closed Arnoux-Yoccoz surface, which is the same as the degree
      of the minimalpolynomial.

    - ``field`` - the number field in which roots of the polynomial
      are searched. The default is RDF, but one want the exact number
      in QQbar.

    OUTPUT:

    - number -- number type depends on the field specified

    EXAMPLES::

        sage: arnoux_yoccoz_factor(3)
        1.83928675521
        sage: arnoux_yoccoz_factor(3, field = QQbar)
        1.839286755214161?

    """
    R.<x> = ZZ[]
    poly = R([-1] * genus + [1])
    return max([abs(x[0]) for x in poly.roots(field)])

class Interval(object):
    """
    An interval of the unit interval $[0,1]$ with opposite sides 
    identified.

    INPUT:

    - ``left_endpoint`` - PointWithCoefficients
    - ``right_endpoint`` - PointWithCoefficients
    - ``left_openness`` - 'closed' or 'open'. The default is 'closed'
    - ``right_openness`` - 'closed' or 'open'. The default is 'open'

    EXAMPLES:

    Any usual interval can be specified::

        sage: Interval(PointWithCoefficients(1/2, [2, 1]), \
                PointWithCoefficients(3/5, [4, -1]))
        [(1/2, (2, 1)), (3/5, (4, -1)))
        sage: Interval(PointWithCoefficients(0, [0]), \
                PointWithCoefficients(0.4, [1]), 'open', 'closed')
        ((0, (0)), (0.400000000000000, (1))]
        
    But intervals wrapping around the endpoints are also considered::

        sage: Interval(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        [(1/2, (0, 1)), (1/4, (1, -2))]

    One can get the endpoints by indexing and the length() using the 
    length()::

        sage: i = Interval(PointWithCoefficients(1/5, (1, 2)), \
                PointWithCoefficients(4/7, (6, -1)))
        sage: i[0]
        (1/5, (1, 2))
        sage: i[1]
        (4/7, (6, -1))
        sage: i.length()
        (13/35, (5, -3))

        sage: j = Interval(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        sage: j.length()
        (3/4, (2, -3))
    
    """
    def __init__(self, left_endpoint, right_endpoint,
            left_openness = 'closed', right_openness = 'open'):
        if not isinstance(left_endpoint, PointWithCoefficients) or \
                not isinstance(right_endpoint,PointWithCoefficients) :
            raise ValueError("The endpoints of the interval must be "
                    "PointWithCoefficients objects.")
        if not 0 <= left_endpoint.value < 1 or not \
                0 <= right_endpoint.value < 1:
                    raise ValueError("The endpoints of the Interval "
                            "must be between 0 and 1.")
        self._lep = left_endpoint
        self._rep = right_endpoint
        if not {'closed', 'open'} >= {left_openness, right_openness}:
            raise ValueError('Openness arguments should be either '
                    '\'closed\' or \'open\'')
        self._left_openness = left_openness
        self._right_openness = right_openness

    @staticmethod
    def _less(x, y, is_equality_allowed):
        """
        Unifies the < and <= operators in a single function.

        INPUT:

        - ``x`` - a number
        - ``y`` - another number
        - ``is_equality_allowed`` - if True, the function returns 
          x<=y, otherwise x<y

        TESTS::

            sage: Interval._less(1, 2, True)
            True
            sage: Interval._less(2, 1, True)
            False
            sage: Interval._less(1, 1, True)
            True
            sage: Interval._less(1, 1, False)
            False

        """
        if is_equality_allowed:
            return x <= y
        else:
            return x < y

    def __repr__(self):
        """
        Returns the representation of self.

        TESTS::

        sage: Interval(PointWithCoefficients(1/2, [0, 1]), \
                PointWithCoefficients(1/4, [1, -2]),'closed','closed')
        [(1/2, (0, 1)), (1/4, (1, -2))]

        """
        s = ''
        if self._left_openness == 'closed':
            s += '['
        else:
            s += '('
        s += str(self._lep) + ', ' + str(self._rep)
        if self._right_openness == 'closed':
            s += ']'
        else:
            s += ')'
        return s

    def __getitem__(self, index):
        """
        Returns the left or right endpoint of the interval.

        INPUT:

        - ``index`` - either 0 or 1, for the left and right endpoints
          respectively

        OUTPUT:

        - PointWithCoefficients - one of the endpoints

        EXAMPLES::

            sage: i = Interval(PointWithCoefficients(1/5, (1, 2)), \
                    PointWithCoefficients(4/7, (6, -1)))
            sage: i[0]
            (1/5, (1, 2))
            sage: i[1]
            (4/7, (6, -1))

        """
        if index == 0:
            return self._lep
        if index == 1:
            return self._rep

    def length(self):
        """
        Returns the length of the interval.

        OUTPUT:

        - PointWithCoefficients - the length of the interval

        EXMAPLES::

            sage: i = Interval(PointWithCoefficients(1/5, (1, 2)), \
                    PointWithCoefficients(4/7, (6, -1)))
            sage: i.length()
            (13/35, (5, -3))

            sage: j = Interval(PointWithCoefficients(0.5, [0, 1]), \
                PointWithCoefficients(0.25,[1, -2]),'closed','closed')
            sage: j.length()
            (0.750000000000000, (2, -3))
    
        """
        return (self._rep - self._lep).mod_one()

    def contains(self, point):
        """
        Decides if a point in contained in self.

        The cooefficients of the point don't matter, only the value.

        INPUT:

        - ``point`` - PointWithCoefficients, must be in the unit
          interval $[0, 1)$

        OUTPUT:

        - boolean - True if ``point`` is contained in self,
          False if not

        EXAMPLES::

            sage: i = Interval(PointWithCoefficients(0.25, [0, 1]), \
                PointWithCoefficients(0.5,[1, -2]),'open','closed')
            sage: i.contains(PointWithCoefficients(0.1, [1,1]))
            False
            sage: i.contains(PointWithCoefficients(0.3, [-2,6]))
            True
            sage: i.contains(PointWithCoefficients(0.25, [-1,5]))
            False
            sage: i.contains(PointWithCoefficients(0.5, [1,7]))
            True

            sage: j = Interval(PointWithCoefficients(0.5, [0, 1]), \
                PointWithCoefficients(0.25,[1, -2]),'closed','open')
            sage: j.contains(PointWithCoefficients(0.1, [1,1]))
            True
            sage: j.contains(PointWithCoefficients(0.3, [-2,6]))
            False
            sage: j.contains(PointWithCoefficients(0.25, [-1,5]))
            False
            sage: j.contains(PointWithCoefficients(0.5, [1,7]))
            True

            sage: k = Interval(PointWithCoefficients(0.3, (1,1)),\
                    PointWithCoefficients(0.3, (1,1)))
            sage: k.contains(PointWithCoefficients(0.3, (1,1)))
            False

            sage: l = Interval(PointWithCoefficients(0.3, (1,1)),\
                    PointWithCoefficients(0.3, (1,1)), 'closed', \
                    'closed')
            sage: l.contains(PointWithCoefficients(0.3, (1,1)))
            True

        """
        if not 0 <= point.value < 1:
            raise ValueError("Only points in the unit interval can be"
                    " tested for containment")
        if self[0].value <= self[1].value:
            if self._less(self[0].value, point.value, 
                    is_equality_allowed = (self._left_openness == 
                        'closed')) and\
                    self._less(point.value, self[1].value, 
                            is_equality_allowed = 
                            (self._right_openness == 'closed')):
                return True
            else:
                return False
        if self._less(self[1].value, point.value, 
                is_equality_allowed = 
                (self._right_openness == 'open')) and\
                self._less(point.value, self[0].value, 
                    is_equality_allowed =
                        (self._left_openness == 'open')):
            return False
        else:
            return True


def has_constant_sign(vec):
    if abs(vec[0]) < epsilon:
        return False
    is_first_positive = (vec[0] > 0)
    for x in vec:
        if abs(x) < epsilon  or (x > 0) != is_first_positive:
            return False
    return True

def get_good_eigendata(transition_matrix, is_twisted):
    ev = transition_matrix.eigenvectors_right()
    ret_list = []
    for x in ev:
        if abs(x[0].imag()) < epsilon and x[0] > 0 \
                and abs(x[0] - 1) > epsilon:
            for vec in x[1]:
                if has_constant_sign(vec):
                    norm = sum([abs(y) for y in vec])
                    if not is_twisted:
                        norm -= abs(vec[-1])
                    normalized_vec = [abs(y / norm) for y in vec]
                    ret_list.append((x[0], normalized_vec))
    return ret_list

class SaddleConnectionError(Exception):
    pass

class RestrictionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

futyi = 'Anyaddal szorakozzal.'

class Foliation(SageObject):
    def __init__(self, involution, lengths, twist = None):
        if not involution.is_bottom_side_moebius():
            if twist == None:
                raise ValueError('The twist must be specified '
                'unless the bottom side is a Moebius band.')

        if isinstance(lengths, list):
            if len(lengths) != len(involution.alphabet()):
                    raise ValueError('Bad number of lengths')
            lcopy = {letter: lengths[involution.index(letter)] 
                    for letter in involution.alphabet()}
                        
        if isinstance(lengths, dict):
            if set(involution.alphabet()) != set(lengths.keys()):
                raise ValueError('Invalid length specification')     
            lcopy = dict(lengths)

        if any(v <= 0 for v in lcopy.values()):
            raise ValueError('Lengths must be positive')

        if involution.is_bottom_side_moebius():
            lcopy['JOKER'] = sum(lcopy.values())
            twist = 0
        if involution.is_torus():
            lcopy['JOKER'] = 1

        totals = [sum(lcopy[letter] for letter in involution[i]) 
                for i in range(2)]

        if abs(totals[0] - totals[1]) > epsilon:
            raise ValueError('The total length on the top and '
                    'bottom are inconsistent.')
        
        #adjusting lengths in case they 
        #differ slightly on top/bottom
        for j in range(len(involution[0])):
            if involution.pair((0, j))[0] == 0:
                lcopy[involution[i][j]] += \
                        (totals[1] - totals[0])/2
                break

        self._lengths = {}
        for letter in lcopy:
            self._lengths[letter] = PointWithCoefficients(\
                    lcopy[letter]/totals[1], 
                    self._basis_vector(len(lcopy) + 1,
                        involution.index(letter)))

        self._divpoints = [[PointWithCoefficients(0,
            [0] * (len(lcopy) + 1))]
            for i in range(2)] 

        for i in range(2):
            for j in range(len(involution[i]) - 1):
                self._divpoints[i].append(self._divpoints[i][-1] +
                        self._lengths[involution[i][j]])

        self._divvalues = [[x.value for x in self._divpoints[i]] 
            for i in range(2)]

        preimage_of_zero = self._mod_one(-twist/totals[1])
        containing_int = bisect_left(self._divvalues[1], 
                preimage_of_zero) % len(involution[1])
        self._involution = involution.rotated(0, -containing_int)
        self._twist = PointWithCoefficients(self._mod_one(\
                self._divpoints[1][containing_int].value - 
                preimage_of_zero), [0] * len(lcopy) + [1])
        self._divpoints[1] = [self._twist]
        for j in range(len(involution[1]) - 1):
            self._divpoints[1].append(self._divpoints[1][-1] +
                    self._lengths[self._involution[1][j]])

        self._divvalues[1] = [x.value for x in self._divpoints[1]] 

        self._length_twist_vector = [0] * len(lcopy)
        self._length_twist_vector.append(self._twist.value)
        for letter in self._lengths:
            self._length_twist_vector[self._involution.index(letter)]\
                    = self._lengths[letter].value
        self._length_twist_vector = vector(self._length_twist_vector)

    @staticmethod
    def _mod_one(x):
        return x - floor(x)

    @staticmethod
    def _basis_vector(n, k):
        l = [0] * n
        l[k] = 1
        return l

    @classmethod
    def orientable_arnoux_yoccoz(self, genus):
        sf = arnoux_yoccoz_factor(genus)
        l = [1 / sf^(i + 1) for i in range(genus)] * 2
        return Foliation(Involution.orientable_arnoux_yoccoz(genus),
                sorted(l, reverse = True), twist = 1)

    @classmethod
    def nonorientable_arnoux_yoccoz(self, genus):
        sf = arnoux_yoccoz_factor(genus - 1)
        return Foliation(Involution.nonorientable_arnoux_yoccoz(\
                genus), [1/sf^i for i in range(genus - 1)])

    @classmethod
    def RP2_arnoux_yoccoz(self):
        sf = arnoux_yoccoz_factor(3)
        return Foliation(Involution.RP2_arnoux_yoccoz(), 
                [1/sf + 1/sf^2, 1/sf^2 + 1/sf^3, 1/sf + 1/sf^3])
                

    def __eq__(self, other):
        return self._involution == other._involution and\
                abs(self._length_twist_vector -
                    other._length_twist_vector) < epsilon

    def __repr__(self):
        if self._involution.is_bottom_side_moebius():
            return "{0}\nLengths:{1}".format(\
                    self._involution, self._length_twist_vector)
        return "{0}\nLengths:{1}\nTwist:{2}".format(\
                self._involution, self._length_twist_vector[:-1], 
                self._twist.value)


    TransitionData = namedtuple('TransitionData', 'tr_matrix,new_inv')

    def new_foliation(self, transition_data):
        new_vector = transition_data.tr_matrix * self._length_twist_vector 
        if self._involution.is_bottom_side_moebius():
            return Foliation(transition_data.new_inv, 
                new_vector.list())
        else:
            return Foliation(transition_data.new_inv, 
                new_vector[:-1].list(), new_vector[-1])

    def _simple_transformation_data(self, new_involution, 
            twist_row = None):
        m = matrix(self._num_coeffs)
        for letter in self._involution.alphabet():
            m[new_involution.index(letter),
                    self._involution.index(letter)] = 1
        if not self._involution.is_bottom_side_moebius():
            m[-1] = twist_row
        return self.TransitionData(tr_matrix = m, 
                new_inv = new_involution)

    def _rotation_data(self, k):
        n = len(self._divpoints[0])
        k = k % n
        twist_row = [0] * self._num_coeffs
        twist_row[-1] = 1
        for i in range(k):
            twist_row[self._involution.index(\
                    self._involution[0][n - 1 - i])] += 1
        return self._simple_transformation_data(\
                self._involution.rotated(k, 0),
                twist_row = twist_row)

    def _reverse_data(self):
        return self._simple_transformation_data(\
                self._involution.reversed(),
                twist_row = [0] * (self._num_coeffs - 1) + [-1])

    def rotated(self, k):
        return self.new_foliation(self._rotation_data(k))

    def reversed(self):
        return self.new_foliation(self._reverse_data())

    def _in_which_interval(self, value, side):
        assert(side == 0 or not self._involution.is_bottom_side_moebius())
        interval = bisect(self._divvalues[side], value)
        if interval == 0:
            to_check = {self._divvalues[side][0]}
        elif interval == len(self._divpoints[side]):
            to_check = {self._divvalues[side][interval - 1],
                    self._divvalues[side][0] + 1}
        else:
            to_check = {self._divvalues[side][k]
                    for k in {interval - 1, interval}}

        if any(abs(value - x) < epsilon for x in to_check):
            raise SaddleConnectionError()
        return (interval - 1) % len(self._divpoints[side])

    def _map(self, point, side, containing_interval = None):
        if side == 1 and self._involution.is_bottom_side_moebius():
            return (1, point.half_added())
        if containing_interval == None:
            containing_interval = self._in_which_interval(\
                    point.value, side)
                    
        flipped = self._involution.is_flipped((side, 
                containing_interval)) 
        new_int= self._involution.pair((side, 
            containing_interval))
        diff = point - self._divpoints[side][containing_interval]
        if not flipped:
            return (new_int[0], 
                    (self._divpoints[new_int[0]][new_int[1]] + 
                        diff).mod_one())
        else:
            return (new_int[0],
                    (self._divpoints[new_int[0]][(new_int[1] + 1) % 
                        len(self._divpoints[side])] - diff).\
                                mod_one())

    IntersectionData = namedtuple('IntersectionData',
            'point, side, is_flipped')

    def _first_intersection(self, side, pos, intervals = None,
            crossings = None):
        assert(not self._involution.is_bottom_side_moebius() or 
                side == 0)
        point = self._divpoints[side][pos]

        if intervals == None:
            def terminate(p):
                return flipped == True
        else:
            def terminate(p):
                return any(interval.contains(p)
                        for interval in intervals)

        flipped = False

        while not terminate(point):
            if crossings != None:
                crossings.append(point.value)
            side = (side + 1) % 2
            if self._involution.is_bottom_side_moebius() and side == 1:
                point = point.half_added()
            else:
                interval = self._in_which_interval(point.value, side)
                if self._involution.is_flipped(side, interval):
                    flipped = not flipped
                (side, point) = self._map(point, side, interval)

        int_data = self.IntersectionData(point = point, 
               side = side, is_flipped = flipped)

        return int_data

    DistanceData = namedtuple('DistanceData', 'side, distance, '
            'is_flipped, orig_side, orig_pos')

    def _get_involution(self, sorted_distances, is_lift):
        done = set()
        flips = set()
        if is_lift:
            remaining_letters = \
                    range((len(sorted_distances[0]) + 
                        len(sorted_distances[1]))/2, 0, -1)
        else:
            remaining_letters = self._involution.alphabet()
        inv_list = [[0] * len(sorted_distances[i]) for i in {0,1}]
        for side in {0, 1}:
            for i in range(len(sorted_distances[side])):
                if (side, i) in done:
                    continue
                letter = remaining_letters.pop()
                on_the_right = True
                d = sorted_distances[side][i]
                new_side = d.orig_side
                new_i = d.orig_pos
                if d.is_flipped: #orintation is reversed
                    on_the_right = False
                    new_i = (new_i - 1) % \
                            len(self._involution[new_side])
                (new_side, new_i) = self._involution.pair(\
                        (new_side, new_i))
                if self._involution.is_flipped(new_side, new_i):
                    on_the_right = not on_the_right
                if not on_the_right:
                    new_i = (new_i + 1) % \
                            len(self._involution[new_side])
                (new_side, new_i) = next((a, b) for a in {0, 1}
                        for b in range(len(sorted_distances[a]))
                        if (sorted_distances[a][b].orig_side,
                        sorted_distances[a][b].orig_pos) == 
                        (new_side, new_i) and (not is_lift or 
                            sorted_distances[a][b].is_flipped !=
                            on_the_right))
                if sorted_distances[new_side][new_i].is_flipped:
                    on_the_right = not on_the_right
                if not on_the_right:
                    new_i = (new_i - 1) % \
                            len(sorted_distances[new_side])
                inv_list[side][i] = inv_list[new_side][new_i] = \
                        letter
                done.add((side, i)); done.add((new_side, new_i))
                if not on_the_right:
                    flips.add(letter)
        if inv_list[1] == []:
            del inv_list[1]
        return Involution(*inv_list, flips = flips)

    @staticmethod
    def _get_sorted_distances(distances):
        return [sorted([x for x in distances if 
            x.side == side], key = lambda s: s.distance.value)
                for side in {0, 1}]

       
    def _create_transition_data(self, intervals, total, 
            distance_getter):
        intersections = self._get_intersections(intervals)

        distances = [distance_getter(intersections, 
            side, pos) for side in {0, 1}
            for pos in range(len(intersections[side]))]

        sorted_distances = self._get_sorted_distances(distances)
        new_involution = self._get_involution(sorted_distances,
                is_lift = False)
        #print intersections
        #print distances
        #print sorted_distances

        m = matrix(self._num_coeffs)
        done = set()
        for side in {0,1}:
            for i in range(len(sorted_distances[side])):
                letter = new_involution[side][i]
                if letter in done:
                    continue
                diff = sorted_distances[side][(i + 1) % 
                        len(sorted_distances[side])].distance -\
                                sorted_distances[side][i].distance
                if diff.value < 0:
                    diff += total
                m[len(done)] = diff.coefficients
                done.add(letter)
        if len(sorted_distances[1]) != 0:
            m[-1] = sorted_distances[1][0].distance.coefficients -\
                    sorted_distances[0][0].distance.coefficients

        return self.TransitionData(tr_matrix = m,
                new_inv = new_involution)

    def foliation_orientable_double_cover(self):
        """
        Returns the orienting double cover of the foliation.

        If the foliation is already orientable (i.e. the surface
        is orientable and the holonomy is trivial, or the surface
        is non-orientable and the holonomy is Z_2, or equivalently
        no interval is flipped), then there is nothing to do.

        So assume that it is not the case and there is at least one
        flipped pair of intervals. Our transverse curve has trivial
        holonomy, so it has two lifts. A leaf segment
        with endpoints in a not flipped pair of intervals lifts to 
        two segments, each having endpoints in one of the curves.
        A leaf segments with endpoints in a flipped pair of intervals
        lifts to two segments, one starting from the first curve, and
        ending in the second, the other the other way around.

        So a leaf segment with endpoints in the first lift
        of our curve is either just a lift of a segment corresponding
        to a not-flippped interval as above, or a concatenation of
        many segments, starting and ending with a segment that comes
        from a flipped interval, and with segments coming from
        not-flipped intervals in the middle. (One first travels from
        the first curve to the second, then returns a bunch of times
        to the second, and finally returns to the first.)

        So division points of the lifted foliation are the old 
        division points union the following points. From each 
        singularity travel along a separatrix, and stop just after
        going through a strip of a flipped pair of intervals.
        These endpoints are also new division points, so one has twice
        as many division points for the lift foliation as for the
        original one.




        """
        if len(self._involution.flips()) == 0:
            return self
        flipped_intersections = self._get_intersections()
        distances = []
        for side in range(2):
            n = len(self._divpoints[side])
            for pos in range(n):
                if self._involution[side][pos] == self._involution\
                        [side][(pos - 1) % n] and self._involution.\
                        is_flipped((side, pos)):
                            continue
                distances.append(self.DistanceData(side = side,
                    distance = self._divpoints[side][pos],
                    is_flipped = False,
                    orig_side = side,
                    orig_pos = pos))
                int_data = flipped_intersections[side][pos]
                distances.append(self.DistanceData(side = int_data.\
                        side,
                        distance = int_data.point,
                        is_flipped = True,
                        orig_side = side,
                        orig_pos = pos))
        sorted_distances = self._get_sorted_distances(distances)

        new_involution = self._get_involution(sorted_distances,
                is_lift = True)


        lengths = {new_involution[side][pos] :(sorted_distances[side]\
                [(pos + 1) % len(sorted_distances[side])].distance -
                sorted_distances[side][pos].distance).mod_one().value
                for side in range(2) for pos in 
                range(len(sorted_distances[side]))}
        return Foliation(new_involution, lengths)
 

    def surface_orientable_double_cover(self):
        """
        Returns the double cover of the foliation that orients the
        surface.

        If the surface is already orientable, there is nothing to do.
        If not, then the our transverse curve, which separates the
        non-orientable surface to a Moebius band and another surface
        has two lifts. The Moebius band lifts to an annulus bounded
        by these two lifts. 
        
        If the other surface is orientable (i.e.
        all intervals are flipped),
        it has two lifts homeomorphic to itself, and they are
        glued to the two sides of the annulus. The folition inside
        the annulus has the effect of a twist of 1/2 between the 
        surfaces on the two sides. All pairs of intervals remain
        flipped.

        If the complement of the Mobius band is non-orientable (i.e.
        there is at least one interval which is not twisted), then
        its double cover has one components with two bounding circles
        that are glued to the two boundary circles of the annulus.
        Intervals that were flipped stay flipped and will appear
        on both sides of the transverse curve. Intervals that were
        not flipped will turn into a pair of intervals on different
        sided, still not flipped (i.e. strips connecting 
        different sides of the annulus).

        In any case, one only has to change the interval pairings for
        the new involution, but the flips are inherited, and add a 
        1/2 twist.

        OUTPUT:

        - Foliation - if self is already orientable, the output is
          itself. Otherwise it is the foliation obtained by taking
          the double cover that orients the surface.

        EXAMPLES::
            
            sage: f = Foliation.nonorientable_arnoux_yoccoz(4)
            sage: g = Foliation.orientable_arnoux_yoccoz(3)
            sage: f.surface_orientable_double_cover() == g
            True
            
            sage: f = Foliation.nonorientable_arnoux_yoccoz(10)
            sage: g = Foliation.orientable_arnoux_yoccoz(9)
            sage: f.surface_orientable_double_cover() == g
            True
        """
        if self._involution.is_surface_orientable():
            return self
        assert(self._involution.is_bottom_side_moebius())
        n = len(self._divpoints[0])
        alphabet = range(n, 0, -1)
        inv_list = [[0] * n for i in range(2)]
        done = set()
        flips = set()
        lengths = {}
        for i in range(n):
            if i in done:
                continue
            j = self._involution.pair((0, i))[1]
            letter = alphabet.pop()
            letter2 = alphabet.pop()
            lengths[letter] = lengths[letter2] = \
                    self._lengths[self._involution[0][i]].value
            if self._involution.is_flipped((0, i)):
                inv_list[0][i] = inv_list[0][j] = letter
                inv_list[1][i] = inv_list[1][j] = letter2
                flips.add(letter)
                flips.add(letter2)
            else:
                inv_list[0][i] = inv_list[1][j] = letter
                inv_list[0][j] = inv_list[1][i] = letter2
            done.add(i); done.add(j)
        return Foliation(Involution(*inv_list, flips = flips), 
                lengths, twist = 1/2)


    def _num_letters(self):
        return len(self._involution.alphabet())

    def _check_not_flipped(self, side, pos):
        if self._involution.is_flipped(side, pos):
            raise RestrictionError('The specified interval should '
                    'not be flipped, but it is.')

    def _check_flipped(self, side, pos):
        if not self._involution.is_flipped(side, pos):
            raise RestrictionError('The specified interval '
                    'should be flipped, but it isn\'t.')

    def _check_same_side(self, side, pos):
        if self._involution.pair((side, pos))[0] != side:
            raise RestrictionError('The interval pair should '
                    'be on the same same side, but it isn\'t.')

    def _check_different_side(self, side, pos):
        if self._involution.pair((side, pos))[0] == side:
            raise RestrictionError('The interval pair should '
                    'be on different sides, but it isn\'t.')

    def _assert_different_pair(self, side1, pos1, side2, pos2):
        if (side1, pos1) in {(side2, pos2),
                self._involution.pair((side2, pos2))}:
            raise RestrictionError('The chosen pairs should be '
                    'different, but they are not.')

    def _get_intersections(self, intervals = None):
        return [[self._first_intersection(side, pos, 
            intervals = intervals)
                for pos in range(len(self._divpoints[side]))]
                for side in {0,1}]

    def _get_right_endpoint(self, side, pos):
        return self._divpoints[side][(pos + 1) %
                len(self._divpoints[side])]

    def _get_left_endpoint(self, side, pos):
        return self._divpoints[side][pos]


    def _restrict_not_flipped_same_side(self, side, pos):
        self._check_not_flipped(side, pos)
        self._check_same_side(side, pos)
        #n = len(self._divpoints[side])
        left_endpoint = self._get_right_endpoint(side, pos)
        side2, pos2 = self._involution.pair((side, pos))
        right_endpoint = self._get_right_endpoint(side2, pos2)

        interval = Interval(left_endpoint, right_endpoint)
        total = interval.length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            distance = (int_data.point - left_endpoint).mod_one()
            if int_data.side == 1:
                distance += total
            return self.DistanceData(side = 0,
                distance = distance,
                is_flipped = int_data.is_flipped,
                orig_side = side,
                orig_pos = pos) 

        return self._create_transition_data([interval], total + total,
                get_distance_data)


    def _restrict_not_flipped_different_side(self, side, pos):
        self._check_not_flipped(side, pos)
        self._check_different_side(side, pos)
        
        left_endpoint = self._get_right_endpoint(side, pos)
        side2, pos2 = self._involution.pair((side, pair))
        right_endpoint = self._get_right_endpoint(side2, pos2)

        interval = Interval(left_endpoint, right_endpoint)

        total = interval.length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            return self.DistanceData(side = int_data.side,
                distance = (int_data.point -left_endpoint).mod_one(),
                is_flipped = int_data.is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data([interval], total)

           
    def _restrict_flipped_two_sided(self, side1, pos1,
            side2, pos2):
        self._assert_different_pair(side1, pos1, side2, pos2)

        Endpoint = namedtuple('Endpoint', 
                'point, end, is_closed, side')

        ints =[[(side1, pos1),
            self._involution.pair((side1, pos1))], [(side2, pos2), 
                self._involution.pair((side2, pos2))]]
        endpoints = []
        for i in {0,1}:
            self._check_flipped(*ints[i][0])
            endpoints.extend([Endpoint(point=self._get_left_endpoint(\
                    *ints[i][0]),end = 'left', is_closed = (i == 1), 
                    side = ints[i][0][0]), 
                    Endpoint(\
                    point = self._get_right_endpoint(*ints[i][1]),
                    end = 'right', is_closed = (i == 1),
                    side = ints[i][1][0])])

        endpoints.sort(key = lambda x: x.point.value)
        while endpoints[0].end == 'right' or\
                endpoints[1].end == 'left':
                    endpoints = endpoints[-1:] + endpoints[:-1]


        def create_interval(i, j):
            return Interval(endpoints[i].point,
                endpoints[j].point,
                left_openness = endpoints[i].is_closed,
                right_openness = endpoints[j].is_closed)


        if endpoints[1].end == 'right':
            intervals = [create_interval(0, 1), create_interval(2, 3)]
        else:
            if endpoints[1].side == endpoints[2].side:
                interval1 = create_interval(1, 2)
                interval2 = create_interval(0, 3)
            else:
                interval2 = create_interval(0, 2)
                interval1 = create_interval(1, 3)
            if endpoints[1].side == 0:
                intervals = [interval1, interval2]
            else:
                intervals = [interval2, interval1]

        total = intervals[0].length() + intervals[1].length() 

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            if int_data.side == 0:
                if intervals[0].contains(int_data.point):
                    distance = (int_data.point -
                            points_with_coeffs[intervals[0][0]]).\
                                    mod_one()
                    new_side = 0
                    is_flipped = int_data.is_flipped
                else:
                    distance = total - (int_data.point - 
                            points_with_coeffs[intervals[1][0]]).\
                                    mod_one()
                    new_side = 1
                    is_flipped = not int_data.is_flipped
            else:
                if intervals[1].contains(int_data.point):
                    distance = total - (int_data.point - 
                            points_with_coeffs[intervals[1][0]]).\
                                    mod_one()
                    new_side = 0
                    is_flipped = not int_data.is_flipped
                else:
                    distance = (int_data.point -
                            points_with_coeffs[intervals[0][0]]).\
                                    mod_one()
                    new_side = 1
                    is_flipped = int_data.is_flipped
            return self.DistanceData(side = new_side,
                distance = distance,
                is_flipped = is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data(intervals, total,
                get_distance_data)

    @staticmethod
    def _get_distances(intersections, getter_function):
        return [getter_function(side, pos) for side in {0, 1}
            for pos in range(len(intersections[side]))]



    def _restrict_flipped_same_side_one_sided(self, side1, pos1,
            side2, pos2):
        self._check_flipped(side1, pos1)
        self._check_flipped(side2, pos2)
        self._assert_different_pair(side1, pos1, side2, pos2)
        assert(self._involution.is_bottom_side_moebius())
        
        n = len(self._divpoints[0])
        center_left = self._divpoints[0][(pos1 + 1) % n]
        center_right = self._divpoints[0][pos2]
        other_left = self._divpoints[0][self._involution.pair(\
                (side1, pos1))[1]]
        other_right = self._divpoints[0][(self._involution.pair(\
                (side2, pos2))[1] + 1) % n]
        center_int = Interval(center_left, center_right,
                left_openness = True, right_openness = True)
        if center_int.contains(other_left) or \
                center_int.contains(other_right):
                    raise RestrictionError('Specified transverse '
                    'curve does not exist. Combinatorically invalid '
                    'choice for the intervals.')
        if (other_right - other_left).mod_one().value <= 1/2:
            raise RestrictionError('Specified transverse curve does '
            'not exist. The length parameters doesn\'t allow '
            'closing the curve to a transverse curve.')
        other_left = other_left.half_added()

        other_int = Interval(other_left, other_right,
                left_openness = False, right_openness = False)

        intervals = [center_int, other_int]

        total = intervals[0].length() + intervals[1].length()

        def get_distance_data(intersections, side, pos):
            int_data = intersections[side][pos]
            if int_data.side == 0:
                if intervals[0].contains(int_data.point):
                    distance = (int_data.point - center_left).\
                            mod_one()
                    is_flipped = int_data.is_flipped
                else:
                    distance = total + total - (int_data.point - 
                            other_left).mod_one()
                    is_flipped = not int_data.is_flipped
            else:
                if intervals[1].contains(int_data.point):
                    distance = total - (int_data.point - 
                            other_left).mod_one()
                    is_flipped = not int_data.is_flipped
                else:
                    distance = (int_data.point - center_left).\
                            mod_one() + total
                    is_flipped = int_data.is_flipped
            return self.DistanceData(side = 0,
                distance = distance,
                is_flipped = is_flipped,
                orig_side = side,
                orig_pos = pos)

        return self._create_transition_data(intervals, total + total,
                get_distance_data)


#    def _find_pseudo_anosov_candidates(self, depth,


#f = Foliation(Involution('a a b b c c'), [1, 2, 3])

