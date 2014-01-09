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

    The mapping class group of the closed torus and the 
    once-punctured torus are the same, and since we are interested in
    constructed pseudo-anosovs, we won't complicate the code with
    treating the case of the closed torus separately. But we can 
    represent once-punctured tori as follows::

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
    def __init__(self, top_letters, bottom_letters = None, 
            flips = []):
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
        if self.is_bottom_side_moebius():
            done = {(1,0), (1,1)}
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


    def _repr_(self):
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

        """
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

    def singularity_partition(self):
        """
        Returns the singularity partition of self.

        OUTPUT:

        - list of lists - each element of the list is a list corresponding to
          a singularity, and each such list contains the tuples of positions
          that are being identified

        EXAMPLES::

            sage: i = Involution('a a b b c c')
            sage: sp = i.singularity_partition()
            sage: len(sp) == 1
            True
            sage: set(sp[0]) == {(0,2), (0,3), (0, 4), (0, 5), (0, 0), (0, 1)}
            True

            sage: i = Involution('a b c d', 'b a d c', flips = 'ac')
            sage: sp = i.singularity_partition()
            sage: len(sp) == 2
            True
            sage: set(sp[0]) == {(1, 1), (0, 2), (1, 0), (0, 1)}
            True
            sage: set(sp[1]) == {(0, 0), (1, 3), (0, 3), (1, 2)}
            True

        """
        return list(self._singularity_partition)

    def which_singularity(self, pos):
        """
        Returns the index of the singularity for the beginning of each
        interval.

        There is a singularity of the foliation on the leaf containing the
        left endpoint of each interval for any suspension. There may be only
        one singularity of all the vertices are identified, or more if not.
        If there are $n$ singularities ($n\ge 1$), we assign 0, 1, ..., $n-1$
        to them in some order, this is called its index. The index is 
        therefore well-defined only up to a permutation of these values.

        INPUT:

        - ``pos`` - a tuple encoding the position. The first
          coordinate is 0 or 1 depending on whether it is a top
          or bottom interval. The second coordinate is the
          index of the interval in that row.
        
        OUTPUT:

        - integer - the index of the specified singularity. 

        EXAMPLES:

        The following Involution has 2 singularities, one has 5 prongs, the
        other 1 prong. The position of the 1-prong singularity is at (1,0).
        Here is a possible output:

            sage: i = Involution('a a b', 'c b c', flips = 'c')
            sage: i.singularity_type()
            (5, 1)
            sage: i.which_singularity((0,0)) 
            0
            sage: i.which_singularity((0,1)) 
            0
            sage: i.which_singularity((0,2)) 
            0
            sage: i.which_singularity((1,0)) 
            1
            sage: i.which_singularity((1,1)) 
            0
            sage: i.which_singularity((1,2)) 
            0

        """
        sp = self._singularity_partition
        for i in range(len(sp)):
            if pos in sp[i]:
                return i
        raise ValueError("Invalid singularity specification.")
         
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
        t = sorted([x for x in map(len, 
            self._singularity_partition)], reverse = True)
        return tuple(t)
    
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
            sage: i.is_bottom_side_moebius()
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

    def _repr_(self):
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

class Interval(SageObject):
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

    def _repr_(self):
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


def make_positive(vec):
    """
    Returns a parallel vector to a given vector with all positive
    coordinates if possible.

    INPUT:

    - ``vec`` - a vector or list or tuple of numbers

    OUTPUT:

    - vector - a vector with positive coordinates if this is possible,
      otherwise -1

    EXAMPLES:

    If all coordinates of the input vector are already positive, the
    same vector is returned::

        sage: v = vector([1, 3, 5])
        sage: make_positive(v) == v
        True

    If all are negative, its negative is returned::

        sage: make_positive((-1, -5))
        (1, 5)

    Even if the coordinates are complex, but have very small imaginary
    part as a result of an approximate eigenvector calculation for 
    example, the coordinates are treated as real::

        sage: make_positive((40.24 - 5.64e-16*I, 1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)
        sage: make_positive((-40.24 - 5.64e-16*I, -1.2 + 4.3e-14*I))
        (40.2400000000000, 1.20000000000000)

    If there is a complex coordinate which is not negligible, -1 is
    returned::

        sage: make_positive((-40.24 - 5.64e-6*I, -1.2 + 4.3e-14*I))
        -1

    If one coordinate is zero, or very close to zero, -1 is returned::

        sage: make_positive((1, 0, 2))
        -1
        sage: make_positive((-40.24e-15*I - 5.64e-16*I, -1.2))
        -1

    If there are both negative and positive coordinates, -1 is 
    returned::

        sage: make_positive((-3, 4, 5))
        -1
        
    """
    if any(abs(x) < epsilon for x in vec) or \
        any(abs(x.imag()) > epsilon for x in vec):
            return -1
    newvec = vector([x.real() for x in vec])
    if vec[0].real() < 0:
        newvec = -newvec
    if any(x < 0 for x in newvec):
        return -1
    return newvec

def get_good_eigendata(transition_matrix, is_twisted):
    """
    


    """
    ev = transition_matrix.eigenvectors_right()
    ret_list = []
    for x in ev:
        if abs(x[0].imag()) < epsilon and x[0].real() > 0 \
                and abs(x[0].real() - 1) > epsilon:
            for vec in x[1]:
                newvec = make_positive(vec)
                if newvec != -1:
                    norm = sum([abs(y) for y in newvec])
                    if is_twisted:
                        norm -= abs(newvec[-1])
                    normalized_vec = [abs(y / norm) for y in newvec]
                    ret_list.append((x[0].real(), normalized_vec))
    return ret_list

class SaddleConnectionError(Exception):
    pass

class RestrictionError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

_tikzcolors = ["red", "green", "blue", "cyan", "magenta", "yellow", 
    "gray", "brown", "lime", "olive", "orange", "pink", "purple", 
    "teal", "violet"]

def _tikzcolor(n):
    return _tikzcolors[n % len(_tikzcolors)]
 
class Foliation(SageObject):
    """
    A measured foliation on a surface of finite type.

    Given a measured foliation we consider a two-sided simple closed
    curve on the surface which is transverse to the foliation. This 
    gives an interval exchange of the circle which is represented by
    an Involutions, length parameters, and a twist parameter. The 
    measure is always normalized so that the length of the curve is 1.

    We consider only surfaces with negative Euler characteristic, so
    foliations on the closed tori and closed Klein bottle are not 
    represented. By the Euler-Poincare theorem the foliation must have
    an even number of separatrices, say $2k$, where $k > 0$. There
    are two cases: all these separatrices hit our simple closed curve
    from above, or there is at least one separatrix hitting each side.

    In the first case the bottom side of the curve is a Moebius band,
    so we don't need a twist parameter. But even though we wouldn't
    need a length parameter for the bottom side to uniquely
    represent the foliation (it can be expressed in terms of the 
    length parameters above the curve), we consider 1/2 as the twist
    paratemeter for the bottom side with interval exchange imagined
    as 'z z'. Thus we have $k + 1$ parameters.

    In the other case there are $k$ length parameters (which again
    may not be independent), but a twist parameter here is 
    necessary which is $k + 1$ parameters again. This consistency may
    be useful for dealing with transition matrices between parameters
    of the same foliation, but looked at from the perspective of
    different simple closed curves.

    INPUT:

    - ``involution`` - an Involution, serving as the combinatorial
      description of the interval exchange

    - ``lenghts`` - as in iet.IntervalExchangeTransformation,
      this is either a list or dictionary of the length parameters.
      If it is a list, then the lengths are assigned to intervals
      in the order of their appearance in the Involution (their
      "index"). If it's a dict, then lengths should be assigned to
      the names of the intervals in the Involution.

    - ``twist`` - a real number, the twist parameter. This can 
      (and recommended to) be omitted if the bottom side is a 
      Moebius band, otherwise it is mandatory, because defaulting to
      zero would lead to an immediate saddle connection.

    EXAMPLES:
    
    Here are two different but equivalent definitions::

        sage: i = Involution('a a b b', 'c c', flips = 'b')
        sage: f = Foliation(i, [1, 2, 3], 1/2); f
        a a -b -b
        c c
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/12

        sage: g = Foliation(i, {'a':1,'b':2,'c':3}, 1/2); f
        a a -b -b
        c c
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/12

        sage: f == g
        True

    Omitting the twist if the bottom side is not a Moebius band
    throws an error::

        sage: Foliation(i, [1, 2, 3])
        Traceback (most recent call last):
        ...
        ValueError: The twist must be specified unless the bottom side is a Moebius band.

    When the bottom side is a Moebius band, it is okay to omit the 
    twist::

        sage: Foliation(Involution('a a b b c c'), [2, 2, 1])
        a a b b c c
        Moebius band
        Lengths: (1/5, 1/5, 1/10)

    The sum of the lengths on top and the sum at the bottom should
    (approximately) equal::

        sage: Foliation(i, [1, 1, 3], 1/2)
        Traceback (most recent call last):
        ...
        ValueError: The total length on the top and bottom are inconsistent.

    If they are just approximately equal, then on top of the
    normalization, the lengths are adjusted so the sums on top and 
    bottom are equal::

        sage: Foliation(i, [1, 2, 3.0000000000001], 1/2)
        a a -b -b
        c c
        Lengths: (0.166666666666661, 0.333333333333322, 0.500000000000017)
        Twist: 0.0833333333333306

    In reality, the same twisted interval exchange transformation can
    be represented by different Involutions and twist parameters. 
    For instance depending on which separatrix is chosen to be zero,
    the top letters can be 'a a b b', 'a b b a', 'b b a a' or
    'b a a b'. But even if one chooses one of these, varying the twist
    and the Involution can result in the same Foliation::

        sage: i = Involution('a b c', 'a c b'); i
        a b c
        a c b
        sage: f = Foliation(i, [1, 2, 3], 1); f
        a b c
        a c b
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/6

        sage: j = Involution('a b c', 'c b a'); j
        a b c
        c b a
        sage: g = Foliation(j, [1, 2, 3], 2); g
        a b c
        a c b
        Lengths: (1/6, 1/3, 1/2)
        Twist: 1/6

        sage: f == g
        True

    So rotate the bottom row of the Involution and change the twist
    if necessary to obtain the smallest possible positive twist. One
    can easily check that the Foliations f and g above have an 
    immediate saddle connection therefore they are not candidates for
    pseudo-anosov stable foliations. We don't check this in the
    constructor, only later when separatrices are lengthened to find
    other curves.

    """ 

    @staticmethod
    def _mod_one(x):
        """
        Returns a number modulo 1.

        INPUT:

        - ``x`` - a real number

        OUTPUT:

        - a real number of the same type as the input

        TESTS::

            sage: Foliation._mod_one(2.5)
            0.500000000000000
            sage: Foliation._mod_one(-1.7)
            0.300000000000000
            sage: Foliation._mod_one(7/6)
            1/6
            sage: Foliation._mod_one(-1/6)
            5/6
            sage: a = QQbar(sqrt(2)); a
            1.414213562373095?
            sage: Foliation._mod_one(a)
            0.4142135623730951?

        """
        return x - floor(x)

    @staticmethod
    def _basis_vector(n, k):
        """
        Returns a standard basis vector of a vector space.

        INPUT:

        - ``n`` - positive integer, the dimension of the vector space

        - ``k`` - 0 <= k < n, the index of the only coordinate which
          is 1.

        OUTPUT:

        - list - the basis vector in as a list

        EXAMPLES::

            sage: Foliation._basis_vector(5, 1)
            [0, 1, 0, 0, 0]
            sage: Foliation._basis_vector(7, 6)
            [0, 0, 0, 0, 0, 0, 1]

        """
        l = [0] * n
        l[k] = 1
        return l

    def __init__(self, involution, lengths, twist = None):
        """
        TESTS::

            sage: i = Involution('a b c', 'a c b'); i
            a b c
            a c b
            sage: f = Foliation(i, [1, 2, 3], 1); f
            a b c
            a c b
            Lengths: (1/6, 1/3, 1/2)
            Twist: 1/6
            sage: f._divpoints
            [[(0, (0, 0, 0, 0)), (1/6, (1, 0, 0, 0)), (1/2, (1, 1, 0, 0))], [(1/6, (0, 0, 0, 1)), (1/3, (1, 0, 0, 1)), (5/6, (1, 0, 1, 1))]]
            sage: f._divvalues
            [[0, 1/6, 1/2], [1/6, 1/3, 5/6]]
            sage: f._lengths['a']
            (1/6, (1, 0, 0, 0))
            sage: f._lengths['b']
            (1/3, (0, 1, 0, 0))
            sage: f._lengths['c']
            (1/2, (0, 0, 1, 0))
            sage: f._twist
            (1/6, (0, 0, 0, 1))
            sage: f._involution
            a b c
            a c b

            sage: i = Involution('a a b b c c'); i
            a a b b c c
            Moebius band
            sage: f = Foliation(i, [1, 2, 3]); f
            a a b b c c
            Moebius band
            Lengths: (1/12, 1/6, 1/4)
            sage: f._divpoints
            [[(0, (0, 0, 0, 0, 0)),
              (1/12, (1, 0, 0, 0, 0)),
                (1/6, (2, 0, 0, 0, 0)),
                  (1/3, (2, 1, 0, 0, 0)),
                    (1/2, (2, 2, 0, 0, 0)),
                      (3/4, (2, 2, 1, 0, 0))],
                       [(0, (0, 0, 0, 0, 1)), (1/2, (0, 0, 0, 1, 1))]]
            sage: f._divvalues
            [[0, 1/12, 1/6, 1/3, 1/2, 3/4], [0, 1/2]]
            sage: f._lengths['a']
            (1/12, (1, 0, 0, 0, 0))
            sage: f._lengths['b']
            (1/6, (0, 1, 0, 0, 0))
            sage: f._lengths['c']
            (1/4, (0, 0, 1, 0, 0))
            sage: f._twist
            (0, (0, 0, 0, 0, 1))
            sage: f._involution
            a a b b c c
            Moebius band

        """
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

    def __eq__(self, other):
        return self._involution == other._involution and\
                abs(self._length_twist_vector -
                    other._length_twist_vector) < epsilon

    def _repr_(self):
        if self._involution.is_bottom_side_moebius():
            return "{0}\nLengths: {1}".format(\
                    self._involution, self._length_twist_vector[:-2])
        return "{0}\nLengths: {1}\nTwist: {2}".format(\
                self._involution, self._length_twist_vector[:-1], 
                self._twist.value)

       

    def _latex_(self, size = 15, color_strength = 50, 
            interval_labelling = True,
            length_labelling = True):
        r"""
        Returns the Latex/Tikz representation of the foliation.

        INPUT:

        - ``size`` - the width of the Tikz picture in cm's. (Usually the height
          as well, though that can be smaller.

        - ``color_strength`` - an integer on a scale of 0 to 100. The lower
          the strength the more white is mixed with the colors. 0 means no
          coloring.

        - ``interval_labelling`` - boolean, whether the intervals are 
          labelled

        - ``length_labelling`` - boolean, whether the lengths of the intervals
          are indicated on the picture. For simpler Foliations this could be
          useful, but if there are really short intervals, it is better to 
          turn this off since there is no room for these decimals.

        OUTPUT:

        - string - the latex representation

        """
        latex.add_to_preamble('\usepackage{tikz}\n')
        s = '\\begin{{tikzpicture}}[scale = {0},'\
            'font=\\tiny]\n'.format(size)
        singularities = ''
        lines = ''
        fillings = ''
        labels = ''
        pos = 'above'
        for i in {0,1}:
            if i == 1:
                pos = 'below'
            if i == 1 and self._involution.is_bottom_side_moebius():
                lines += '\\draw (0,-0.2) [dotted] -- (1,-0.2);\n'
                fillings += '\\fill[yellow!{cs}!white] (0,0) rectangle '\
                        '(1,-0.2);\n'.format(cs = color_strength / 2)
                labels += '\\node[font = \large] at (0.5, -0.1) '\
                        '{{Moebius band}};\n'
                break
            for j in range(len(self._divvalues[i])):
                begin_percent = color_strength
                end_percent = 0

                signed_letter = letter = self._involution[i][j]
                if self._involution.is_flipped((i,j)):
                    signed_letter = '-' + letter
                    if (i, j) > self._involution.pair((i,j)):
                        begin_percent, end_percent = end_percent, begin_percent

                x1 = self._divvalues[i][j]
                x2 = self._divvalues[i][j] + \
                        self._lengths[letter].value
                midx = (x1 + x2)/2
                y = (-1)^i*0.5
                p1 = '({0},0)'.format(x1)
                p2 = '({0},{1})'.format(x1, y)
                p3 = '({0},{1})'.format(x2, y)
                p4 = '({0},0)'.format(x2)

                color = _tikzcolor(self._involution.index(letter))

                lines += '\\draw {0} -- {1};\n'.format(p1, p2)
                if i == 0 or j < len(self._divvalues[i]) - 1:
                    lines += '\\draw[dashed] {0} -- {1};\n'.format(p2, p3)
                    fillings += '\\shade[left color = {0}!{bp}!white, '\
                            'right color = {0}!{ep}!white] {1} rectangle '\
                            '{2};\n'.format(color, p1, p3, 
                                    bp = begin_percent, ep = end_percent)
                else:
                    p3 = '({0},{1})'.format(x2 - 1, y)
                    if x2 - 1 > 1 - x1:
                        midx -= 1
                    lines += '\\draw[dashed] {0} -- {1};\n'.format(p2, 
                        '(1,-0.5)')
                    lines += '\\draw[dashed] {0} -- {1};\n'.format(
                            '(0,-0.5)', p3)
                    middle_percent = color_strength * (x2 - 1) / (x2 - x1)
                    fillings += '\\shade[left color = {0}!{bp}!white, '\
                            'right color = {0}!{mp}!white] {1} rectangle '\
                            '{2};\n'.format(color, p1, '(1,-0.5)', 
                                    mp = middle_percent, bp = begin_percent)
                    fillings += '\\shade[left color = {0}!{mp}!white, '\
                            'right color = {0}!{ep}!white] {1} rectangle '\
                            '{2};\n'.format(color, '(0,0)', p3, 
                                    mp = middle_percent, ep = end_percent)

                sing_color = _tikzcolor(self._involution.\
                        which_singularity((i,j)))
                singularities += '\\filldraw[fill={col}, draw = black] {0} '\
                        'circle (0.005);\n'.format(p2, col = sing_color)
                if length_labelling:
                    labels += '\\node at ({0},0) [{1}] {{{2}}};\n'.\
                            format(midx, pos, signed_letter)
                if interval_labelling:
                    labels += '\\node at ({0},{1}) [{2}] {{{3}}};\n'.\
                            format(midx, y, pos, round(x2 - x1, 4))

        lines += '\\draw (1,0) -- (1,0.5);\n'\
                '\\draw[very thick] (0,0) -- (1,0);\n'
        s += fillings + lines + singularities + labels
        s += '\\end{tikzpicture}'
        return s


    @classmethod
    def RP2_arnoux_yoccoz(self):
        sf = arnoux_yoccoz_factor(3)
        return Foliation(Involution.RP2_arnoux_yoccoz(), 
                [1/sf + 1/sf^2, 1/sf^2 + 1/sf^3, 1/sf + 1/sf^3])
                

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

