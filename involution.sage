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






