from collections import deque, namedtuple
epsilon = 1e-10

def _mod_one(a):
    return a - floor(a)

def get_alphabet(alphabet_length):
    if alphabet_length <= 0:
        raise ValueError('The length of the alphabet must be '
                    'at least 1')
    if alphabet_length > 26:
        return range(1, alphabet_length + 1)
    else:
        return map(chr, range(97, 97 + alphabet_length))

class Involution(object):
    def __init__(self, pairing, flips = None):
        if isinstance(pairing, str):
            l = pairing.split()
        elif isinstance(pairing, list):
            l = pairing[:]
        letters = set(l)
        for letter in letters:
            if l.count(letter) != 2: 
                raise ValueError('Letters must reappear twice')

        if flips == None:
            flips = set()
        else:
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
        return self._pairing[index]

    def __eq__(self, other):
        if self._permutation != other._permutation:
            return False
        return map(self.index, self.flips()) == map(other.index,
                other.flips())

    @classmethod
    def random(cls, alphabet_length, without_flips = True):
        import random
        l = 2 * get_alphabet(alphabet_length)
        random.shuffle(l)
        if without_flips:
            flips = []
        else:
            k = random.randint(1, alphabet_length)
            flips = random.sample(set(l), k)
        return Involution(l, flips)

    @classmethod
    def random_for_pseudo_anosov(cls, alphabet_length,
            without_flips = True):
        if alphabet_length <= 2:
            raise ValueError('A foliation without saddle connnection'
                    ' needs an involution with an alphabet of length'
                    ' at least three')
        if alphabet_length == 3:
            if without_flips:
                return Involution('a a b b c c')
            else:
                raise ValueError('A non-orientable pseudo-anosov '
                        'foliation must come from an involution '
                        'with alphabet of length at least four.')
        while True:
            inv = Involution.random(alphabet_length, without_flips)
            s = inv.singularity_type()
            if 1 not in s and 2 not in s and\
                    not inv.has_saddle_connection(): 
                return inv 

    @classmethod
    def simplest_involution(cls, alphabet_length, without_flips=True):
        a = get_alphabet(alphabet_length)
        if without_flips:
            flips = []
        else:
            flips = a
        return Involution(sorted(2 * a), flips)

    def alphabet(self):
        return set(self._pairing)

    def flips(self):
        return set(self._flips)

    def index(self, letter):
        return self._index[letter]

    def to_permutation(self):
        return self._permutation

    def rotate(self, k):
        new_list = deque(self._pairing)
        new_list.rotate(k)
        return Involution(list(new_list), self._flips)

    def reverse(self):
        new_list = deque(self._pairing)
        new_list.reverse()
        return Involution(list(new_list), self._flips)

    def singularity_partition(self):
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
                    if self._pairing[j] not in self._flips:
                        j = (j + 1) % n
                        to_left = False
                else:
                    j = self._permutation[j] - 1
                    if self._pairing[j] not in self._flips:
                        to_left = True
                    else:
                        j = (j + 1) % n

                partition[-1].append(j)
                done.add(j)
                if j == i:
                    break
        return partition

    def singularity_type(self):
        return tuple(sorted(map(len, self.singularity_partition()), 
            reverse = True))

    def has_saddle_connection(self):
        n = len(self._pairing) / 2
        for i in range(n):
            if len(set(self._pairing[i:i+n])) == n:
                return True
        return False

class PointNonOr(object):
    def __init__(self, point, coefficients):
        self.point = point
        self.coefficients = vector(coefficients)

    def __repr__(self):
        return repr((self.point, self.coefficients))

    def __add__(self, other):
        return PointNonOr(self.point + other.point, 
                self.coefficients + other.coefficients)

    def __sub__(self, other):
        return PointNonOr(self.point - other.point, 
                self.coefficients - other.coefficients)

    def adjusted(self):
        if self.point == 0: #need to check beacuse of a bug
            # with algebraic numbers
            n = 0
        else:
            n = floor(self.point)
        return PointNonOr(self.point - n, [x - 2 * n 
            for x in self.coefficients])

    def half_added(self):
        new_point = PointNonOr(self.point + 1/2,
                [x + 1 for x in self.coefficients])
        return new_point.adjusted()

def basis_vector(n, k):
    l = [0] * n
    l[k] = 1
    return l

def arnoux_yoccoz_factor(genus, field = RDF):
    R.<x> = ZZ[]
    poly = R([-1] * genus + [1])
    return max([abs(x[0]) for x in poly.roots(field)])

def _less(x, y, is_equality_allowed):
        if is_equality_allowed:
            return x <= y
        else:
            return x < y

class Interval(object):
    def __init__(self, left_endpoint, right_endpoint,
            is_closed_on_left = True, is_closed_on_right = False):
        self._lep = left_endpoint
        self._rep = right_endpoint
        self._is_closed_on_left = is_closed_on_left
        self._is_closed_on_right = is_closed_on_right

    def __repr__(self):
        s = ''
        if self._is_closed_on_left:
            s += '['
        else:
            s += '('
        s += str(self._lep) + ',' + str(self._rep)
        if self._is_closed_on_right:
            s += ']'
        else:
            s += ')'
        return s

    def contains(self, point):
        if self._lep <= self._rep:
            if _less(self._lep, point, self._is_closed_on_left) and\
                    _less(point, self._rep, self._is_closed_on_right):
                return True
            else:
                return False
        if _less(self._rep, point, not self._is_closed_on_right) and\
                _less(point, self._lep, not self._is_closed_on_left):
            return False
        else:
            return True




def is_positive(vec):
    if abs(vec[0]) < epsilon:
        return False
    is_first_positive = (vec[0] > 0)
    for x in vec:
        if abs(x) < epsilon  or (x > 0) != is_first_positive:
            return False
    return True

def get_good_eigendata(transition_matrix, is_surface_orientable):
    ev = transition_matrix.eigenvectors_right()
    ret_list = []
    for x in ev:
        if abs(x[0].imag()) < epsilon and x[0] > 0 \
                and abs(x[0] - 1) > epsilon:
            for vec in x[1]:
                if is_positive(vec):
                    norm = sum([abs(y) for y in vec])
                    if is_surface_orientable:
                        norm -= abs(vec[-1])
                    else:
                        norm *= 2
                    normalized_vec = [abs(y / norm) for y in vec]
                    ret_list.append((x[0], normalized_vec))
    return ret_list
 

class SaddleConnectionError(Exception):
    pass

class RestrictionError(Exception):
    def __init__(self, value):
        self.value = value


class FoliationNonOrientableSurface(object):
    def __init__(self, involution, lengths):
        sing_part = involution.singularity_partition()
        self._singularity_type = \
                tuple(sorted(map(len, sing_part), reverse = True))
        if 1 in self._singularity_type:
            raise ValueError('The foliation has a one-pronged '
                    'singularity')
        if 2 in self._singularity_type:
            raise ValueError('The same interval exchange can be '
                    'specified on a smaller number of intervals.')
        self._involution = involution
        if isinstance(lengths, list):
            if len(lengths) != len(involution.alphabet()):
                raise TypeError('Bad number of lengths')
            l = {}
            done = set()
            count = 0
            for x in involution:
                if x not in done:
                    done.add(x)
                    l[x] = lengths[count]
                    count += 1
            lengths = l

        if isinstance(lengths, dict):
            if involution.alphabet() != set(lengths.keys()):
                raise TypeError('Invalid length specification')

        for v in lengths.values():
            if v <= 0:
                raise ValueError('Lengths must be positive')

        total = sum(lengths.values()) 
        self._lengths = {}
        n = len(lengths)
        vec_list = [0] * n
        for k in lengths:
            self._lengths[k] = PointNonOr(lengths[k]/2/total,\
                    basis_vector(n, involution.index(k)))
            vec_list[self._involution.index(k)] =\
                    self._lengths[k].point

        self._length_vector = vector(vec_list)

        self._divpoints = [PointNonOr(0, [0] * n)]
        for i in range(2 * n):
            self._divpoints.append(self._divpoints[-1] + \
                    self._lengths[involution[i]])

        self._divpoints_justpoints =[x.point for x in self._divpoints]

    def __eq__(self, other):
        return self.involution() == other.involution() and \
                abs(self._length_vector - other._length_vector) <\
                epsilon

    def __repr__(self):
        return repr(self._involution) + '\n' +\
                repr([self._lengths[x].point for x in self._involution])

    @classmethod
    def arnoux_yoccoz_foliation(cls, genus, field = RDF):
        inv = Involution.simplest_involution(genus - 1)
        sf = arnoux_yoccoz_factor(genus - 1, field)
        lengths = [1/sf**i for i in range(genus - 1)]
        return FoliationNonOrientableSurface(inv, lengths)

    @classmethod
    def random(cls, alphabet_length, involution = None, 
            orientable = True):
        import random
        if involution == None:
            involution = Involution.random_for_pseudo_anosov(\
                    alphabet_length,
                    without_flips = orientable)
        lengths = [random.random() for i in range(alphabet_length)]
        return FoliationNonOrientableSurface(involution, lengths)

    def alphabet_length(self):
        return len(self._lengths) 

    def singularity_type(self):
        return self._singularity_type

    def euler_char(self):
        return sum([2 - k for k in self._singularity_type]) / 2

    def genus(self):
        return 2 - self.euler_char()

    def is_orientable(self):
        return len(self._involution.flips()) == 0

    def involution(self):
        return self._involution

    def permutation(self):
        return self._involution.to_permutation()

    TransitionData = namedtuple('TransitionData', 'tr_matrix,new_inv')

    def new_foliation(self, transition_data):
        new_vector = transition_data.tr_matrix * self._length_vector 
        return FoliationNonOrientableSurface(transition_data.new_inv, new_vector.list())


    def _simple_transformation_data(self, new_involution):
        m = matrix(self.alphabet_length())
        for letter in self._involution.alphabet():
            m[new_involution.index(letter), 
                    self._involution.index(letter)] = 1
        return self.TransitionData(tr_matrix = m, 
                new_inv = new_involution)


    def _rotation_data(self, k):
        return self._simple_transformation_data(\
                self._involution.rotate(k))

    def _reverse_data(self):
        return self._simple_transformation_data(\
                self._involution.reverse())

    def _in_which_interval(self, point):
        import bisect
        interval =bisect.bisect(self._divpoints_justpoints, point) - 1
        if interval == len(self._divpoints_justpoints) - 1 or \
                True in [abs(point - 
                    self._divpoints_justpoints[i]) < epsilon
                    for i in {interval, interval + 1}]:
            raise ValueError('Containing interval cannot be found '
                    'because the point is too close too division '
                    'point.')
        return interval

    def _map(self, point, containing_interval):
        flipped = (self._involution[containing_interval] 
                in self._involution.flips())
        new_interval = self._involution.to_permutation()[\
                containing_interval] - 1
        diff = point - self._divpoints[containing_interval]
        if not flipped:
            return self._divpoints[new_interval] + diff
        else:
            return self._divpoints[new_interval + 1] - diff
        
    IntersectionData = namedtuple('IntersectionData',
            'point, is_from_above, orientation_reversing')

    def _first_intersection(self, interval_list,index_of_singularity):
        point = self._divpoints[index_of_singularity]
        count = 0
        
        def is_contained(p):
            for interval in interval_list:
                if interval.contains(p.point):
                    return True
            return False

        from_above = True
        reversing = False
        while not is_contained(point):
            count += 1
            if count > 10000:
                print self, interval_list, index_of_singularity, point
            point = point.half_added()

            if is_contained(point): 
                from_above = False
                break
            try:
                interval = self._in_which_interval(point.point)
            except ValueError:
                raise SaddleConnectionError()
            point = self._map(point, interval)
            if self._involution[interval] in self._involution.flips():
                reversing = not reversing

        return self.IntersectionData(point, from_above, reversing)


    def _create_transition_data(self, distances):
        tuples = list(enumerate(distances))
        tuples.sort(key = lambda x: x[1][0].point)
        n = len(distances)
        
        if self.is_orientable():
            new_involution = Involution([self._involution[x[0]]
                for x in tuples])
        else:
            done = set() 
            flips = set()
            remaining_letters = self._involution.alphabet()
            inv_list = [0] * n
            for i in range(n):
                if i in done:
                    continue
                if len(remaining_letters) == 0:
                    print self
                    print distances
                    print tuples
                letter = remaining_letters.pop()
                right_side = True
                k = tuples[i][0]
                if tuples[i][1][1]: #the orientation is reversed
                    right_side = False
                    k = (k - 1) % n
                k = self.permutation()[k] - 1
                if self._involution[k] in self._involution.flips():
                    right_side = not right_side
                if not right_side:
                    k = (k + 1) % n
                k = next(i for i in range(n) if tuples[i][0] == k)
                if tuples[k][1][1]:
                    right_side = not right_side
                if not right_side:
                    k = (k - 1) % n
                inv_list[i] = inv_list[k] = letter
                done.add(i)
                done.add(k)
                if not right_side:
                    flips.add(letter)
            new_involution = Involution(inv_list, flips)
        
        m = matrix(self.alphabet_length())
        done = set()
        for i in range(2*self.alphabet_length()):
            letter = new_involution[i]
            if letter in done:
                continue
            done.add(letter)
            m[len(done) - 1] =\
                    tuples[i + 1][1][0].coefficients - \
                    tuples[i][1][0].coefficients
 
        return self.TransitionData(tr_matrix = m,
                new_inv = new_involution)

    def _restrict_to_not_flipped(self, interval):
        if self._involution[interval] in self._involution.flips():
            raise RestrictionError('The specified interval should '
                    'not be flipped, but it is.')
        n = 2 * self.alphabet_length()
        left_endpoint = self._divpoints[(interval + 1) % n]
        right_endpoint = self._divpoints[\
                self._involution.to_permutation()[interval] % n]

        intersections = [self._first_intersection(\
                [Interval(left_endpoint.point, 
                    right_endpoint.point)], i) 
            for i in range(n)]
            

        total = (right_endpoint - left_endpoint).adjusted()

        def distance_from_left_endpoint(intersection_data):
            diff =(intersection_data.point - left_endpoint).adjusted()
            if intersection_data.is_from_above:
                return diff
            else:
                return total + diff
        
        distances = [(distance_from_left_endpoint(x), 
            x.orientation_reversing) for x in intersections]

        return self._create_transition_data(distances)



    def _restrict_to_flipped(self, center_left_index, 
            center_right_index):
        if not self._involution.flips().issuperset(\
                {self._involution[center_left_index],
                    self._involution[center_right_index]}):
                    raise RestrictionError('The specified interval '
                            'should be flipped, but it isn\'t.')

        n = len(self._divpoints) - 1
        diff = (center_right_index - center_left_index) % n
        diff2 = (self.permutation()[center_right_index] - 1 -
                center_left_index) % n
        if diff2 < diff or diff == 0 or diff2 == 0:
            raise RestrictionError('Specified transverse curve does not '
                    'exist. Combinatorically invalid choice for '
                    'the intervals.')

        other_left = self._divpoints[self.permutation()[\
                center_left_index] - 1]
        other_right = self._divpoints[self.permutation()[\
                center_right_index] % n]
        if (other_right - other_left).adjusted().point <= 0.5:
            raise RestrictionError('Specified transverse curve does not '
                    'exist. The length parameters doesn\'t allow '
                    'closing the curve to a transverse curve.')
        other_left = other_left.half_added()
        
        center_left = self._divpoints[(center_left_index + 1) % n]
        center_right = self._divpoints[center_right_index]
        
        center_int = Interval(center_left.point, center_right.point,
                is_closed_on_left = True, is_closed_on_right = True)
        
        other_int = Interval(other_left.point, other_right.point, 
                is_closed_on_left = False,
                is_closed_on_right = False)

        intersections = [self._first_intersection(\
                [center_int, other_int], i)
                for i in range(2 * self.alphabet_length())]

        total = (center_right - center_left).adjusted() +\
                (other_right - other_left).adjusted()

        def distance_from_left_endpoint(intersection_data):
            p = intersection_data.point
            or_ = intersection_data.orientation_reversing
            if intersection_data.is_from_above:
                if center_int.contains(p.point):
                    return ((p - center_left).adjusted(), or_)
                return (total + total - (p - other_left).adjusted(),
                        not or_)
            if other_int.contains(p.point):
                return (total - (p - other_left).adjusted(), not or_)
            return (total + (p - center_left).adjusted(), or_)

        distances = [distance_from_left_endpoint(x) 
            for x in intersections]

        return self._create_transition_data(distances)


    def _find_pseudo_anosov_candidates(self, depth,
            transition_matrix_so_far = None,
            initial_foliation = None,
            encoding_sequence_so_far = None):
        if transition_matrix_so_far == None:
            transition_matrix_so_far = matrix.identity(
                    self.alphabet_length())
            initial_foliation = self
            encoding_sequence_so_far = []
        #print self
        #print encoding_sequence_so_far

        ret_list = []

        n = 2 * self.alphabet_length()

        for i in range(n):
            tr_data = self._rotation_data(i)
            final_foliation = self.new_foliation(tr_data)
            final_matrix = tr_data.tr_matrix *transition_matrix_so_far
            is_rev = False
            for j in range(2):
                if j == 1:
                    tr_data = final_foliation._reverse_data()
                    final_matrix = tr_data.tr_matrix * final_matrix
                    is_rev = True
                #print final_matrix, final_matrix.eigenvectors_right()

                if tr_data.new_inv == initial_foliation.involution():
                    for eigen_data in get_good_eigendata(\
                            matrix(RDF, final_matrix),
                            is_surface_orientable = False):
                        fol = FoliationNonOrientableSurface(\
                                initial_foliation.involution(),
                                eigen_data[1])
                        ret_list.append(PseudoAnosov(\
                                foliation = fol,
                                encoding_seq = \
                                encoding_sequence_so_far,
                                rotatedby = i,
                                is_reversed = is_rev,
                                tr_matrix = final_matrix,
                                good_eigen_data = eigen_data))


        if depth > 0:
            for i in range(n):
                if self._involution[i] in self._involution.flips():
                    done = set()
                    j = (i + 1) % n 
                    while self._involution[j] != self._involution[i]:
                        if self._involution[j] in done or\
                                self._involution[j] not in\
                                self._involution.flips():
                            j = (j + 1) % n 
                            continue
                        done.add(self._involution[j])
                        try:
                            tr_data = self._restrict_to_flipped(i, j)
                        except (SaddleConnectionError,
                                RestrictionError):
                            j = (j + 1) % n 
                            continue
                        ret_list += self.new_foliation(tr_data).\
                                _find_pseudo_anosov_candidates(\
                                depth - 1,
                                tr_data.tr_matrix * 
                                transition_matrix_so_far,
                                initial_foliation,
                                encoding_sequence_so_far + [(i,j)])
                else:
                    try:
                        tr_data = self._restrict_to_not_flipped(i)
                    except SaddleConnectionError:
                        continue
                    ret_list += self.new_foliation(tr_data).\
                                _find_pseudo_anosov_candidates(\
                                depth - 1,
                                tr_data.tr_matrix * 
                                transition_matrix_so_far,
                                initial_foliation,
                                encoding_sequence_so_far + [i])

        return ret_list
                           
    def find_pseudo_anosovs(self, depth):
        candidates = self._find_pseudo_anosov_candidates(depth)
        return filter(PseudoAnosov.verify, candidates)


class PseudoAnosov(namedtuple('PseudoAnosov', 'foliation, encoding_seq,'
    'rotatedby, is_reversed, tr_matrix, good_eigen_data')):
    def __repr__(self):
        return "Pseudo-anosov on the genus {0} non-orientable surface "\
               "with stretch factor {1} and polynomial {2}".format(\
               self.foliation.genus(), self.stretch_factor(),
               self.polynomial())



    def stretch_factor(self):
        s = self.good_eigen_data[0]
        if s > 1:
            return s
        else:
            return 1/s

    def polynomial(self):
        return self.tr_matrix.charpoly()

    def is_exact(self):
        return self.good_eigen_data[0] in QQbar

    def has_totally_real_trace_field(self):
        if self.is_exact():
            sf = self.stretch_factor()
        else:
            for ev in self.tr_matrix.eigenvalues():
                if abs(ev - self.good_eigen_data[0]) < epsilon:
                    sf = ev
                    break
        trace = sf + 1/sf
        p = trace.minpoly()
        nroots = sum([root[1] for root in p.roots(RDF)])
        return nroots == p.degree()

    def get_exact_pseudo_anosov(self):
        if self.is_exact():
            return self
        for x in get_good_eigendata(self.tr_matrix,
                is_surface_orientable = False):
            if abs(vector(x[1]) - self.foliation._length_vector) < epsilon:
                f = FoliationNonOrientableSurface(\
                        self.foliation.involution(),
                        x[1])
                return self._replace(foliation = f,
                        good_eigen_data = x)
        raise RuntimeError('Something bad is going on with '
                'exact and approximate eigenvector calculations.')

    def verify(self):
        fol = self.foliation
        tr_matrix = matrix.identity(self.tr_matrix.ncols())

        for x in self.encoding_seq:
            if isinstance(x, tuple):
                try:
                    tr_data = fol._restrict_to_flipped(*x)
                except (SaddleConnectionError, 
                        RestrictionError):
                    return False
            else:
                try:
                    tr_data = fol._restrict_to_not_flipped(x)
                except (SaddleConnectionError,
                        RestrictionError):
                    return False
            tr_matrix = tr_data.tr_matrix * tr_matrix
            try:
                fol = fol.new_foliation(tr_data)
            except ValueError:
                return False

        tr_data = fol._rotation_data(self.rotatedby)
        tr_matrix = tr_data.tr_matrix * tr_matrix
        fol = fol.new_foliation(tr_data)

        if self.is_reversed:
            tr_data = fol._reverse_data()
            tr_matrix = tr_data.tr_matrix * tr_matrix
            fol = fol.new_foliation(tr_data)

        return tr_matrix == self.tr_matrix

    @classmethod
    def find_stretch_factors(cls, alphabet_length, depth = 2, 
            timelimit = 10, is_orientable = True):
        import time
        start = time.time()
        pas = []
        while time.time() < start + timelimit:
            pas = FoliationNonOrientableSurface.random(\
                    alphabet_length, orientable = is_orientable)\
                    .find_pseudo_anosovs(depth) 
        tuples = set([(pa.foliation.genus(), pa.stretch_factor(),
            pa.polynomial()) for pa in pas])
        return sorted(tuples)



class PseudoAnosovDatabase(object):
    def __init__(self):
        self._depth_of_search = {} 

    def _pick_a_pseuso_anosov(self):
        return min(self._depth_of_search,
                key = self._depth_of_search.get)

    def start_search(self):
        do_random_search = False
        if len(self._entries) == 0:
            do_random_search = True
        else:
            pick = self._pick_a_pseuso_anosov()
            if self._depth_of_search[pick] >= 2:
                do_random_search = True





h = FoliationNonOrientableSurface.arnoux_yoccoz_foliation(4)
h8 = FoliationNonOrientableSurface(Involution('a b a c b c d d'),[0.0944422795324, 0.178500974036, 0.074245928173, 0.152810818258])
h8_2 = FoliationNonOrientableSurface(Involution('a b a c b c d d'), [0.0820003383086, 0.194387134927, 0.065494592731, 0.158117934033])    

f = FoliationNonOrientableSurface(Involution('c c b d b a d a', ['b', 'a']), [0.0962912017836, 0.14222527893, 0.119258240357, 0.14222527893]) 
g = FoliationNonOrientableSurface(Involution('a a d c d b c b', ['b', 'd', 'c']), [0.0944422795324, 0.00178500974036, 0.0074245928173, 0.00152810818258]) 
g._restrict_to_flipped(3,4)
