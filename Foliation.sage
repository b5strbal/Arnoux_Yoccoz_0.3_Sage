from sage.structure.sage_object import SageObject
from bisect import bisect
from collections import namedtuple

epsilon = 1e-10

_global_storage = []

def _integer_mod(a, b):
    x = a % b
    if x < 0:
        return x + b
    else:
        return x

def _mod_one(a):
    return a - floor(a)


def rotating_permutation(n, k):
    return Permutation([(i + k) % n + 1 for i in range(n)])

def rotating_matrix(n, k):
    return rotating_permutation(n, k).to_matrix()

def is_minimal(permutation):
    n = permutation.size()
    if n < 3:
        return True
    for i in range(n):
        if permutation[(i + 1) % n] - permutation[i] in [1, 1 - n]:
            return False
    return True


def singularity_permutation(permutation):
    p = permutation
    p_inv = p.inverse()
    n = p.size()
    return Permutation([p_inv[p[(i - 1) % n] % n] for i in range(n)])


def _canonical_form(intex, twist):
    one_minus_twist = _mod_one(1 - twist)
    rotateby = bisect(intex.range_singularities(), one_minus_twist)
    #for i in [rotateby - 1, rotateby]:
    #    if abs(one_minus_twist - intex.range_singularities()[i]) < 0.000000001:
    #        raise ValueError, 'The foliation has a closed leaf. It should not.'
    newtwist = intex.range_singularities()[rotateby] - one_minus_twist 
    n = len(intex.lengths())
    return (iet.IntervalExchangeTransformation(rotating_permutation(n, rotateby) *\
            intex.permutation().to_permutation(), intex.lengths()), newtwist)

def is_between(left_endpoint, right_endpoint, point):
    if left_endpoint < right_endpoint:
        if point <= left_endpoint or point > right_endpoint:
            return False
        else:
            return True
    else:
        if point <= right_endpoint or point > left_endpoint:
            return True
        else:
            return False













class PointWithCoordinates(SageObject):
    def __init__(self, point, coordinates):
        self.point = point
        self.coordinates = coordinates

    def __repr__(self):
        return repr((self.point, self.coordinates))

    def __add__(self, other):
        new_point = PointWithCoordinates(self.point + other.point, 
                self.coordinates + other.coordinates)
        if new_point.point >= 1:
            new_point.point -= 1
            for i in range(len(self.coordinates) - 1):
                new_point.coordinates[i] -= 1
        return new_point

    def __sub__(self, other):
        new_point = PointWithCoordinates(self.point - other.point, 
                self.coordinates - other.coordinates)
        if new_point.point < 0:
            new_point.point += 1
            for i in range(len(self.coordinates) - 1):
                new_point.coordinates[i] += 1
        return new_point       

 
def is_positive(vec):
    if abs(vec[0]) < epsilon:
        return False
    is_first_positive = (vec[0] > 0)
    for x in vec:
        if abs(x) < epsilon  or (x > 0) != is_first_positive:
            return False
    return True
        
def get_pseudo_anosov_candidates(transition_matrix, is_orientable):
    global _global_storage
    ev = transition_matrix.eigenvectors_right()
    ret_list = []
    for x in ev:
        if x[0].imag() == 0.0 and x[0] > 0 and abs(x[0] - 1.0) > epsilon:
            for vec in x[1]:
                if is_positive(vec):
                    norm = sum([abs(y) for y in vec])
                    if is_orientable:
                        norm -= abs(vec[-1])
                    else:
                        norm *= 2
                    normalized_vec = [abs(y / norm) for y in vec]
                    ret_list.append((x[0], normalized_vec))
                else:
                    _global_storage.append(('no positive eigenvector:', vec))
        else:
            _global_storage.append(('no good eigenvalue:', x[0]))
    return ret_list
            



class Foliation(SageObject):
    r"""
    Oriented projective measured foliation on an orientable surface
    """
    def __init__(self, permutation, lengths, twist):
        if not is_minimal(permutation):
            raise ValueError, 'The same interval exchange can be specified on less intervals.'
        intex = iet.IntervalExchangeTransformation(permutation, lengths)
        normalized_twist = twist / intex.length()
        intex = intex.normalize()
        ret = _canonical_form(intex, twist)
        self._intex = ret[0]
        self._twist = ret[1]

        all_singularities = sorted(self._intex.domain_singularities() \
                + [x + self._twist for x in self._intex.range_singularities()])
        for i in range(2 * self.num_intervals()):
            if abs(all_singularities[i + 1] - all_singularities[i]) < epsilon:
                raise ValueError, ("The foliation and immediate saddle connection "
                        "or one of the intervals of the exchange is very short.")
        self._singularity_cycles = singularity_permutation(self._intex.permutation().to_permutation()).to_cycles()
        self._translation_matrix = matrix(self.num_intervals(), self.num_intervals() + 1)
        inverse_perm = self.permutation().inverse()
        for i in range(self.num_intervals()):
            for j in range(self.num_intervals()):
                if i < j and inverse_perm[i] > inverse_perm[j]:
                    self._translation_matrix[i, j] = 1
                if i > j and inverse_perm[i] < inverse_perm[j]:
                    self._translation_matrix[i, j] = -1
            self._translation_matrix[i, -1] = 1
        self._translations = self._translation_matrix * vector(self._intex.lengths() + [self._twist])

    def __repr__(self):
        return "Oriented foliation on a genus {0} surface.\n".format(self.genus()) +\
            "Permutation: {0}\n".format(self.permutation()) +\
            "Lengths: {0}\n".format(self._intex.lengths()) +\
            "Twist: {0}".format(self._twist)

    def euler_char(self):
        return -sum([len(x) - 1 for x in self._singularity_cycles])

    def singularity_type_of_Abelian_diff(self):
        return sorted([len(x) - 1 for x in self._singularity_cycles if len(x) > 1])

    def singularity_type_prongs(self):
        return [(x + 1) * 2 for x in self.singularity_type_of_Abelian_diff()]

    def genus(self):
        return (2 - self.euler_char()) / 2
     
    def num_intervals(self):
        return len(self._intex.lengths())

    def permutation(self):
        return self._intex.permutation().to_permutation()

    def rotation_data(self, k):
        normalized_k = _integer_mod(k, self.num_intervals())
        rotating_dist = sum(self._intex.lengths()[self.num_intervals() - normalized_k:])
        containing_interval = bisect(self._intex.range_singularities(), _mod_one(-rotating_dist- self._twist)) 
        new_permutation = rotating_permutation(self.num_intervals(), containing_interval) *\
                self.permutation() * \
                rotating_permutation(self.num_intervals(), k)
        m = matrix(self.num_intervals() + 1)
        for i in range(self.num_intervals()):
            m[i, _integer_mod(i - k, self.num_intervals())] = 1
        for j in range(normalized_k):
            m[-1, self.num_intervals() - j - 1] += 1
        for j in range(self.num_intervals() - containing_interval):
            m[-1, self.permutation()[self.num_intervals() - j - 1] - 1] -= 1
        m[-1, -1] = 1
        return (m, new_permutation)
        
    def new_foliation(self, transition_matrix, new_permutation):
        old_vector = vector(self._intex.lengths() + [self._twist])
        new_vector = transition_matrix * old_vector
        return Foliation(new_permutation, new_vector.list()[:-1], new_vector[-1])


    def first_intersection(self, left_endpoint, right_endpoint, index_of_singularity, upwards = True):
        point = self._intex.domain_singularities()[index_of_singularity]
        interval_intersection_count = [0] * self.num_intervals()
        if upwards:
            while True:
                i = bisect(self._intex.domain_singularities(), point) - 1
                point = _mod_one(point + self._translations[i])
                interval_intersection_count[i] += 1
                if is_between(left_endpoint, right_endpoint, point):
                    break
        else:
            while not is_between(left_endpoint, right_endpoint, point):
                i = bisect(self._intex.range_singularities(), _mod_one(point - self._twist)) - 1
                domain_interval = self.permutation()[i] - 1
                point = _mod_one(point - self._translations[domain_interval])
                interval_intersection_count[domain_interval] -= 1
        coefficients = vector([1] * index_of_singularity + [0] * (self.num_intervals() - index_of_singularity + 1))
        coefficients += vector(interval_intersection_count) * self._translation_matrix
        calculated_point = coefficients * vector(self._intex.lengths() + [self._twist])
        coefficients -= floor(calculated_point) * vector([1] * self.num_intervals() + [0])
        return PointWithCoordinates(point, coefficients) 

    def restrict_data(self, singularity, wrapping):
        if wrapping == 0:
            raise ValueError, 'wrapping cannot be 0'
        left_endpoint_coeff = vector([1] * singularity + [0] * (self.num_intervals() - singularity + 1))
        right_endpoint_coeff = left_endpoint_coeff + self._translation_matrix[singularity]
        left_endpoint = PointWithCoordinates(self._intex.domain_singularities()[singularity], left_endpoint_coeff)
        right_endpoint = PointWithCoordinates(self._intex.range_singularities()[self.permutation().inverse()[singularity] - 1]\
                + self._twist, right_endpoint_coeff)
        if wrapping < 0:
            left_endpoint, right_endpoint = right_endpoint, left_endpoint
        if wrapping in [-1, 1]:
            le, re = left_endpoint.point, right_endpoint.point
        else:
            le, re = -1, 1
        upper_intersections = [self.first_intersection(le, re, i, False) \
                for i in range(self.num_intervals())]
        lower_intersections = [self.first_intersection(le, re, i, True) \
                for i in range(self.num_intervals())]
        upper_tuples = list(enumerate(upper_intersections))
        if wrapping > 0:
            sort_key = lambda x: -_mod_one(left_endpoint.point - x[1].point)
        else:
            sort_key = lambda x: -_mod_one(right_endpoint.point - x[1].point)
        upper_tuples.sort(key = sort_key) 
     
        total_coordinates = right_endpoint.coordinates - left_endpoint.coordinates +\
                vector([1] * self.num_intervals() + [0]) * (abs(wrapping) - 1)
        if right_endpoint.point < left_endpoint.point:
            total_coordinates += vector([1] * self.num_intervals() + [0]) 
        total_length = total_coordinates * (vector(self._intex.lengths() + [self._twist]))
        
        transition_matrix = matrix(self.num_intervals() + 1)
        for i in range(self.num_intervals() - 1):
            transition_matrix[i] = (upper_tuples[i + 1][1] - upper_tuples[i][1]).coordinates
        transition_matrix[self.num_intervals() - 1] = total_coordinates -\
                (upper_tuples[-1][1] - upper_tuples[0][1]).coordinates

        def distance_from_left_endpoint(point_with_coordinates, is_close_to_left):
            if is_close_to_left:
                diff = point_with_coordinates - left_endpoint
                if diff.point == 0.0:
                    return (1.0, vector([1] * self.num_intervals() + [0]))
                else:
                    return (diff.point, diff.coordinates)
            else:
                diff = right_endpoint - point_with_coordinates
                return (total_length - diff.point, total_coordinates - diff.coordinates)

        reference_point_distance = distance_from_left_endpoint(upper_tuples[0][1], wrapping > 0)

        lower_distances_from_left_endpoint = sorted(list(enumerate([distance_from_left_endpoint(x, wrapping < 0) 
            for x in lower_intersections])), key = lambda x: x[1][0])

        position = bisect([x[1][0] for x in lower_distances_from_left_endpoint], reference_point_distance[0])

        transition_matrix[-1] = lower_distances_from_left_endpoint[position % self.num_intervals()][1][1] -\
            reference_point_distance[1]
        if position == self.num_intervals():
            transition_matrix[-1] += total_coordinates

        pairing = [0] * self.num_intervals()
        for i in range(self.num_intervals()):
            pairing[upper_tuples[i][0]] = i 
        perm = [pairing[lower_distances_from_left_endpoint[(position + i) % self.num_intervals()][0]] + 1
                for i in range(self.num_intervals())]

        return (transition_matrix, Permutation(perm))



class PointWithCoordinatesNonOr(SageObject):
    def __init__(self, point, coordinates):
        self.point = point
        self.coordinates = coordinates

    def __repr__(self):
        return repr((self.point, self.coordinates))

    def __add__(self, other):
        new_point = PointWithCoordinatesNonOr(self.point + other.point, 
                self.coordinates + other.coordinates)
        if new_point.point >= 1:
            new_point.point -= 1
            for i in range(len(self.coordinates)):
                new_point.coordinates[i] -= 2
        return new_point

    def __sub__(self, other):
        new_point = PointWithCoordinates(self.point - other.point, 
                self.coordinates - other.coordinates)
        if new_point.point < 0:
            new_point.point += 1
            for i in range(len(self.coordinates)):
                new_point.coordinates[i] += 2
        return new_point       

def _random_involution(num_intervals):
    import random
    if num_intervals % 2 != 0:
        raise ValueError, 'The number of intervals must be even.'
    if num_intervals < 6:
        raise ValueError, ('The involution corresponding to a foliation without saddle connection '
    'must act on at least 6 intervals.')
    if num_intervals == 6:
        return Permutation([2, 1, 4, 3, 6, 5])
    mix = range(num_intervals)
    perm_list = [0] * num_intervals
    while True:
        random.shuffle(mix)
        for i in range(num_intervals // 2):
            perm_list[mix[2 * i]] = mix[2 * i + 1] + 1
            perm_list[mix[2 * i + 1]] = mix[2 * i] + 1
        perm = Permutation(perm_list)
        if _is_involution_good(perm) and is_minimal(perm):
            return perm





class FoliationNonOrientableSurface(SageObject):
    r"""
    Normalized measured foliation with \mathbb{Z}_2 holonomy on a non-orientable surface.
    """

    def __init__(self, permutation, lengths):
        if permutation.number_of_fixed_points() > 0:
            raise ValueError, "The permutation should not have fixed points."
        if permutation.inverse() != permutation:
            raise ValueError, 'The permutation must be an involution.'
        if permutation.size() != 2 * len(lengths):
            raise ValueError, 'The list of lengths should be half as long as the size of the permutation'

        count = 0
        self._small_index = 2 * lengths
        for x in range(2 * len(lengths)):
            if permutation[x] - 1 > x:
                self._small_index[x] = self._small_index[permutation[x] - 1] = count
                count += 1
        
        total = sum(lengths)
        self._lengths = [x/total/2 for x in lengths]

        detailed_lengths = [self._lengths[self._small_index[i]] for i in range(2 * len(lengths))]
        self._intex = iet.IntervalExchangeTransformation(permutation, detailed_lengths)
        self._lifted_foliation = self.to_foliation()

        def get_next(prev, distance, index_to_increase):
            new = PointWithCoordinatesNonOr(prev.point + distance, prev.coordinates + vector([0] * len(lengths)))
            new.coordinates[index_to_increase] += 1
            return new


        self._divpoints = [PointWithCoordinatesNonOr(0, vector([0] * len(lengths)))] 
        for i in range(2 * len(lengths) - 1):
            self._divpoints.append(get_next(self._divpoints[-1], detailed_lengths[i], self._small_index[i]))

        self._translation_matrix = matrix(2 * len(lengths), len(lengths))
        inverse_perm = self.permutation().inverse()
        for i in range(2 * len(lengths)):
            for j in range(2 * len(lengths)):
                if i < j and inverse_perm[i] > inverse_perm[j]:
                    self._translation_matrix[i, self._small_index[j]] += 1
                if i > j and inverse_perm[i] < inverse_perm[j]:
                    self._translation_matrix[i, self._small_index[j]] += -1
        self._translations = self._translation_matrix * vector(self._lengths)

        self._index_of_interval = []
        count = 0
        for i in range(self.num_intervals()):
            j = self.permutation()[i] - 1 
            if j > i:
                self._index_of_interval.append(count)
                count += 1
            else:
                self._index_of_interval.append(self._index_of_interval[j])

    @classmethod
    def random(cls, num_intervals, permutation = None):
        import random
        if permutation == None:
            permutation = _random_involution(num_intervals)
        lengths = [random.random() for i in range(num_intervals // 2)]
        return FoliationNonOrientableSurface(permutation, lengths)


    def __repr__(self):
        return 'Foliation with Z_2 holomony on a non-orientable surface of genus {0}.\n'\
                'Permutation: {1}\nLengths: {2}'.format(self.genus(), 
                        self._intex.permutation().to_permutation(), self._intex.lengths())

    def num_intervals(self):
        return len(self._intex.lengths())
    
    def permutation(self):
        return self._intex.permutation().to_permutation()

    def to_foliation(self):
        return Foliation(self._intex.permutation().to_permutation(), self._intex.lengths(), 1/2)

    def euler_char(self):
        return self.to_foliation().euler_char() / 2

    def genus(self):
        return 2 - self.euler_char()

    def singularity_type_prongs(self):
        prongs_of_lift = self._lifted_foliation.singularity_type_prongs()
        return [prongs_of_lift[i] for i in range(len(prongs_of_lift)) if i % 2 == 0] 

    TransitionData = namedtuple('TransitionData', 'tr_matrix, new_perm')

    def rotation_data(self, k):
        new_permutation = rotating_permutation(self.num_intervals(), -k) *\
                self.permutation() *\
                rotating_permutation(self.num_intervals(), k)
        m = matrix(self.num_intervals() / 2)
        count = 0
        for i in range(self.num_intervals()): 
            if new_permutation[i] - 1 > i:
                m[count, self._index_of_interval[(i - k) % self.num_intervals()]] = 1
                count += 1

        return self.TransitionData(tr_matrix = m, new_perm = new_permutation)
    
    def reverse_data(self):
        reversing_permutation = Permutation(range(self.num_intervals(), 0, -1))
        new_permutation = reversing_permutation * self.permutation() * \
                reversing_permutation
        m = matrix(self.num_intervals() / 2)
        count = 0
        for i in range(self.num_intervals()): 
            if new_permutation[i] - 1 > i:
                m[count, self._index_of_interval[self.num_intervals() - 1 - i]] = 1
                count += 1

        return self.TransitionData(tr_matrix = m, new_perm = new_permutation)
    


    def get_interval(self, point):
        interval = self._intex.in_which_interval(point)- 1
        if abs(point - self._intex.domain_singularities()[interval]) < epsilon or\
                abs(point - self._intex.domain_singularities()[interval + 1]) < epsilon:
                    raise ValueError, 'Saddle connection found.'
        return interval
    
    def first_intersection(self, left_endpoint, right_endpoint, index_of_singularity):
        point = self._intex.domain_singularities()[index_of_singularity]
        interval_intersection_count = [0] * self.num_intervals()
        from_above = True
        while not is_between(left_endpoint, right_endpoint, point):
            point = _mod_one(point + 1/2)
           
            if is_between(left_endpoint, right_endpoint, point):
                from_above = False
                break
            interval = self.get_interval(point)
            point += self._translations[interval]
            interval_intersection_count[interval] += 1
            self.get_interval(point)

        coefficients = self._divpoints[index_of_singularity].coordinates
        coefficients += vector(interval_intersection_count) * self._translation_matrix
        calculated_point = coefficients * vector(self._lengths)
        coefficients -= floor((calculated_point - point) * 2 + 0.5) * vector([1] * (self.num_intervals() // 2))
        #print calculated_point, point
        calculated_point = coefficients * vector(self._lengths)
        #print calculated_point
        return (PointWithCoordinates(point, coefficients), from_above)

    def restrict_data(self, singularity):
        left_endpoint = self._divpoints[singularity]
        right_endpoint = self._divpoints[self.permutation().inverse()[singularity] - 1]

        intersections = [self.first_intersection(left_endpoint.point, right_endpoint.point, i)
                for i in range(self.num_intervals())]
        total = right_endpoint - left_endpoint

        def distance_from_right_endpoint(point_with_coordinates, is_from_above):
            diff = right_endpoint - point_with_coordinates
            if is_from_above:
                return diff 
            else:
                return PointWithCoordinatesNonOr(diff.point + total.point, 
                        diff.coordinates + total.coordinates)

        intersections = [(x, distance_from_right_endpoint(x[0], x[1]))
                for x in intersections]

        tuples = list(enumerate(intersections))
        tuples.sort(key = lambda x: -x[1][1].point)
     
        pairing = [0] * self.num_intervals()
        for i in range(self.num_intervals()):
            pairing[tuples[i][0]] = i 
        new_permutation = Permutation([pairing[self.permutation()[tuples[i][0]] - 1] + 1
                for i in range(self.num_intervals())])
        
        m = matrix(self.num_intervals() // 2)
        count = 0
        for i in range(self.num_intervals()):
            if new_permutation[i] - 1 > i:
                m[count] = (tuples[i][1][1] - tuples[i + 1][1][1]).coordinates
                count += 1

        return self.TransitionData(tr_matrix = m, new_perm = new_permutation)

    def new_foliation(self, transition_data):
        new_vector = transition_data.tr_matrix * vector(self._lengths)
        return FoliationNonOrientableSurface(transition_data.new_perm, new_vector.list())

    def find_pseudo_anosov_candidates(self, depth, transition_matrix_so_far = [], 
            permutation_so_far = [], 
            encoding_sequence_so_far = [],
            initial_permutation = []): 
        if transition_matrix_so_far == []:
            transition_matrix_so_far = matrix.identity(self.num_intervals() // 2)
            permutation_so_far = self.permutation()
            initial_permutation = self.permutation()
            _global_storage = []
            
        ret_list = []
                
        for i in range(self.num_intervals()):
            for j in range(2):
                if j == 0:
                    tr_data = self.rotation_data(i) 
                    final_encoding = encoding_sequence_so_far + [i, 'not reversed']
                    final_foliation = self.new_foliation(tr_data) 
                    final_matrix = tr_data.tr_matrix * transition_matrix_so_far
                    final_permutation = tr_data.new_perm 
                else:
                    tr_data = final_foliation.reverse_data()
                    final_encoding = encoding_sequence_so_far + [i, 'reversed']
                    final_matrix = tr_data.tr_matrix * final_matrix
                    final_permutation = tr_data.new_perm

                if final_permutation == initial_permutation:
                    ret_list.extend([(x, final_encoding)
                        for x in get_pseudo_anosov_candidates(matrix(RDF, final_matrix), False)])
                else:
                    global _global_storage
                    _global_storage.append('not matching permutation')
            
        if depth > 0:
            for i in range(self.num_intervals()):
                try:
                    tr_data = self.restrict_data(i)
                    ret_list += [pa for pa in self.new_foliation(tr_data).\
                            find_pseudo_anosov_candidates(depth - 1,\
                        tr_data.tr_matrix * transition_matrix_so_far,\
                        tr_data.new_perm, encoding_sequence_so_far + [i], \
                        initial_permutation)] 
                except ValueError:
                    pass
            
        return ret_list

    def verify_candidate(self, pseudo_anosov_candidate, exact_check = False):
        try:
            new_fol = foliation = FoliationNonOrientableSurface(self.permutation(),\
                    pseudo_anosov_candidate[0][1])
        except ValueError:
            #print 'Candidate with immediate saddle connection:', pseudo_anosov_candidate
            return False

        try:
            transition_matrix = matrix.identity(self.num_intervals() // 2)
            for i in range(len(pseudo_anosov_candidate[1]) - 1):
                if i == len(pseudo_anosov_candidate[1]) - 2:
                    tr_data = new_fol.rotation_data(pseudo_anosov_candidate[1][i])
                else:
                    tr_data = new_fol.restrict_data(pseudo_anosov_candidate[1][i])
                transition_matrix = tr_data.tr_matrix * transition_matrix
                new_fol = new_fol.new_foliation(tr_data)
            if pseudo_anosov_candidate[1][-1] == 'reversed':
                tr_data = new_fol.reverse_data()    
                transition_matrix = tr_data.tr_matrix * transition_matrix
                new_fol = new_fol.new_foliation(tr_data)
        except ValueError:
            #print 'Candidate with non-immediate saddle connection:', new_fol, pseudo_anosov_candidate
            return False

        ret_list = []

        if pseudo_anosov_candidate[0][0] in QQbar:
            for cand in get_pseudo_anosov_candidates(matrix(QQ, transition_matrix), False):
                if pseudo_anosov_candidate[0] == cand:
                    #print transition_matrix.eigenvalues()
                    #print transition_matrix.charpoly()
                    #print pseudo_anosov_candidate
                    return True
        elif pseudo_anosov_candidate[0][0] in CDF:
            for cand in get_pseudo_anosov_candidates(matrix(RDF, transition_matrix), False):
                if pseudo_anosov_candidate[0] == cand:
                    if not exact_check:
                        print matrix(RDF, transition_matrix).eigenvalues()
                        print transition_matrix.charpoly()
                        return True
                    qqbar_candidates = get_pseudo_anosov_candidates(matrix(QQ, transition_matrix), False)
                    for x in qqbar_candidates:
                        if abs(x[0] - cand[0]) < epsilon and sum([abs(x[1][k] - cand[1][k]) for k in range(len(x[1]))]) < epsilon:
                            #print pseudo_anosov_candidate
                            return self.verify_candidate((x, pseudo_anosov_candidate[1]))
                else:
                    return False
        else:
            raise ValueError, 'Unexpected number field'
    
    def find_pseudo_anosovs(self, depth):
        return [cand for cand in self.find_pseudo_anosov_candidates(depth) if self.verify_candidate(cand)]

def _is_involution_good(involution):
    half_size = involution.size() // 2
    i = 0
    while i < half_size:
        for j in range(half_size):
            if i < involution[i + j] - 1 < i + half_size: 
                i = i + j + 1
                break
        else:
            return False
    return True

#def matrix_with_charpoly(polynomial):
#    m = matrix(polynomial.degree())
#    coeffs = polynomial.coeffs()
#    for i in range(m.nrows() - 1):
#        m[i, i + 1] = 1
#    for i in range(m.nrows()):
#        m[-1, i] =  -coeffs[i] / coeffs[-1]
#    return m

def arnoux_yoccoz_factor(genus, field = CDF):
    R.<x> = ZZ[]
    poly = R([-1] * genus + [1])
    return max([abs(x[0]) for x in poly.roots(CDF)])



def arnoux_yoccoz_permutation(genus):
    l = []
    for i in range(genus):
        l.append(2 * i + 2)
        l.append(2 * i + 1)

    return Permutation(l)

def arnoux_yoccoz_foliation(genus, field = CDF):
    p = arnoux_yoccoz_permutation(genus)
    sf = arnoux_yoccoz_factor(genus, field)
    lengths = [1/sf**i for i in range(genus)]
    return FoliationNonOrientableSurface(p, lengths)

f = arnoux_yoccoz_foliation(3).to_foliation()
g = Foliation(Permutation([1]), [1], 0.2)
h = arnoux_yoccoz_foliation(3)
hh = FoliationNonOrientableSurface(arnoux_yoccoz_permutation(3), [0.261016378495, 0.1766049821, 0.0623786394048])
hhh = FoliationNonOrientableSurface(arnoux_yoccoz_permutation(3), [0.262714, 0.155554, 0.0817322]) 
h8 = FoliationNonOrientableSurface(Permutation([3,5,1,6,2,4,8,7]),[0.0944422795324, 0.178500974036, 0.074245928173, 0.152810818258])
h8_2 = FoliationNonOrientableSurface(Permutation([3,5,1,6,2,4,8,7]), [0.0820003383086, 0.194387134927, 0.065494592731, 0.158117934033])    
