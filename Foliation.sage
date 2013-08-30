from sage.structure.sage_object import SageObject
from bisect import bisect

def _mod(a, b):
    x = a % b
    if x < 0:
        return x + b
    else:
        return x

#def _post_shift(permutation, k):
#    n = permutation.size()
#    return Permutation([_mod(permutation[i] - 1 + k, n) + 1 for i in range(n)])

def rotating_permutation(n, k):
    return Permutation([(i + k) % n + 1 for i in range(n)])

def rotating_matrix(n, k):
    return rotating_permutation(n, k).to_matrix()

def singularity_permutation(permutation):
    p = permutation
    p_inv = p.inverse()
    n = p.size()
    return Permutation([p_inv[p[(i - 1) % n] % n] for i in range(n)])



def _canonical_form(intex, twist):
    one_minus_twist = _mod(1 - twist, 1)
    rotateby = bisect(intex.range_singularities(), one_minus_twist)
    for i in [rotateby - 1, rotateby]:
        if abs(one_minus_twist - intex.range_singularities()[i]) < 0.000000001:
            raise ValueError, 'The foliation has a closed leaf. It should not.'
    newtwist = intex.range_singularities()[rotateby] - one_minus_twist 
    n = len(intex.lengths())
    return (iet.IntervalExchangeTransformation(rotating_permutation(n, rotateby) *\
            intex.permutation().to_permutation(), intex.lengths()), newtwist)

class Foliation(SageObject):
    r"""
    Oriented projective measured foliation on an orientable surface
    """
    def __init__(self, permutation, lengths, twist):
        intex = iet.IntervalExchangeTransformation(permutation, lengths)
        normalized_twist = twist / intex.length()
        intex = intex.normalize()
        ret = _canonical_form(intex, twist)
        self._intex = ret[0]
        self._twist = ret[1]
        self._singularity_cycles = singularity_permutation(self._intex.permutation().to_permutation()).to_cycles()

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
        normalized_k = _mod(k, self.num_intervals())
        rotating_dist = sum(self._intex.lengths()[self.num_intervals() - normalized_k:])
        containing_interval = bisect(self._intex.range_singularities(), _mod(-rotating_dist- self._twist, 1)) 
        new_permutation = rotating_permutation(self.num_intervals(), containing_interval) *\
                self.permutation() * \
                rotating_permutation(self.num_intervals(), k)
        m = matrix(self.num_intervals() + 1)
        for i in range(self.num_intervals()):
            m[i, _mod(i - k, self.num_intervals())] = 1
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

    def rotate(self, k):
        x = self.rotation_data(k)
        return self.new_foliation(x[0], x[1])


class FoliationNonOrientableSurface(SageObject):
    r"""
    Normalized measured foliation with \mathbb{Z}_2 holonomy on a non-orientable surface.
    """

    def __init__(self, permutation, lengths):
        if permutation.number_of_fixed_points() > 0:
            raise ValueError, "The pormutation should not have fixed points."
        if permutation.inverse() != permutation:
            raise ValueError, 'The permutation must be an involution.'
        if permutation.size() != 2 * len(lengths):
            raise ValueError, 'The list of lengths should be half as long as the size of the permutation'
        total = sum(lengths)
        self._lengths = [x/total/2 for x in lengths]

        detailed_lengths = 2 * lengths
        count = 0
        for x in range(2 * len(lengths)):
            if permutation[x] - 1 > x:
                detailed_lengths[x] = detailed_lengths[permutation[x] - 1] = self._lengths[count]
                count += 1
        self._intex = iet.IntervalExchangeTransformation(permutation, detailed_lengths)

    def __repr__(self):
        return 'Foliation with Z_2 holomony on a non-orientable surface of genus {0}.\n'\
                'Permutation: {1}\nLengths: {2}'.format(self.genus(), 
                        self._intex.permutation().to_permutation(), self._intex.lengths())

    def num_intervals(self):
        return self._intex.lengths().size()
    
    def to_foliation(self):
        return Foliation(self._intex.permutation().to_permutation(), self._intex.lengths(), 0.5)

    def euler_char(self):
        return self.to_foliation().euler_char() / 2

    def genus(self):
        return 2 - self.euler_char()

    def rotation_data(self, k):
        p = Permutation([(i + k) % self.num_intervals() + 1 for i in range(self.num_intervals())])
        m = p.to_matrix()
        return 

def _ay_entry(i,j,n):
    if (j - i) % n == n - 1 or j == n - 1:
        return 1
    else:
        return 0

def arnoux_yoccoz_factor(genus):
    m = matrix(genus, lambda i, j: _ay_entry(i, j, genus))
    return max([abs(x) for x in m.eigenvalues()])

def arnoux_yoccoz_foliation(genus):
    l = []
    for i in range(genus):
        l.append(2 * i + 2)
        l.append(2 * i + 1)

    p = Permutation(l)
    sf = arnoux_yoccoz_factor(genus)
    lengths = [1/sf**i for i in range(genus)]
    return FoliationNonOrientableSurface(p, lengths)

f = arnoux_yoccoz_foliation(3).to_foliation()
#f = Foliation(Permutation([2, 1, 4, 3, 6, 5]), [0.2718, 0.2718, 0.1478, 0.1478, 0.0803, 0.0803], 0.5)
g = Foliation(Permutation([1]), [1], 0.2)


