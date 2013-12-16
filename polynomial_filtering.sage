from itertools import product

epsilon = 1e-8

def passes_newton_test(polynomial, largest_root_bound):
    derivative = polynomial.derivative()
    x = largest_root_bound
    prev_x = 100000000
    count = 0
    while abs(x - prev_x) >= epsilon:
        count += 1
        if count > 100:
            print x, y, m
        prev_x = x
        y = polynomial(x)
        if y <= 0:
            return False
        m = derivative(x)
        if m <= epsilon:
            return False
        #print x, y, m
        x = x - y/m
        if x <= 1.00001:
            return False
    return True
        
def _add_next_possible_coeffs(ret_list, first_coeffs, first_traces, 
        remaining_trace_ranges):
    #print first_coeffs, first_traces, remaining_trace_ranges
    if len(remaining_trace_ranges) == 0:
        ret_list.append(first_coeffs)
        return
    next_coeff = 0
    n = len(first_coeffs)
    for j in range(n):
        next_coeff += (-1)^j * first_traces[j] * first_coeffs[n - 1 - j]
    try:
        next_trace = next(x for x in remaining_trace_ranges[0] 
                if (next_coeff + (-1)^n * x) % (n + 1) == 0)
    except StopIteration:
        return
    next_coeff += (-1)^n * next_trace
    next_coeff /= n + 1
    while next_trace in remaining_trace_ranges[0]:
        _add_next_possible_coeffs(ret_list, first_coeffs + [next_coeff],
                first_traces + [next_trace], remaining_trace_ranges[1:])
        next_trace += n + 1
        next_coeff += (-1)^n

def get_possible_coefficients(trace_ranges):
    #print trace_ranges
    ret_list = []
    _add_next_possible_coeffs(ret_list, [], [], trace_ranges)
    return ret_list

def find_polynomial_candidates(degree, largest_root_bound, is_orientable):
    if is_orientable:
        assert(degree % 2 == 0)
        trace_bounds = [floor(degree / 2 * (largest_root_bound^i +\
            1/largest_root_bound^i)) for i in range(1, degree / 2 + 1)]
        choice_for_constant = 1
    else:
        trace_bounds = [floor((degree // 2) *(largest_root_bound**i + \
            1/largest_root_bound**i)) + degree % 2 
            for i in range(1, degree)]
        choice_for_constant = 2

    trace_ranges = [range(trace_bounds[0] + 1)] + \
        [range(-trace_bounds[i], trace_bounds[i] + 1)
                for i in range(1, len(trace_bounds))]
    print "Number of possibilities for traces: ", \
            round(exp(sum([ln(len(s)) for s in trace_ranges])))
    possible_coeffs = get_possible_coefficients(trace_ranges)
    print "Of these, polynomials with integer coefficients: ", \
            len(possible_coeffs)

    ret_list = []
    R.<x> = ZZ[] 
    for middle_coeffs in possible_coeffs:
        if is_orientable:
            coeffs = [1] + middle_coeffs + middle_coeffs[-2::-1] + [1]
        else:
            coeffs = [1] + middle_coeffs + [1]
            if (coeffs[1] + coeffs[-2]) % 2 != 0:
                continue
        for i in range(2):
            if i == 1:
                if coeffs[1] > 0:
                    coeffs = [coeffs[k] * (-1)^k for k in range(len(coeffs))]
                else:
                    continue
            for j in range(choice_for_constant):
                coeffs[-1] = (-1)^j
                poly = R(coeffs[::-1])
                if not passes_newton_test(poly, largest_root_bound):
                    continue
                if set(poly.coeffs()[1::2]) == {0}:
                    continue
                roots = poly.roots(CDF)
                abs_roots = sorted([abs(x[0]) for x in roots])
                for root_tuple in roots:
                    root = root_tuple[0]
                    if root_tuple[1] == 1 and root.imag() == 0 and \
                            1 + epsilon < root < largest_root_bound and \
                            root > abs_roots[-2] + epsilon and \
                            1/root - epsilon < abs_roots[0]:
                                ret_list.append(poly)
                                ret_list.append(roots)
    return ret_list


#def find_min_dilatation_candidates(genus, upper_bound, is_orientable):



