from itertools import product

epsilon = 1e-8

def find_polynomial_candidates(degree, largest_root_bound):
    trace_bounds = [floor((degree - 1) * largest_root_bound**i + \
            1/largest_root_bound**i) for i in range(1, degree)]
    R.<x> = ZZ[] 
    range_list = [range(trace_bounds[0])] + [range(-trace_bounds[i], trace_bounds[i] + 1) \
            for i in range(1, degree - 1)]
    for traces in product(*range_list):
        coeffs = [1]
        for i in range(len(traces)):
            next_coeff = 0 
            for j in range(i + 1):
                next_coeff += (-1)^j * traces[j] * coeffs[i - j]
            if next_coeff % (i + 1) != 0:
                break
            else:
                coeffs.append(next_coeff / (i + 1))
        else:
            if (coeffs[1] + coeffs[-1]) % 2 == 0:
                coeffs.append(0)
                for i in range(2):
                    coeffs[-1] = (-1)^i
                    poly = R(coeffs[::-1])
                    roots = poly.roots(CDF)
                    abs_roots = sorted([abs(x[0]) for x in roots])
                    for root_tuple in poly.roots(CDF):
                        root = root_tuple[0]
                        if root_tuple[1] == 1 and abs(root.imag()) < epsilon and \
                                1 + epsilon < abs(root) < largest_root_bound and \
                                abs(root) == abs_roots[-1] and \
                                1/abs(root) + epsilon < abs_roots[0]:
                            print traces, poly, poly.roots(CDF)





