good_discriminants = {148, 621, 1492, 2917, 564, 8173, 12436, 18117, 25492, 34861}

def search_surfaces_with_given_trace_field(good_discriminants,
        repeat = 10):
    count = 0
    ret_set = set()
    while count < repeat:
        m = matrix.random(ZZ, 3, 3, x=0, y=3)
        S = m*m.transpose()
        try:
            N.<a> = NumberField(S.charpoly())
        except:
            continue
        d = N.disc()
        mu = max(x[0] for x in (S.charpoly()).roots(ring =CDF))
        m.set_immutable()
        if d in good_discriminants:
            ret_set.add((m, d, S.charpoly(), mu))
            count += 1
    return ret_set

def get_products(A, B, how_much_to_go, current_list = []):
    if how_much_to_go == 0:
        return current_list
    if current_list == []:
        new_list = [A, A^(-1), B, B^(-1)]
    else:
        new_list = flatten([[A*M, A^(-1)*M, B*M, B^(-1)*M] for M in current_list[-1]])
    return get_products(A, B, how_much_to_go - 1, current_list + [new_list])


def get_traces(surface, depth = 5):
    traces = set()
    mu = surface[-1]
    A = matrix([[1,mu^(1/2)],[0,1]])
    B = matrix([[1,0],[mu^(1/2),1]])
    for m in flatten(get_products(A, B, depth)):
        t = m.trace()
        if 2.0001 < abs(t):
            traces.add(t.n(prec = 20))
    return sorted(list(traces), key = lambda x: abs(x))
    
def search_pseudo_anosov_traces(good_discriminants, repeat):
    surfaces = search_surfaces_with_given_trace_field(\
            good_discriminants, repeat)
    for s in surfaces:
        print s, get_traces(s, depth = 6)


search_pseudo_anosov_traces([148], repeat = 10)

