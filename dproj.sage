separator = '' 
for i in range(80):
    separator += '#'

def print_matrix(m):
    for row in m:
        print row

# Change the representation of a weight vector 
# from the linear combiation of fundamental weights
# to a vector in the ambient orthonormal basis. 
def weight_to_ambient(w, ambient_space):
    vec = 0
    omega = ambient_space.fundamental_weights()
    for i, c_i in enumerate(w.to_vector(), start=1):
        vec += c_i * omega[i]
    return vec

# return the quotient group W/W_P. 
# W is the Weyl group.
# W_P is the stabilizer of weight v_0.
def get_W_W_P(v_0, W):
    #print 'W = {}'.format(W)

    W_P_gens = [W.unit()]
    for i in v_0.weyl_stabilizer():
        W_P_gens.append(W.simple_reflection(i))
    W_P = W.subgroup(W_P_gens)
    #print 'W_P = {}'.format(W_P)

    W_W_P = []
    for w_i in W:
        coset = []
        for w_j in W_P:
            coset.append(w_i * WeylGroupElement(W, w_j))

        sorted_coset = sorted(coset, reverse=True)
        if sorted_coset in W_W_P:
            continue
        else:
            W_W_P.append(sorted_coset)

    print 'W/W_P = ['
    for coset in W_W_P:
        coset_in_string = '['
        for e in coset:
            coset_in_string += '{}.W_P, '.format(e.reduced_word())
        coset_in_string += ']'
        print coset_in_string
    print ']'

    return W_W_P


# Find the result of the Donagi projection.
# v_0 is the weight fixed by W_P.
# w_i is a representative of the coset W/W_P
# W is the Weyl group.
def dproj(v_0, w_i, W, W_W_P):
    print 'Find the Donagi projection of w_i.W_P = {}.W_P '.format(
        w_i.reduced_word()
    )
    print 'with respect to v_0 = {}.'.format(v_0)

    c_i = []
    proj = 0
    proj_str = ''

    for j, coset in enumerate(W_W_P):
        w_j = coset[0]
        v_j = v_0.weyl_action(w_j)
        c_i_j = v_0.weyl_action(w_i).symmetric_form(v_j)
        c_i.append(c_i_j)
        proj_str += '({}){}.W_P'.format(
            c_i_j, 
            w_j.reduced_word(),
        )
        if j < len(W_W_P) - 1:
            proj_str += ' + '

    print '{}W_P |-> {}'.format(
        w_i.reduced_word(), proj_str
    ) 
    
    return c_i

def get_dproj_matrix(root_system, n_of_v_0):
    print separator
    R = RootSystem(root_system)
    print 'Consider {}'.format(R)
    L = R.weight_space()
    A = R.ambient_space()
    W = WeylGroup(R)

    print separator
    print 'Representation of Weyl group elements '
    print 'as a word `[i_1,i_2,\ldots,i_k]` of minimal length'
    print 'such that `s_{i_1} s_{i_2} \cdots s_{i_k}=self`,'
    print 'where `s` are the simple reflections.'
    for i, w in enumerate(W):
        print 'W[{}] = {}'.format(i, w.reduced_word())

    v_0 = L.fundamental_weight(n_of_v_0)
    print 'v_0 = {} = {},'.format(v_0, weight_to_ambient(v_0, A))
    print 'where Lambda[i] is the i-th fundamental weight, '
    print 'and the vector representation is in the orthonomal basis.'

    print separator

    print 'Find W/W_P.'
    W_W_P = get_W_W_P(v_0, W)

    print separator

    c = []

    for i, coset in enumerate(W_W_P):
        w_i = coset[0]
        print 'w_{} = {}'.format(i, w_i.reduced_word())
        v_i = v_0.weyl_action(w_i)
        print '(w_{}).(v_0) = {} = {}'.format(
            i, v_i, weight_to_ambient(v_i, A)
        )

        # double-check
        for w_i_2 in coset[1:]:
            v_i_2 = v_0.weyl_action(w_i_2)
            if v_i != v_i_2:
                print 'Error! the action of W/W_P on v_0 is ill-defined.'

        c_i =  dproj(v_0, w_i, W, W_W_P) 
        c.append(c_i)

        print separator

    return matrix(c)

