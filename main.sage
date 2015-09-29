load('dproj.sage')

import sys, getopt

# Default root system
root_system = 'A2'
# Use the standard representation by default
n_of_v_0 = 1


optlist, args = getopt.getopt(sys.argv[1:], 'R:n:')
if len(optlist) == 0:
    print 'Usage: sage dproj.sage -R root_system -n n_of_v0'
    print 'By default, use root_system = "A2" and n_of_v0 = 1.'

for opt, arg in optlist:
    if(opt == '-R'):
        root_system = arg
    elif(opt == '-n'):
        n_of_v_0 = int(arg)

c = get_dproj_matrix(root_system, n_of_v_0)

print "Donagi's projection matrix:"
print_matrix(c)
print "\nEigenvectors and Eigenvalues:"
print matrix(c).eigenvectors_right()
