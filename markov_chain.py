import math
import numpy as np
import scipy as spy

# helper function that adds two matrices
def matrix_add(A, B):
    if (not isinstance(A, list) or not isinstance(B, list)):
        raise Exception("Error: cannot add matrices A and B because they have different dimensions")
    for i in range(len(A)):
        if (not isinstance(A[i], list) or not isinstance(B[i], list)):
            raise Exception("Error: cannot add matrices A and B because they have different dimensions")
    
    result = A.copy()

    for i in range(len(result)):
        for j in range(len(result[0])):
            result[i][j] + B[i][j]
    
    return result

# helper function that computes the result of multiplying a 
# scalar x with a matrix A
def scalar_mult(A, x):
    if (not isinstance(A, list) or not(isinstance(x, int) or isinstance(x, float))):
        raise Exception("Error: invalid input for scalar mult")
    for i in range(len(A)):
        if (not isinstance(A[i], list)):
            raise Exception("Error: invalid input for scalar mult")

    Ax = A.copy()
    for i in range(len(A)):
        for j in range(len(A)):
            Ax[i][j] = x * Ax[i][j]

    if (len(Ax) != len(Ax[0])):
        raise Exception("Error: scalar mult resulted in wrong dimensions")

    return Ax

# helper function that computes the dot product of two vectors
def dot_product(a, b):
    if (not isinstance(a, list) or not isinstance(b, list)):
        raise Exception("Error: invalid input for dot product")
    if (len(a) != len(b)):
        raise Exception("Error: invalid length for vectors in dot product")

    sum = 0
    for i in range(len(a)):
        sum += a[i]*b[i]

    return sum
        
# helper function that computes the result of matrix multiplying A and v
def matrix_multiply(A, v):
    if not isinstance(A, list) or not isinstance(v, list):
        raise Exception("Error: invalid input for matrix multiply")
    if (len(A[0]) != len(v)):
        raise Exception("Error: invalid dimensions matrix multiplication")

    Av = []
    for i in range(len(A)):
        Av.append(dot_product(A[i], v))
        
    return Av

# checks that a matrix P is a markov transition matrix
# takes in a matrix represented as a list of lists
# 1. checks that matrix square
# 2. checks that elements of P are nonnegative and not greater than 1
# 3. checks that the sum the elements in each column
#    is equal to 0
# 4. checks that P has 1 as an eigenvalue
def is_markov_transition_matrix(P):
    # if P is not a list, raise Exception
    if not isinstance(P, list):
        raise Exception("Error: P is not a matrix")

    # if an element of P is not a list, raise Exception
    for i in range(len(P)):
        if not isinstance(P[i], list):
            raise Exception("Error: P is not a matrix")

    dim = len(P) # dimension of P

    # if P is not square, raise Exception
    for i in range(dim):
        if len(P[i]) != dim:
            raise Exception("Error: P is not square")

    # if P contains a negative value or a value greater than 1, raise Exception
    for i in range(dim):
        for j in range(dim):
            if (P[i][j] < 0 or 1 < P[i][j]):
                raise Exception("Error: P contains a negative value or a value greater than 1")

    # if elements of each column of P do not add up to 1, raise Exception
    for i in range(dim):
        sum = 0
        for j in range(dim):
            sum += P[i][j]
        # if the sum of the elements of a column is 0, make each element 1/N
        if sum == 0:
            for j in range(dim):
                P[i][j] = 1/dim
        if sum != 1:
            raise Exception("Error: the elements in a column of P do not add up to 1")

    # if 1 is not an eigenvalue of P, raise Exception
    eig = np.linalg.eig(np.array(P))[0]
    tolerance = 0.1
    eig1_found = False
    for i in range(len(eig)):
        if (abs(eig[i] - 1) < tolerance):
            eig1_found = True
    if (eig1_found == False):
        raise Exception("Error: 1 is not an eigenvalue of P")

    return True

# takes in a pagerank vector x that has been multiplied by A many times
# outputs the ranks of the elements of the original list
# Ex: might take the form ranks=[2, 1, 3, 0], indicating that the most
# important node is 2, and the least important is 0
def calculate_rank(x):
    rank = []
    # print(x)
    for i in range(len(x)):
        max_index = 0
        max = x[max_index]
        for j in range(len(x)):
            if x[j] > max:
                max_index = j
                max = x[j]
        rank.append(max_index)
        x[max_index] = 0

    return rank

# main function
# takes in a markov chain, "links", and d a damping factor between 0 and 1
# outputs a set of rankings
def rank(links, d=.85):
    # count number of nodes
    num_nodes = len(links)

    # initialize markov transition matrix with 0s
    A = [[0 for i in range(num_nodes)] for j in range(num_nodes)]

    for i in range(num_nodes):
        # count number of edges for each node
        num_edges = len(links[i])
        for j in range(len(links[i])):
            # add the probabilities for each "path" from i to its neighbors
            neighbor = links[i][j]
            A[i][neighbor] = 1/num_edges

    # checking that A is a valid markov transition matrix
    print("Verifying markov transition matrix...")
    if (not is_markov_transition_matrix(A)):
        raise Exception("Error: rank verification failed")

    # x is the pagerank vector with elements initialized to 1
    print("Creating pagerank vector...")
    x = [1 for i in range(num_nodes)]

    # WORKS UP TO THIS POINT
    # POSSIBLE ERRORS WITH CUSTOM MATRIX MATH FUNCTIONS AHEAD

    # create a modified matrix to account for reducible web graphs
    # (reducible web graphs do not have bidirectional paths from each node t0
    # each other node)
    # 1. make an all 1s square matrix, all_ones, of size num_nodes
    # 2. then perform M = dA + (1-d)/dim(all_ones)
    print("Creating modified matrix M...")
    all_ones = [[1 for i in range(num_nodes)] for j in range(num_nodes)]
    M = matrix_add(scalar_mult(A, d), scalar_mult(all_ones, (1-d)/num_nodes))         # SOMETHING COULD BE WRONG HERE
    # perform x = Mx for an arbitrarily large number of times
    # until M converges
    print("Repeatedly muliplying M and x...")
    k = 20
    for i in range(k):
        x = matrix_multiply(M, x)                                                     # SOMETHING COULD BE WRONG HERE

    # the resulting x after k iterations gives the rank of each of the original elements
    print("Analyzing rank...")
    result = calculate_rank(x)

    return result