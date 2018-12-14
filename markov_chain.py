import math
import numpy as np
import scipy as spy

# helper function that computes the result of multiplying a 
# scalar x with a matrix A
def scalarMult(A, x):
    if (not isinstance(A, list) or not(isinstance(x, int) or isinstance(x, float))):
        raise Exception("Error: invalid input for scalar mult")
    for i in range(len(A)):
        if (not isinstance(A[i], list)):
            raise Exception("Error: invalid input for scalar mult")

    Ax = []
    for i in range(len(A)):
        Ax.append(x*A[i])

    return Ax

# helper function that computes the dot product of two vectors
def dotProduct(a, b):
    if (not isinstance(a, list) or not isinstance(b, list)):
        raise Exception("Error: invalid input for dot product")
    if (len(a) != len(b)):
        raise Exception("Error: invalid length for vectors in dot product")

    sum = 0
    for i in range(len(a)):
        sum += a[i]*b[i]

    return sum
        
# helper function that computes the result of matrix multiplying A and v
def matrixMultiply(A, v):
    if not isinstance(A, list) or not isinstance(v, list):
        raise Exception("Error: invalid input for matrix multiply")
    if (len(A[0]) != len(v)):
        raise Exception("Error: invalid dimensions matrix multiplication")

    Av = []
    for i in range(len(A)):
        Av.append(dotProduct(A[i], v))
        
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
        if not isinstance(P[i]):
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
    if (1 in eig) == False:
        raise Exception("Error: 1 is not an eigenvalue of P")

# takes in a pagerank vector x that has been multiplied by A many times
# outputs the ranks of the elements of the original list
# Ex: might take the form ranks=[2, 1, 3, 0], indicating that the most
# important node is 2, and the least important is 0
def calculate_rank(x):
    rank = []
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

    # create a modified matrix to account for reducible web graphs
    # (reducible web graphs do not have bidirectional paths from each node t0
    # each other node)
    # 1. make an all 1s square matrix, all_ones, of size dim
    # 2. then perform M = dA + (1-d)/dim(all_ones)
    print("Creating modified matrix M...")
    all_ones = [[1 for i in range(num_nodes)] for i in range(num_nodes)]
    M = d*A + ((1-d)/num_nodes)*all_ones
        
    # perform x = Mx for an arbitrarily large number of times
    # until M converges
    print("Repeatedly muliplying M and x...")
    k = 10
    for i in range(k):
        x = matrixMultiply(M, x)

    # the resulting x after k iterations gives the rank of each of the original elements
    print("Analyzing rank...")
    result = calculate_rank(x)

    return result