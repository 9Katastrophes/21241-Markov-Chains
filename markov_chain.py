import math
import numpy as np
import scipy as spy

# checks that a matrix P is a markov transition matrix
# takes in a matrix represented as a list of lists
# 1. checks that matrix square
# 2. checks that elements of P are nonnegative and not greater than 1
# 3. checks that the sum the elements in each column
#    is equal to 0
# 4. checks that P has 1 as an eigenvalue
def is_markov_transition_matrix(P):
    # if P is not a list, raise Exception
    if (not isinstance(P, list)):
        raise Exception("Error: P is not a matrix")

    # if an element of P is not a list, raise Exception
    for i in range(len(P)):
        if (not isinstance(P[i], list)):
            raise Exception("Error: P is not a matrix")

    dim = len(P) # dimension of P

    # if P is not square, raise Exception
    for i in range(dim):
        if (len(P[i]) != dim):
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
        if (sum == 0):
            for j in range(dim):
                P[i][j] = 1/dim
        if (sum != 1):
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

# takes in a probability matrix A and a pagerank vector x
# outputs the ranks of the elements of the original list
# Ex: might take the form ranks=[2, 1, 3, 0], indicating that the most
# important node is 2, and the least important is 0
def calculate_rank(A, x):
    A = np.array(A).T
    x = np.array(x)
    evects = np.linalg.eig(A)[1]

    print(A)

    D = [[0 for i in range(len(A))] for j in range(len(A))]
    D[0][0] = 1
    D = np.array(D)

    #x = np.matmul(np.matmul(np.matmul(evects, D), np.linalg.inv(evects)), x)
    x = evects[0]

    print(x)

    rank = []
    for i in range(len(x)):
        max_index = 0
        max = 0
        for j in range(len(x)):
            if (abs(x[j]) > max):
                if (j in rank):
                    continue
                max_index = j
                max = abs(x[j])
        rank.append(max_index)

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
    # POSSIBLE ERRORS WITH FUNCTIONS AHEAD

    # analyzing rank using diagonalization
    print("Analyzing rank...")
    result = calculate_rank(A, x)

    return result