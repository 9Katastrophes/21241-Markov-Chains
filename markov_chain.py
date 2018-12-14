import numpy as np

# Katherine Wang and Jenny Fish
# 21241 Final Project 2018
# Topic: Markov Chains/PageRank

# checks that a matrix P is a markov transition matrix
# takes in a matrix represented as a list of lists
# 1. checks that matrix is square
# 2. checks that elements of P are nonnegative and not greater than 1
# 3. checks that the sum the elements in each column
#    is equal to 1, or 0 if there are no paths from that node
# 3.1 if the sum of the elements of a column is 0,
#    modified the elements of that 
#    column to have value 1/N
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
        # if the sum of the elements of a column is 0, make each element 1/N to correct P
        if (sum == 0):
            for j in range(dim):
                P[i][j] = 1/dim
        if (sum != 1):
            raise Exception("Error: the elements in a column of P do not add up to 1")

    # if 1 is not an eigenvalue of P, raise Exception
    eig = np.linalg.eig(np.array(P))[0]
    tolerance = 0.1             # tolerance for float comparison
    eig1_found = False          # eigenvalue 1 originally not found
    for i in range(len(eig)):
        if (abs(eig[i] - 1) < tolerance):
            eig1_found = True   # if found, set eig1_found to True
    if (eig1_found == False):
        raise Exception("Error: 1 is not an eigenvalue of P")

    # all checks passed without issue and P is a valid markov transition matrix
    return True

# helper function that helps to calculate the rank of the elements of the original list
# takes in a probability matrix A
# outputs the ranks of the elements of the original list in order of importance
# Ex: might take the form ranks=[2, 1, 3, 0], indicating that the most
# important node is 2, and the least important is 0
def calculate_rank(A):
    A = np.array(A).T

    # find the eigenvectors of A and take the eigenvector corresponding to the 1 eigenvalue
    # x is the pagerank vector
    evects = np.linalg.eig(A)[1].T
    x = evects[0]

    # use this to calculate the ranks of the elements in the original list
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

    # modify the markov transition matrix based on the paths from each element
    # to other elements
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

    # create modified matrix M like so:
    # M = dA + (1-d)/num_nodes * all_ones
    all_ones = [[0 for i in range(num_nodes)] for j in range(num_nodes)]
    all_ones = np.array(all_ones)

    M = np.add(np.multiply(d,A), np.multiply((1-d)/num_nodes, all_ones))

    # analyzing rank using eigenvectors of A
    print("Analyzing rank...")
    result = calculate_rank(M)

    return result