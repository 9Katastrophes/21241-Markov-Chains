import markov_chain

#----------------------#
# Markov Chain Samples #
#----------------------#

# INPUT: 
print("***Begin test 1***\n")

input = [[1, 5], [2, 5], [1, 3, 5], [4], [1, 5], [2, 6], [0, 1]]


# OUTPUT:
output = markov_chain.rank(input)

print("Expected output")
print([5, 2, 1, 6, 3, 4, 0])

print("Actual output")
print(output)

  
#INPUT:
print("***Begin test 2***\n")

input = [[1,3,4],[0,2,4],[3,6],[2,4,6],[5,8],[4,6,8],[0,7,9],[0,6,8],[2,9],[0,2,8]]


#OUTPUT:
output = markov_chain.rank(input)

print("Expected output")
print([2,6,8,0,3,9,4,5,7,1])

print("Alternative output")
print([2,6,8,0,4,3,9,5,7,1])

print("Actual output")
print(output)

#(or, with a slightly different strategy, perhaps [2,6,8,0,4,3,9,5,7,1]