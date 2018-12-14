# Markov Chains

We will be building an algorithm to rank states in a Markov chain based
on importance. This is similar to a Google PageRank algorithm, using random walks to
estimate the importance of a website.

The structure of this kind of method looks something as follows:
• We begin with a Markov chain on a finite state space. You can think of this is a list
of websites, with the links that go between them begin indicated.
• From here, we generate a matrix P that indicates the probabilities of moving from one
site to another. That is, the ij entry of the matrix is P(Xn = i | Xn−1 = j), the chance
of moving from the j
th element to the i
th. You should be able to verify that if you
multiply P by a vector indicating your current position (represented as the probability
you are in a given state) you should get a vector representing the probabilities of
landing in each of the states after one random step.
• Using the Perron-Frobenius Theorem, we can observe that the largest eigenvalue of
this matrix, even up to absolute value, is 1, and it has a positive eigenvector. All the
other eigenvalues will have absolute value less than 1.
• Now, suppose we begin our random walk at some vector x. If we write x as a linear
combination of eigenvectors and consider the limit limn→∞ P
nx, notice what happens
to the eigenvalues... only 1 does not go to 0. (Specifics, of course, can be found in the
reading!)
• Hence, the limit of the random walk is always the same. We use this as our ranking
vector, and rank nodes higher if we are more likely to visit them on our walk, in the
limit.
