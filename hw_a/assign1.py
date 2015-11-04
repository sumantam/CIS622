import sys
import numpy
import os
import time

def diagonalize(vec):

    len = numpy.shape(vec);
    no_of_ele = len[0];

    X = numpy.zeros(shape=(no_of_ele, no_of_ele));

    # Put down the elements along the diagonals

    i = 0;

    while not (i >= no_of_ele) :
        X[i, i] = vec[i];
        i += 1;

    return X

# The following are the global variables

n = 0;
epsilon = 0.0;
beta = 0.23;

# Start of the main program

#print ("Hello World");
n = int(input("What is the value of number of variables ?"));

while True:
    epsilon = float(input("What is the value of the epsilon ?"));

    if (epsilon >= 0.5) :
        print("The value of the epsilon should be a positive number less than 0.5");
        continue
    else:
        break

while True:
    beta = float(input("What is the value of the beta ?"));

    if (epsilon >= 1) :
        print("The value of the epsilon should be a positive number less than 1");
        continue
    else:
        break

print("The value of n is", n);

#
# This is the matrix which is fully zeros.
#

A = numpy.zeros(shape=(2*n, 3*n));

#print ("This is the first print");
#print(A);


row = 0;
while  not (row >= n) :
    col = row;
    #while not (col > row) :
    A[2*row + 1, 3*col] = 1;
    A[2*row + 1, 3*col + 2] = 1;

    A[2*row , 3*col] = 1;
    A[2*row , 3*col + 1 ] = -1

    if (row  > 0) :
        A[2*row + 1, 3*col -3] = epsilon;
        A[2*row , 3*col -3] = -epsilon;

    #    col += 1;
    row += 1;

#print("This is the second print");
#print (A)

# populate the value of b

b = numpy.zeros(shape=(2*n, 1));

b[0] = epsilon;

counter = 1;
while not (counter > (2*n-1)) :
    b[counter] = counter%2;
    counter +=1;

#print("The value of b is being printed");
#print(b);

#
# The affine scaling algorithm
#
# Since the objective function is max x_n so it is ok to do
#
# min (-x_n)
# subject to the constraints
#

# The initialization of the affine scaling algorithm is done in the following way

e = numpy.zeros(shape=(3*n, 1));

e.fill(1);

#print(e);

ex = b - numpy.dot(A, e);

print(ex);

# Extend the matrix to a larger matrix and insert the elements

C = numpy.concatenate((A, ex), axis=1);

#
# Hence now the problem is the following
# min (-x_n) + 2000000*x_(n+1)
#
# subject to the following conditions.
# Cx = b
#

#
# Initilaize X
#

elements_C = numpy.shape(C)
size = elements_C[1];

#
# In this case the objective function is being hardcoded.
# One should provide a way of having the objective function not hardcoded.
#

c = numpy.zeros(shape=(size, 1));

c[size-2, 0] = -1;
c[size-1, 0] = 20000;

print(c);



vec = numpy.ones(shape=(size,1));
ele = numpy.ones(shape=(size,1));

#print (X);

while True:
    X = diagonalize(vec);

    print(X)
    print("***********");

    p1 = numpy.linalg.inv(numpy.dot(numpy.dot(C, numpy.dot(X,X)), numpy.transpose(C)));
    p2 = numpy.dot(numpy.dot(C, numpy.dot(X, X)), c);

    p = numpy.dot(p1, p2);
    r = c - numpy.dot(numpy.transpose(C), p);

    error = numpy.dot(numpy.dot(numpy.transpose(ele), X), r);

    print ("The value of the vector r is ", r)

    print ("The error is ", error);

    if ((r>=0).all() and (not (error >= beta))):
        break;
    else:
        # Unboundedness check
        if ((-numpy.dot(numpy.dot(X, X), r) > 0).all()):
            break;
        else:
            vec = vec - beta*numpy.dot(numpy.dot(X, X), r)/(numpy.linalg.norm(numpy.dot(X,r)));
            print("The value of the vector is ", vec);
            time.sleep(2);

    if not (vec >= 0).all():
        print ("There is a HUGE PROBLEM WITH YOUR CODE");
        break;

    print("-------------##############################################----------------------");