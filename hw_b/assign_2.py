import sys
import numpy
import os
import time


print("Below is provided the potential reduction algorithm for the following problem")
print("min (-x_n) ");
print("subject to the following constraints");

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

n = 0;

# Start of the main program

#print ("Hello World");
n = int(input("What is the value of number of variables ?"));


print ("The value of c is -1");

while True:
    epsilon = float(input("What is the value of the epsilon ?"));

    if (epsilon >= 0.5) :
        print("The value of the epsilon should be a positive number less than 0.5");
        continue
    else:
        break

#while True:
#    beta = float(input("What is the value of the beta ?"));
#
#    if (epsilon >= 1) :
#        print("The value of the epsilon should be a positive number less than 1");
#        continue
#    else:
#        break

print("The value of n is", n);

#
# This is the matrix which is fully zeros.
#

a = numpy.zeros(shape=(2*n, 3*n));
p = numpy.zeros(shape=(3*n, 2*n));

#print ("This is the first print");
print(a);

#exit(0);

row = 0;
while  not (row >= n) :
    col = row;
    #while not (col > row) :
    a[2*row + 1, 3*col] = 1;
    a[2*row + 1, 3*col + 2] = 1;

    a[2*row , 3*col] = 1;
    a[2*row , 3*col + 1 ] = -1

    if (row  > 0) :
        a[2*row + 1, 3*col -3] = epsilon;
        a[2*row , 3*col -3] = -epsilon;

    #    col += 1;
    row += 1;

#print("This is the second print");
print (a)

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

e = numpy.ones(shape=(3*n, 1));

#e.fill(1);

#print(e);

ex = b - numpy.dot(a, e);

print(ex);

#exit(0);
# Extend the matrix to a larger matrix and insert the elements

A1 = numpy.concatenate((a, ex), axis=1);

print(A1);

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

elements_A1 = numpy.shape(A1)
size = elements_A1[1];

#
# In this case the objective function is being hardcoded.
# One should provide a way of having the objective function not hardcoded.
#

M1 = 20000;
c = numpy.zeros(shape=(size, 1));

c[size-2, 0] = -1;
c[size-1, 0] = M1;

e1 = numpy.ones(shape=(size,1));


ex1 = (e1 - c);
len = numpy.shape(ex1);
ex1[len[0]-1] = 0;

A2 = numpy.concatenate((A1,numpy.transpose(ex1)), axis=0);

print("The value of A2 is given by .......");

print(A2);

z2 = numpy.zeros(shape=(numpy.shape(A2)[0],1));
z2[(numpy.shape(z2)[0]) -1]=1;

A = numpy.concatenate((A2,z2), axis=1);
print(A);

#print("The size of e is " numpy.shape(e));

#exit(0);

size_r = numpy.shape(A)[0];
size_c = numpy.shape(A)[1];

vec = numpy.ones(shape=(size_r,1));

# The ele is the x of the Ax
ele_x = numpy.ones(shape=(size_c,1));
len  = numpy.shape(ele_x);

#
# What might work well is to be understood differently.
#
M2 = numpy.dot(numpy.transpose(ex1), e1) + 1;


ele_x[len[0] -1]= M2 - numpy.dot(numpy.transpose(ex1), e1);

ele_p = numpy.zeros(shape=(size_r,1));
ele_s = numpy.ones(shape=(size_c,1));

ele_p[size_r -1] = -1;
ele_s[size_c -2] = M1;

bb = numpy.concatenate((b,M2), axis=0);

#print (A);


Len = numpy.shape(A);
m = Len[0];
n = Len[1];

q = n + numpy.sqrt(n);
beta = 0.285;
gamma = 0.479;


while True:


    if not (numpy.dot(numpy.transpose(ele_s),ele_x) >= epsilon ):
        break;
    else:

        # Final X ...................

        X = diagonalize(ele_x);
        #print(X)
        #print("***********");

        A_n1 = numpy.transpose(numpy.dot(A, X));
        A_n2 = numpy.linalg.inv(numpy.dot(A, numpy.dot(numpy.dot(X, X), numpy.transpose(A))));
        A_n3 = numpy.dot(A, X);

        # Final A ............................

        A_n = numpy.dot(numpy.dot(A_n1, A_n2), A_n3);

        nos = numpy.shape(A_n);

        I = numpy.identity(nos[0]);
        e2 = numpy.ones(shape=(nos[0],1));

        #print("The value of X is");
        #print(X);

        #print("The element ele_s is ");
        #print(ele_s);

        #print("The element e2");
        #print(e2);

        u0 = (numpy.dot(X, ele_s) - e2);
        #print("The value of the element u0");
        #print(u0);

        u1 = q/(numpy.dot(numpy.transpose(ele_s), ele_x));
        print("The value of the element u1");
        print(u1);


        # Final u ........................

        u = numpy.dot(( I - A_n),u1*u0);
        print("The value of the element u is ")
        print(u);

        # Final d ............................

        d = - beta * numpy.dot(X, u)/numpy.linalg.norm(u);

        print("The value of d is ")
        print(d);

        #break;
        # Primal step

        if (numpy.linalg.norm(u)>= gamma):

            print ("Entering the first condition")
            ele_x = ele_x + d;
            ele_s = ele_s;
            ele_p = ele_p;

        else:
            print ("Entering the second condition")

            ele_x = ele_x;

            p1 = numpy.linalg.inv(numpy.dot(numpy.dot(A, numpy.dot(X, X)), numpy.transpose(A)));
            p2 = numpy.dot(A, X);
            p3 = (numpy.dot(X, ele_s) - numpy.dot(numpy.transpose(ele_s), ele_x)/q * e2);

            ele_p = ele_p + numpy.dot(numpy.dot(p1, p2), p3);

            print ("Entering the second condition")

            s1 = (u + e2);
            s2 = numpy.linalg.inv(X);

            print("The value of s1");
            print(s1);

            print("The value of s2");
            print(s2);

            print("The value of the element s and x")
            print(ele_s);
            print(ele_x);

            s3 = numpy.dot(s2,s1);

            print("The value of s3 is ")
            print(s3);

            ele_s = (1/q)*numpy.dot(numpy.transpose(ele_s), ele_x)* s3;

            #break;


    print("The value is of the ele_x ");
    print(ele_x);
    print("The value is of the ele_p");
    print(ele_p);
    print("The value is of the ele_s");
    print(ele_s);

    print("The value is of the condition ");
    print(numpy.dot(numpy.transpose(ele_s), ele_x));

    if not (ele_s >= 0).all():
        print(" U are making a s BLUNDER ");
        break;

    if not (ele_x >= 0).all():
        print(" U are making a x BLUNDER ");
        break;

        #time.sleep(2);
