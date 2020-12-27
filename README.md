---
author:
- Group 6
title: "**GTDMO - 2nd Assignment**"
---

# Introduction to the Problem

A survey is submitted to 15 customers of a restaurant asking them if 5
different aspects of the restaurant should be improved, namely

-   M1 = variety of starters

-   M2 = variety of first dishes

-   M3 = variety of cakes

-   M4 = speed of service

-   M5 = availability of parking

The collected data are stored in the file data6.csv. Higher grades
correspond to an advice of bigger improvement, while low grades mean
that the customer is satisfied of the present situation.

The collected data are represented in the following table

::: {.center}
        M1   M2   M3   M4   M5
  ---- ---- ---- ---- ---- ----
   1    0    0    1    13   7
   2    8    8    2    0    2
   3    12   7    14   0    1
   4    0    2    1    15   11
   5    0    0    0    13   9
   6    1    0    1    9    10
   7    0    0    0    14   14
   8    9    10   10   1    1
   9    0    1    0    11   10
   10   0    0    0    10   6
   11   1    0    0    13   14
   12   0    1    0    11   11
   13   0    0    2    7    11
   14   13   7    8    0    1
   15   4    12   15   0    0
:::

1.  Can you interpret the results in terms of "concepts" behind the
    evaluation?

2.  Can you group and classify the customers with respect to their
    attention to such concepts?

# Solution

In order to answer the problem we can start running a preliminary
analysis drawing an heatmap and looking for clusters in the columns of
the dataset to better visualize the data:

As the clusterization suggests, we can easily notice that the customers
seems to be split in 2 different groups:

-   First group suggests strongly to improve the aspects regarding the
    speed of service and the availability of parking (M4 and M5).

-   Second group conveys that is the variety of dishes, that needs more
    improvement (M1, M2 and M3).

Looking at some descriptive statistics we can highlight that

-   on average the suggestion of improvement is of 5 points (median)

-   the point that the clients have given more often is 0 (moda)

-   the criteria 1 to 3 have lower necessity of improvement compared to
    the global average.

-   Instead, M4 and M5 are the aspects that are required of major
    improvements.

As shown in the table below:

## Part 1

In order to answer the questions we begin by computing Singular Value
Decomposition (SVD) of our matrix in *R*. The SVD technique decomposes a
matrix M of rank(n) into three components, three new matrices of
rank(n):

$$M = U\Sigma V^T$$

where:

-   M is a nxd matrix

-   $\Sigma$ is a diagonal $r\times r$ matrix containing the singular
    values of M, of length min(n, p), in a decreasing order.

-   U is an orthogonal $nxr$ matrix whose columns contain the left
    singular vectors of M, present if nu $>$ 0.

-   $V^T$ is an orthogonal $rxd$ matrix whose columns contain the right
    singular vectors of M, present if nv $>$ 0. The original V matrix
    has, of course, dimension $dxr$.

The rank of our original Matrix $M$ is 5 as we can easily discover
thanks to a few lines of code in *R*.\
We are initially interested just in $\Sigma$ among the three Matrices we
decomposed $M$ in. $\Sigma$ is a diagonal Matrix containing the Singular
Values on its main diagonal and $0$ everywhere else. The Singular Values
represent the weight of the \"concepts\" of the whole information of a
Matrix but when these weights are small enough they can be ignored and
cut off to both reduce the matrix size and to highlight the remaining
ones. From the preview obtained in the heatmap (Figure 1), we can
already guess what the results of our $\Sigma$ will be. $$\Sigma = 
\begin{bmatrix}
$50.045$ & $0$ & $0$ & $0$ & $0$ \\ 
$0$ & $36.444$ & $0$ & $0$ & $0$ \\ 
$0$ & $0$ & $9.879$ & $0$ & $0$ \\ 
$0$ & $0$ & $0$ & $6.600$ & $0$ \\ 
$0$ & $0$ & $0$ & $0$ & $6.254$ \\ 
\end{bmatrix}$$

As a matter of fact, we can state that just the first $2$ values of them
seem to account for more than 96% of the information energy contained in
our dataset as it can be better visualized in the figures below
$$\frac{\sum_{j=1}^2{\sigma_j^2}}{\sum_{j=1}^5{\sigma_j^2}} = \frac{50.04^2 + 36.44^2}{50.04^2 + 36.44^2 + 9.87^2 + 6.6^2 + 6.25^2}  = 0.96$$

In Figure 2 it is possible to observe the graph rapidly tending to $0$
right after the second singular value, which is exactly what we seek in
SVD.

From Figure 3 we can observe how the first singular value alone capture
more than half of the total information energy while from the third on,
every jump is always ticker than the previous one. Note that in Figure
$3$, with \"Cumulative sum of sigmas over the sum of sigmas\" it is
actually meant $$\frac{\sum_{j=1}^r{\sigma_j}}{\sum_{j=1}^m{\sigma_j}}$$

Truncating the Singular Values and their corresponding columns in U and
rows in $V^T$ after 2 positions, we can obtain 3 new matrices
$\tilde{U}$, $\tilde{\Sigma}$, $\tilde{V}$ of rank 2 which are the best
approximation of M with a rank 2.

This can be formalized in the Eckart-Young Theorem which, in fact,
states that the best absolute approximation to the Matrix M that has
rank r is defined as:
$$\underset{\tilde{M} s.t. rank(\tilde{M})=r}{\arg\min} ||M - \tilde{M}||_F$$
which it turns out to be exactly
$\tilde{U}$$\tilde{\Sigma}$$\tilde{V^T}$.

The notation in the previous equation refers to the Frobenius norm,
which is the square root of the summation of the squares of the
differences between the individual matrix entries. It can be represented
by the equation:

$$||A-B||_F = \sqrt{\sum_{ij}(A_{ij}-B_{ij})^2}$$

This is how we are able to decompose a matrix into lower rank matrices
without losing much of the relevant data.

$$M \approx \tilde{U} \tilde{\Sigma} \tilde{V}^T$$

We can now conclude that keeping just the first 2 singular values allow
us to retain more than 96% of the informational energy of our data set
but reducing the rank of M from 5 to 2. What we are left with, then, are
the equations (8) and (9). The $\tilde{U}$ matrix contains the left
singular vectors representing the *row to concept* similarity matrix
while the $\tilde{V}$ matrix containing the right singular vectors
represents the *columns to concept* similarity matrix. $\tilde{\Sigma}$
contains of course the Singular Values, representing the weight of the
*concepts*. $$\tilde{M} = \tilde{U} \tilde{\Sigma} \tilde{V}^T$$
$$%\begin{align*}
\tilde{M} =
\begin{bmatrix}
$-0.286$ & $-0.040$ \\ 
$-0.054$ & $0.267$ \\ 
$-0.066$ & $0.524$ \\ 
$-0.371$ & $-0.025$ \\ 
$-0.310$ & $-0.061$ \\ 
$-0.268$ & $-0.016$ \\ 
$-0.391$ & $-0.073$ \\ 
$-0.074$ & $0.449$ \\ 
$-0.296$ & $-0.041$ \\ 
$-0.227$ & $-0.045$ \\ 
$-0.378$ & $-0.055$ \\ 
$-0.309$ & $-0.043$ \\ 
$-0.252$ & $-0.009$ \\ 
$-0.056$ & $0.434$ \\ 
$-0.051$ & $0.494$ \\ 
\end{bmatrix} 
\begin{bmatrix}
$50.045$ & $0$ \\ 
$0$ & $36.444$ \\ 
\end{bmatrix}%
\begin{bmatrix}
$-0.069$ & $-0.080$ & $-0.088$ & $-0.737$ & $-0.662$ \\ 
$0.549$ & $0.525$ & $0.635$ & $-0.126$ & $-0.065$ \\ 
\end{bmatrix}
%\end{align*}$$

Now we can finally interpret the results in terms of \"concepts\"
dividing the 5 initial questions of our survey in 2 main groups:

1.  the first group, i.e. the group whose concern is greatest, and
    therefore the one the restaurant should be most aware of, is mainly
    interested with what can be conceptualized as *services quality*:

    -   M4 = speed of the service

    -   M5 = parking availability

2.  The second group is composed by the first 3 aspects, respectively:

    -   M1 = variety of starters

    -   M2 = variety of first dishes

    -   M3 = variety of cakes

    These 3 aspects seems instead to have a slightly minor weight in
    customer's concerns. They can be summarized in the concept of
    *courses variety* .

## Part 2

The second question can be addressed just by looking at the $\tilde{U}$
matrix that, as stated above, proposes the strength of the similarity
between customers and the concepts.

$$\tilde{U} =\begin{bmatrix}
$-0.286$ & $-0.040$ \\ 
$-0.054$ & $0.267$ \\ 
$-0.066$ & $0.524$ \\ 
$-0.371$ & $-0.025$ \\ 
$-0.310$ & $-0.061$ \\ 
$-0.268$ & $-0.016$ \\ 
$-0.391$ & $-0.073$ \\ 
$-0.074$ & $0.449$ \\ 
$-0.296$ & $-0.041$ \\ 
$-0.227$ & $-0.045$ \\ 
$-0.378$ & $-0.055$ \\ 
$-0.309$ & $-0.043$ \\ 
$-0.252$ & $-0.009$ \\ 
$-0.056$ & $0.434$ \\ 
$-0.051$ & $0.494$ \\ 
\end{bmatrix}$$

We can further analyze this classification by computing the
*\"similarity to the concepts\" score* of every single user and
comparing them together. If we represent the customer with her/his new
surveys preferences vector $\Vec{m}$, we can compute a score vector as

$$\Vec{s}_j = \Vec{u}_j\tilde{V}$$

where $\Vec{u}_j$ is the vector related to customer $j$ and $V$ is the
matrix of our SVD mapping the aspects of the survey to the concepts.
Note that we could do this for every new customer and check for their
affiliation to the concepts. If we want this for every customer, we can
simply multiply

$$S = \tilde{M}\tilde{V}$$

obtaining

Where the $j-th$ row represent the score vector of user $j$. If we
assign a number to every customer we can easily classify them in 2
groups:

1.  Users $\{1,4,5,6,7,9,10,11,12,13\}$ are heavy related to the concept
    represented by the first singular value, therefore are more
    concerned with improving *services quality* of the restaurant.

2.  Users $\{2,3,8,14,15\}$ are instead highly related to the concept
    represented by the second singular value and they are indeed more
    interested in improving the *courses variety* of the restaurant.

## Appendix

For the sake of Reproducibility Research, here I attach the full
commented R code used to obtain all Figures, Tables, equations and
results.\

``` {.r language="R"}
```
