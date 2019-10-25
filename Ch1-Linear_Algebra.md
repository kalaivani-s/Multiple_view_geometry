
# Chapter 1 - Linear Algebra

## 1. Vector spaces

### 1.1 Introduction
- A set $V$ is called a vector space over field $R$ if it is closed under summation and scalar multiplication, $\textit{i.e.,}$ $\alpha v_1 + \beta v_2 \in V\ \forall v_1, v_2 \in V$ and $ \alpha, \beta \in R$.
- A set $W \subset V$ is called a subspace if $0 \in W$ and $W$ is closed under summation and scalar multiplication
- **Span(S)** - The subspace formed by linear combination of a set $S$ of vectors $\{v_1, v_2, \ldots, v_k\} \in V$
$$Span(S) = \{ v \in V | v = \sum_{i=1}^k \alpha_i v_i \}$$
- Set $S = \{v_1, v_2, \ldots, v_k\}$ is **linearly independent** when $\sum_{i=1}^k \alpha_i v_i = 0 \iff \alpha_i = 0\ \forall i$, $\textit{i.e.,}$ none of the vectors in $S$ can be expressed as a linear combination of the other vectors in $S$.
- **Basis vectors** - maximal set of linearly independent vectors $S = \{v_1, v_2, \ldots, v_k\}$ that span the vector space $Span(S)$

### 1.2 Properties of basis
A vector space $V$ can have infinitely many bases. For any two bases $B$ and $B'$,
- $B$ and $B'$ contain same number of vectors called dimension of $V$
- Any vector $v \in V$ can be uniquely expressed as a linear combination of basis vectors in $B$, $\textit{i.e.,}$ $ v = \sum_i=1^k \alpha_i v_i$.
- Basis vectors in $B'$ can be expressed as a linear combination of basis vectors in $B$, $\textit{i.e.,}$ $b'_i = \sum_j=1^k \alpha_{ij} v_j$. Coefficients $\alpha_{ij}$ form the *basis transformation matrix A*.

### 1.3 Inner product or Dot product or Scalar product
a) Inner product on vector space $V$ is defined as $<.,.> : V \times V \in R$ such that
- $<u, \alpha v + \beta w> = \alpha <u,v> + \beta <u,w>$ (linear)
- $<u,v> = <v,u>$ (symmetric)
- $<v,v>\ \ge\ 0$ and $<v,v> = 0 \iff v = 0$ (positive definite)

b) **Norm**: Length of vector $v$, $|v| = \sqrt<v,v>$

c) **Metric**: For measuring lengths and distances, $d(v,w) = |v-w| = \sqrt<v-w, v-w>$

d) Two vectors $v$ and $w$ are orthogonal $\iff <v,w> = 0$. (Basis vectors do not have to be necessarily orthogonal.)

e) For $ V \in R^n$ and canonical basis $B = I_n$, the inner product is defined as $<x,y> = x^Ty = \sum_{i=1}^n x_i y_i$. L2-norm or Euclidean norm is given by $<x,x> = x^Tx = \sqrt(x_1^2 + \ldots + x_n^2)$ 

***Note***
- A vector space with a metric (to measure lengths / distances) is called *metric space*.
- A vector space whose metric is induced by inner product is called *Hilbert space*.

### 1.4. Kronecker product and stack of matrix
1) **Kronecker product** of two matrices $A \in R^{m \times n}$ and $B \in R^{k \times l}$ yields a matrix $C \in R^{mk \times nl}$ where $B$ is multiplied by each coefficient of $A$ and stacked.

2) **Stack** of matrix $A \in R^{m \times n}$ yields $C \in R^{mn}$ by stacking columns of $A$ vertically

## 2. Linear transformations and matrices

- Linear transformation $L$ between two vector spaces $V$ and $W$ is a map $L : V \rightarrow W$ such that linearity (summation and scaling) holds.
$$ L(x+y) = L(x) + L(y) \forall x,y \in V $$
$$L(\alpha x) = \alpha L(x)$$
- Linear transformation is defined by a matrix $A$ which consists of $L$ applied to all basis vectors, $\textit{i.e.,}\ A =[L(e_1), \ldots, L(e_n)] \in R^{m \times n}$
- Since linear transformations are represented as matrices, linear algebra studies properties of matrices
- Unlike inner product of vectors which doesn't yield another vector, product of matrices yield another matrix

### 2.1 Groups
A set of linear transformations form a group $G$ with an operation $\circ : G \times G \rightarrow G$ such that
- application of operation results in a matrix within the group, $g_1 \circ g_2 \in G\ \forall g_1, g_2 \in G$
- associativity holds, $(g_1 \circ g_2) \circ g_3 = g_1 \circ (g_2 \circ g_3)\ \forall g_1, g_2, g_3 \in G$
- an identity transformation $e$ exists, $\exists e \in G, e \circ g = g \circ e = g\ \forall g \in G$
- inverse transformation $g^{-1}$ exists, $\exists g^{-1} \in G, g^{-1} \circ g = g \circ g^{-1} = e\ \forall g \in G$

All invertible or non-singular real matrices, $A$, have $det(A) \neq 0$ and form a group (general linear group) under matrix multiplication.

Groups can also have a matrix representation.

### 2.2 Affine group A(n)
Affine transformation $L : R^n \rightarrow R^n$ is defined by a matrix $A$ such that $L(x) = Ax + b$ where $A$ must be invertible. Note that this is not a linear transformation unless $b=0$. 

More precisely, represent $x$ as homogeneous coordinates, $x \in R^{n+1}$, to make $L$ a linear map, such that,
$$ 
L ( \left[ {\begin{array}{c}
   x \\
   1 \\
  \end{array} } \right] ) = 
  \left[ {\begin{array}{cc}
   A & b \\
   0 & 1 \\
  \end{array} } \right]  
  \left[ {\begin{array}{c}
   x \\
   1 \\
  \end{array} } \right]  
$$

and $L : R^{n+1} \rightarrow R^{n+1}$

### 2.3 Orthogonal group R(n)
A matrix $ R \in \mathcal{M}(n) $ is called orthogonal if it preserves the inner product ($\mathcal{M}(n)$ denotes all $n \times n$ real matrices).
$$ <Rx,Ry>\ =\ <x,y>\ \forall x,y \in R^n$$
$$ <Rx,Ry>\ =\ xR^TRy\ =\ <x,y>\ \forall x,y \in R^n$$
Hence, we must have $R^TR = RR^T = I$. Since, the scalar product of the transformed vectors equals the scalar product of the actual vectors, orthogonal transformations preserves the angles. Note that scalar product between two vectors gives the angle between the vectors. 

Orthogonal transformations generally represent rotations and mirroring.

### 2.4 Euclidean group E(n)
A Euclidean transformation $E$ is defined by a rotation $R$ (orthogonal matrix) and a translation $T \in R^n$ (vector).
$$E(x) = Rx + T$$
In homogeneous coordinates, 
$$ 
E(n) = \left[ {\begin{array}{cc}
   R & T \\
   0 & 1 \\
  \end{array} } \right]
$$

### 2.5 Range and Null space
- Range or span of a transformation $A$ is defined as subspace of $R^m$ that can be reached by $A$. Range is given by the span of the column vectors of $A$.
- Null space or kernel of $A$ are the subset of vectors $\in R^n$ that are mapped to $0$ by $A$.
- Range and kernel are complementary. $dim(range(A)) + dim(kernel(A)) = n$.
- Useful to study solution of linear equations $Ax = b$
    * $Ax = b$ has a solution only if $b \in range(A)$
    * It has unique solution only if $kernel(A) = {0}$, $\textit{i.e.,}$ only zero vector is in the null space
    * If $kernel(A)$ has other vectors besides the zero vector, then there are many solutions for $Ax = b$

### 2.6 Rank of a matrix
Rank of a matrix $A \in R^{m \times n}$
- $Rank(A) = dim(range(A))$
- $Rank(A) = n - dim(kernel(A))$
- $ 0 \leq rank(A) \leq min(m,n)$
- $Rank(A)$ = maximum number of linearly independent row (or column) vectors of $A$

### 2.7 Eigenvalues and Eigenvectors


## 3. Singular value decomposition


```

```
