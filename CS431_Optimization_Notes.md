# CS 431 — Optimization: Theory and Algorithms
### Step-by-step study notes (Dr. Divya Padmanabhan, IIT Goa)

> **How these notes work**
> - Every lecture is one self-contained block between `<!-- LECTURE START: id -->` and `<!-- LECTURE END: id -->`, so blocks can be cut and moved without breaking anything.
> - Lecture numbers are filled in once the correct order is fixed (until then: *Lecture ?*).
> - Each lecture is explained from the basics: every symbol is defined and every derivation is written line by line, including steps skipped in the handwritten notes.
> - 🔍 **Missing step** = something not written in the original notes, added here. ⚠️ = common confusion. ✅ = key takeaway.
> - Maths is written in LaTeX (`$...$`). Open this file in a Markdown viewer that renders maths (VS Code preview, Obsidian, Typora, GitHub).

## Table of contents

- [Lecture 1: Introduction: what optimization is, modelling, classes of problems, easy vs hard](#lecture1)
- [Lecture 2: Linear Algebra Review (Aug 3 class: first half only)](#la-aug3)
- [Lecture 3: Linear Algebra Review (Part 1): vectors, subspaces, rank, eigenvalues, definiteness, norms](#la)
- [Lecture 4: Calculus Review: limits, partial derivatives, gradient, contours, chain rule, Jacobian, directional derivative](#calculus-aug7)
- [Lecture 5: Calculus Review (Part 2): Hessian, continuity, Big-O / small-o, Taylor series](#aug10)
- [Lecture 6: Unconstrained Optimization: problem setup, FONC & SONC (partial notes)](#unconstrainedopt)
- [Lecture 7: Basic Topology: open/closed/bounded/compact sets, sup & inf, Weierstrass theorem](#topology)
- [Lecture 8: Optimization Problems & Optimality Conditions (FONC, SONC, SOSC)](#aug17)
- [Lecture 9: Linear Regression via Optimality Conditions (Normal Equations)](#linearreg)
- [Lecture 10: One-Dimensional Optimization: bisection search (Aug 21 class, partial)](#aug21)
- [Lecture 11: One-Dimensional Optimization: bisection, Newton's method, root finding, secant method](#onedim)
- [Lecture 12: Descent Methods: descent directions, gradient & Newton directions, step sizes, stopping rules](#apr28-uncon)
- [Lecture 13: Steepest Descent with Exact Line Search on a Quadratic](#graddescent-exactls)
- [Lecture 14: Rates of Convergence: definitions only (partial notes)](#rates-of-convergence)
- [Lecture 15: Rates of Convergence (full lecture): order $p$, factor $\gamma$, the limit-ratio test, examples](#sep2)
- [Lecture 16: Convergence Results for Gradient Descent and Newton's Method](#convergenceresults)
- [Lecture 17: Backtracking Line Search (Armijo condition)](#backtracking-line-search)
- [Lecture 18: Convexity (Part 1): convex sets, convex hull, convex functions, why convex optimization is "easy"](#convexitysep7)
- [Lecture 19: Convexity (Part 2): Jensen's inequality, first- & second-order tests, examples, convexity-preserving operations, epigraph](#convexity-sep-9)
- [Lecture 20: Convex Hulls, Carathéodory's Theorem, and "every problem is convex"](#sep11)
- [Lecture 21: Separating Hyperplane Theorem (partial notes)](#sephyperplane)
- [Lecture 22: Linear Programming: standard form and applications (max flow, transportation, diet)](#oct7)
- [Lecture 23: Separating Hyperplane Theorem and Farkas' Lemma (theorem of alternatives)](#fl)
- [Lecture 24: LP Solutions: standard form, slack/surplus/free variables, basic (feasible) solutions, Fundamental Theorem](#lp-bfs)
- [Lecture 25: Linear Programs: the geometrical (graphical) method](#lp-geometric)
- [Lecture 26: Basic Feasible Solutions = Extreme Points (with examples and proof)](#bfs-oct12)
- [Lecture 27: Proof of the Fundamental Theorem of Linear Programming](#oct14-fundatheorem)
- [Lecture 28: LP Duality: introduction (Oct 14 class, partial)](#oct14)
- [Lecture 29: LP Duality: deriving the dual, primal–dual pairs, weak duality, transportation & diet duals](#oct16)
- [Lecture 30: Strong Duality Theorem for LPs (with proof via Farkas)](#dualitytheoremlp)
- [Lecture 31: Complementary Slackness, certificates of optimality, and fractional knapsack](#compslack-lp)
- [Lecture 32: Sensitivity Analysis: dual variables as shadow prices](#sensitivity-lp)
- [Lecture 33: Constrained Optimality (Part 1): feasible & descent directions, active constraints, toward Lagrange](#co-optimality)
- [Lecture 34: General Constrained Optimization: examples, special classes (LP, QP, QCQP), eliminating constraints](#co-oct26)
- [Lecture 35: Constrained Optimality (Part 2): Lagrange condition for equality constraints, derived via Farkas](#oct28)
- [Lecture 36: KKT Conditions: statement, regular points, worked examples (and why regularity matters)](#nov2)
- [Lecture 37: KKT Examples: KKT points need not be minimizers, and minimizers need not be KKT points](#nov4)
- [Lecture 38: KKT for Convex Problems: Slater's condition and KKT sufficiency](#nov4-convex)
- [Lecture 39: Convex KKT Applications: entropy maximisation, maximum likelihood, LP via KKT, projected gradient descent](#nov5)
- [Lecture 40: Practical Algorithms for Constrained Optimization: projected gradient, feasible directions (SLP), Frank–Wolfe](#nov6)
- [Lecture 41: Constrained Optimization Algorithms I (Yinyu Ye, Stanford MS&E211 #12): gradient projection, feasible directions, reduced gradient, Newton on KKT](#yinyuye-lecture12)
- [Lecture 42: The Frank–Wolfe Algorithm (reading: S. Pokutta, "The Frank–Wolfe algorithm: a short introduction")](#frankwolfe)
- [Lecture 43: Lagrangian Duality (Part 1): dual function, dual problem, examples, LP & QP duals, a duality gap](#duality1)
- [Lecture 44: Lagrangian Duality (Part 2): strong duality for convex problems, KKT from strong duality, projection onto a half-space](#duality2)
- [Lecture 45: Algorithms Using the Dual: dual ascent and dual decomposition](#dualascent)
- [Lecture 46: Duality Theory of Constrained Optimization (reading: R. M. Freund, MIT, 2004)](#lec18-duality-thy)
- [Lecture 47: Newton's Method for Equality-Constrained Optimization](#newton-equality)
- [Lecture 48: ADMM: Alternating Direction Method of Multipliers (reading: R. Tibshirani, CMU 10-725, Lecture 22)](#lecture-22)
- [Lecture 49: Dual Methods (full notes): dual ascent → method of multipliers → ADMM](#dualmethods)

---

<!-- LECTURE START: lecture1 -->
<a id="lecture1"></a>

## Lecture 1: Introduction: what optimization is, modelling, classes of problems, easy vs hard

*Source file:* `Lecture1.pdf`

### 1. What this lecture is about

This is the opening lecture of **CS 431, Optimization: Theory and Algorithms**. It answers:
1. Where do optimization problems come from? (daily life, engineering, ML)
2. What does a mathematical optimization problem look like? (objective, constraints, variables)
3. How do we **model** a word problem as maths? (max flow, vertex cover, mess location; the slides only pose these, the formulations are written out below)
4. What are the main **classes** (LP, QP, integer, SDP, …)?
5. Which problems are **easy** (polynomial time) vs **hard** (NP-hard)? Is convex always easy?

---

### 2. Background you need first

- A **graph** $G=(V,E)$: a set of nodes $V$ joined by edges $E$ (roads, pipes, links).
- $|V|$ = number of nodes. The number of subsets of $V$ is $2^{|V|}$ (each node is either in or out).
- "Polynomial time" = running time grows like $n^k$ for some fixed $k$ (fast). "NP-hard" = no polynomial-time algorithm is known, and finding one would solve thousands of other hard problems.

---

### 3. Main content

#### 3.1 Optimization in daily life

- **Planning a trip** (a week in Berlin): *minimise* travel cost, *maximise* places covered.
- **Choosing a university:** minimise travel risk; constraints such as top-ten ranking, fees within budget, good opportunities.
- **Preparing for exams** (6 courses): plan time so all subjects are learnt, assignments done, and you still eat and sleep.

Every one of these has **things you choose** (variables), **a goal** (objective), and **rules** (constraints).

#### 3.2 What the course covers

How optimization problems arise in engineering; modelling them mathematically; solution techniques; software (MATLAB, Python); recognising easy vs hard problems; approximating hard ones.

**Examples:** recommendation engines (maximise chance of purchase), influence maximisation (choose $k$ people within budget to spread info to as many as possible), cab services (match drivers to trips: driver near source, no trip unmatched, short routes).

#### 3.3 Standard form

$$\min_x\ f(x)\qquad\text{s.t.}\quad x\in\Omega$$
- $f:\mathbb{R}^n\to\mathbb{R}$: **objective function** (single objective in this course),
- $\Omega\subseteq\mathbb{R}^n$: **constraint set / feasible set**,
- $x$: **decision variable**,
- $x^*=\arg\min f(x)$: the minimiser.

#### 3.4 Modelling example 1: Maximum flow

*Problem:* send as much oil as possible from source $S$ to sink $T$ in a network $G=(V,E)$. Each link $(i,j)$ has capacity $u_{ij}$; nothing can be stored at intermediate nodes.

🔍 *Formulation (not given on the slides).*
- Variables: $x_{ij}$ = amount of oil sent along edge $(i,j)$.
- Objective: total flow leaving $S$: $\max\sum_{j:(S,j)\in E}x_{Sj}$.
- Capacity: $0\le x_{ij}\le u_{ij}$ for every edge.
- Flow conservation (no storage): for every node $v\neq S,T$: $\sum_{i:(i,v)\in E}x_{iv}=\sum_{j:(v,j)\in E}x_{vj}$ (in = out).

All functions are linear, so this is a **Linear Program**.
**Certificate of optimality:** any **cut** (a set of edges whose removal separates $S$ from $T$) has total capacity that is an **upper bound** on the max flow, since all oil must cross it. In the slides' example a flow of 11 is found and a cut of value 11 exists, so 11 is optimal (max-flow = min-cut).

#### 3.5 Modelling example 2: Vertex cover ("safe IIT campus")

*Problem:* place guards at intersections so that every road is watched (a guard at a junction sees all roads meeting there); minimise the number of guards.

🔍 *Formulation.*
- Variables: $x_v\in\{0,1\}$ for each intersection $v$ (1 = guard placed).
- Objective: $\min\sum_{v\in V}x_v$.
- Constraints: every road $(u,v)\in E$ must be covered: $x_u+x_v\ge1$.

This is an **integer program**. The slides find 6 guards and ask: is it the best? Here, yes. But for a general network, brute force must check up to $2^{|V|}$ subsets. That is hopeless for a town (with $|V|=100$, $2^{100}\approx10^{30}$).

#### 3.6 Modelling example 3: Mess location (exercise from the slides)

5 hostels at locations $h_1,\dots,h_5\in\mathbb{R}^2$ with $n_j$ students; 2 messes with capacities $c_1,c_2$. Choose mess locations and an allocation that minimises the total distance students travel.

🔍 *A formulation.*
- Variables: mess locations $m_1,m_2\in\mathbb{R}^2$; $y_{ij}\ge0$ = number of students from hostel $j$ assigned to mess $i$.
- Objective: $\min\sum_{i=1}^2\sum_{j=1}^5y_{ij}\,\|m_i-h_j\|_2$.
- Every student is assigned: $\sum_{i=1}^2y_{ij}=n_j$ for each $j$.
- Capacity: $\sum_{j=1}^5y_{ij}\le c_i$ for each $i$ (the slides say "a total of $c_i$ can go"; use $=c_i$ if exactly).
- (Optionally $y_{ij}$ integer.)

Note the product $y_{ij}\|m_i-h_j\|$ of two variables: this is a **non-convex** problem, even though each piece looks simple.

#### 3.7 Classes of optimization problems

- **Unconstrained:** $\Omega=\mathbb{R}^n$.
- **Constrained:** $\Omega\subset\mathbb{R}^n$, usually
$$\Omega=\{x:\ g_i(x)\le0,\ i=1,\dots,m;\ \ h_j(x)=0,\ j=1,\dots,l\}.$$
- By function type:
  - **Linear Programming (LP):** everything linear, e.g. $\min\sum_ic_ix_i$ s.t. $Ax=b$ (usually with $x\ge0$).
  - **Quadratic Programming (QP):** at least one function is quadratic.
- By variable domain:
  - **Discrete / Integer Programming:** some variables must be integers (e.g. vertex cover).
  - **Semidefinite (SDP)** and **second-order cone programming (SOCP)**: variables are matrices constrained to be PSD, or vectors in a "cone".

#### 3.8 Easy vs hard

- Max flow: solvable in polynomial time (**easy**).
- Vertex cover: **NP-hard**, unlikely to be polynomial.
- Strategy: if in P, find efficient algorithms; if not, find good **approximations** or special solvable cases.

**Is discrete easier than continuous?** A discrete feasible set is countable (finite), a continuous one uncountable. Yet the answer is **NO**:

| Class | Easy examples | Hard examples |
|---|---|---|
| Discrete / combinatorial | Linear assignment, network flow | Knapsack, vertex cover |
| Continuous | Linear programs, SOCP | Polynomial optimization |

**Is convex always easy?** Traditionally believed so, but **not in general**. LP, SOCP, SDP are polynomial-time; **copositive optimization** is convex but hard. Convexity is still special: it has very desirable properties (local = global; see the Convexity blocks), efficient methods for many classes, and it can even help solve some discrete problems (relaxations).

#### 3.9 Research venues (for interest)

SIAM J. Optimization, Mathematical Programming, Operations Research, Math of OR, INFORMS journals; applied: Management Science, Transportation Science, Networks, NeurIPS; preprints: Optimization Online, arXiv.

---

### 4. Summary

✅ An optimization problem: decision variables $x$, objective $f$, feasible set $\Omega$.
✅ Modelling = turning words into variables, objective and constraints (max flow = LP; vertex cover = integer program).
✅ Classes: unconstrained/constrained; LP, QP, IP, SDP, SOCP.
✅ Easy/hard does not split along discrete vs continuous, nor exactly along convex vs non-convex.

<!-- LECTURE END: lecture1 -->

---

<!-- LECTURE START: la-aug3 -->
<a id="la-aug3"></a>

## Lecture 2: Linear Algebra Review (Aug 3 class: first half only)

*Source file:* `LA-Aug3.pdf`

> 📌 **This file is the first half of the full Linear Algebra notes (`LA.pdf`).** Everything in it is explained in detail, with all missing steps, in the block **"Linear Algebra Review (Part 1)"**. If you keep both blocks, read that one. This block is a short checklist of what was covered on Aug 3.

### What was covered on this day

1. **Linear combination:** $\alpha_1v_1+\alpha_2v_2+\alpha_3v_3$ for vectors $v_i\in\mathbb{R}^n$ and scalars $\alpha_i\in\mathbb{R}$. (Example vectors $(1,-1,2)^T,(2,4,-3)^T,(7,10,6)^T$.)
2. **Linear independence:** $\sum_i c_iv_i=0\Rightarrow c_i=0\ \forall i$. For two vectors this just means neither one is a multiple of the other.
3. **Subspace** $S\subseteq\mathbb{R}^n$: closed under scaling ($\alpha v\in S$) and addition ($v_1+v_2\in S$). Example: $S=\{v\in\mathbb{R}^3:v_1=0\}$.
4. **Span, basis, dimension:** span = all linear combinations. Basis = smallest independent spanning set, e.g. $B_1=\{(0,1,0)^T,(0,0,1)^T\}$ or $B_2=\{(0,1,0)^T,(0,1,1)^T\}$ for the $S$ above. Dimension = size of any basis ($=2$ here).
5. **$A\in\mathbb{R}^{m\times n}$ as a map** $x\in\mathbb{R}^n\mapsto b=Ax\in\mathbb{R}^m$.
6. **Four subspaces:** $\mathcal N(A)$, $\mathcal N(A^T)$, column space $\mathcal R(A)$, row space $\mathcal R(A^T)$. Also: $\mathcal N(A)\perp\mathcal R(A^T)$ inside $\mathbb{R}^n$, and $\mathcal N(A^T)\perp\mathcal R(A)$ inside $\mathbb{R}^m$. Orthogonal means $v^Tz=0$.
7. **Solving $Ax=b$:** 0, 1 or ∞ solutions. Every solution is $x_p+q$ with $Aq=0$. Worked example with row operations $R_2-7R_1$, then $R_2\times(-1/18)$.
8. **Inconsistent example** started at the end of the class: $A=\begin{pmatrix}1&4&3\\2&8&6\end{pmatrix}$, $b=(-2,15)^T$. It is finished in the full LA block: row 2 becomes $0=19$, so there is no solution.

✅ Read the full **Linear Algebra Review (Part 1)** block for the derivations.

<!-- LECTURE END: la-aug3 -->

---

<!-- LECTURE START: la -->
<a id="la"></a>

## Lecture 3: Linear Algebra Review (Part 1): vectors, subspaces, rank, eigenvalues, definiteness, norms

*Source file:* `LA.pdf`

### 1. What this lecture is about

Optimization is mostly done with vectors and matrices. Before we can minimise a function of many variables, we need the "grammar" of linear algebra. This lecture reviews:

1. Linear combinations and linear independence
2. Subspaces, span, basis, dimension
3. The four fundamental subspaces of a matrix (null space, column space/range, row space, left null space)
4. Solving $Ax=b$ (0, 1 or infinitely many solutions) and row reduction
5. Rank and the Rank–Nullity theorem
6. Eigenvalues and eigenvectors (and why eigenvectors of a symmetric matrix are orthogonal)
7. Positive (semi)definite matrices and Sylvester's criterion
8. Vector norms, inner products, Cauchy–Schwarz inequality

**Convention used in the whole course:** every vector is a **column vector**. $x\in\mathbb{R}^n$ means
$$x=\begin{pmatrix}x_1\\x_2\\\vdots\\x_n\end{pmatrix},$$
and $x^T=(x_1,\dots,x_n)$ is the corresponding row vector.

---

### 2. Background you need first

- $\mathbb{R}$ = the set of real numbers. $\mathbb{R}^n$ = all lists of $n$ real numbers (vectors with $n$ entries).
- $A\in\mathbb{R}^{m\times n}$ means $A$ is a matrix with $m$ rows and $n$ columns.
- **Matrix–vector product.** If $A\in\mathbb{R}^{m\times n}$ and $x\in\mathbb{R}^n$, then $Ax\in\mathbb{R}^m$ and its $i$-th entry is
  $$(Ax)_i=\sum_{j=1}^n A_{ij}x_j .$$
  Two equally useful ways to see it:
  - **Row view:** $(Ax)_i = a_i^T x$, where $a_i^T$ is the $i$-th row of $A$.
  - **Column view:** $Ax = x_1 A_{:,1}+x_2A_{:,2}+\dots+x_nA_{:,n}$, i.e. $Ax$ is a *combination of the columns* of $A$ with weights $x_j$. (This view is the key to the "column space".)
- **Transpose:** $(A^T)_{ij}=A_{ji}$. Rules: $(AB)^T=B^TA^T$, $(A^T)^T=A$.
- **Dot product / inner product:** $x^Ty=\sum_i x_iy_i$ (a single number).
- Symbols: $\forall$ = "for all", $\exists$ = "there exists", $\Rightarrow$ = "implies", $\Leftrightarrow$ = "if and only if", $\in$ = "belongs to", $\subseteq$ = "is a subset of".

---

### 3. Main content with full derivations

#### 3.1 Linear combinations

Given vectors $v_1,v_2,\dots,v_m\in\mathbb{R}^n$ and real numbers ("scalars") $\alpha_1,\dots,\alpha_m\in\mathbb{R}$, the vector
$$\alpha_1v_1+\alpha_2v_2+\dots+\alpha_mv_m$$
is called a **linear combination** of the $v_i$.

**Example from the notes ($n=3$):**
$$v_1=\begin{pmatrix}1\\-1\\2\end{pmatrix},\;v_2=\begin{pmatrix}2\\4\\-3\end{pmatrix},\;v_3=\begin{pmatrix}7\\10\\6\end{pmatrix},\quad \alpha_1,\alpha_2,\alpha_3\in\mathbb{R}.$$
🔍 *Missing step — computing one.* Take $\alpha_1=1,\alpha_2=2,\alpha_3=0$:
$$1\cdot\begin{pmatrix}1\\-1\\2\end{pmatrix}+2\begin{pmatrix}2\\4\\-3\end{pmatrix}+0\cdot v_3=\begin{pmatrix}1+4\\-1+8\\2-6\end{pmatrix}=\begin{pmatrix}5\\7\\-4\end{pmatrix}.$$
You multiply each vector entry-by-entry by its scalar and then add entry-by-entry.

#### 3.2 Linear independence

**Two vectors.** $v_1,v_2$ are *linearly dependent* if one is a scalar multiple of the other: $v_2=\alpha v_1$ for some $\alpha$. They are *linearly independent* if this is impossible — i.e. there is **no** $\alpha$ with $v_2-\alpha v_1=0$.

The general definition (works for any number of vectors):
$$c_1v_1+c_2v_2=0\;\Rightarrow\;c_1=c_2=0 .$$
🔍 *Why this is the same thing.* Suppose $c_1v_1+c_2v_2=0$ with $c_2\neq0$. Divide by $c_2$: $v_2=-\frac{c_1}{c_2}v_1$, so $v_2$ is a multiple of $v_1$ (dependent). So "independent" means the *only* way to combine them into zero is the trivial way (all coefficients zero).

**General collection.** $v_1,\dots,v_m$ are linearly independent (l.i.) if
$$\sum_{i=1}^m c_iv_i = 0_n \;\Rightarrow\; c_i=0\ \ \forall i .$$
($0_n$ is the zero vector with $n$ entries.) Meaning: no vector in the collection can be built out of the others.

⚠️ Matrix form: putting the $v_i$ as columns of $V=[v_1\ \cdots\ v_m]$, independence means $Vc=0\Rightarrow c=0$, i.e. the null space of $V$ (defined below) is just $\{0\}$.

#### 3.3 Subspace

A set $S\subseteq\mathbb{R}^n$ is a **subspace** if it is closed under the two vector operations:

1. **Scaling:** $v\in S\Rightarrow \alpha v\in S$ for every $\alpha\in\mathbb{R}$.
2. **Addition:** $v_1,v_2\in S\Rightarrow v_1+v_2\in S$.

🔍 Consequence: taking $\alpha=0$ in (1) shows $0\in S$ always. So any set not containing the origin (e.g. a line not through the origin) is **not** a subspace.

**Example from the notes:** $S=\{v\in\mathbb{R}^3: v_1=0\}$ (all vectors whose first entry is zero), i.e. vectors of the form $\begin{pmatrix}0\\v_2\\v_3\end{pmatrix}$.
- Scaling: $\alpha\begin{pmatrix}0\\v_2\\v_3\end{pmatrix}=\begin{pmatrix}0\\\alpha v_2\\\alpha v_3\end{pmatrix}$ — first entry still $0$ ✓.
- 🔍 Addition (not written): $\begin{pmatrix}0\\v_2\\v_3\end{pmatrix}+\begin{pmatrix}0\\w_2\\w_3\end{pmatrix}=\begin{pmatrix}0\\v_2+w_2\\v_3+w_3\end{pmatrix}$ — first entry $0+0=0$ ✓.
So $S$ is a subspace (geometrically: the $v_2$–$v_3$ plane inside 3-D space).

#### 3.4 Span, basis, dimension

**Span** of $v_1,\dots,v_m$:
$$\operatorname{span}\{v_1,\dots,v_m\}=\Big\{v:\ v=\sum_{i=1}^m c_iv_i,\ c_i\in\mathbb{R}\Big\},$$
the set of *all* linear combinations. A span is always a subspace.

**Basis for $S$:** the smallest linearly independent set of vectors from $S$ that spans $S$.

**Example** $S=\{v\in\mathbb{R}^3:v_1=0\}$. Two different bases:
$$B_1=\Big\{\begin{pmatrix}0\\1\\0\end{pmatrix},\begin{pmatrix}0\\0\\1\end{pmatrix}\Big\},\qquad B_2=\Big\{\begin{pmatrix}0\\1\\0\end{pmatrix},\begin{pmatrix}0\\1\\1\end{pmatrix}\Big\}.$$
🔍 *Check $B_2$ spans $S$:* any $(0,a,b)^T = (a-b)\,(0,1,0)^T + b\,(0,1,1)^T$. Verify: $(0,\ (a-b)+b,\ b)=(0,a,b)$ ✓. And they are independent (neither is a multiple of the other).

The notes also show: $\{(0,1,0)^T,(0,1,1)^T,(0,2,1)^T\}$ is **not** a basis — it spans $S$ but is not independent, because $(0,2,1)^T=(0,1,0)^T+(0,1,1)^T$. Drop the redundant vector and you get the basis $\{(0,1,1)^T,(0,2,1)^T\}$ (also fine).

**Dimension of $S$** = number of vectors in *any* basis of $S$ (all bases have the same size). Here $\dim S=2$.

#### 3.5 A matrix as a map, and its four fundamental subspaces

$A\in\mathbb{R}^{m\times n}$ takes $x\in\mathbb{R}^n$ to $b=Ax\in\mathbb{R}^m$. The picture in the notes: a region of $\mathbb{R}^n$ gets mapped into a (shaded) region of $\mathbb{R}^m$.

| Name | Definition | Lives in |
|---|---|---|
| **Null space** $\mathcal N(A)$ | $\{x\in\mathbb{R}^n: Ax=0\}$ | $\mathbb{R}^n$ |
| **Left null space** $\mathcal N(A^T)$ | $\{x\in\mathbb{R}^m: A^Tx=0\}$ | $\mathbb{R}^m$ |
| **Column space / Range** $\mathcal R(A)$ | $\{b\in\mathbb{R}^m: b=Ax \text{ for some } x\in\mathbb{R}^n\}$ | $\mathbb{R}^m$ |
| **Row space** $=\mathcal R(A^T)$ | $\{b\in\mathbb{R}^n: b=A^Tx\text{ for some }x\in\mathbb{R}^m\}$ | $\mathbb{R}^n$ |

- The null space is **trivial** if $\mathcal N(A)=\{0_n\}$ (only $x=0$ is sent to zero).
- Why is it called "column space"? By the column view, $Ax=\sum_j x_j(\text{column } j)$, so the set of all $Ax$ is exactly the span of the columns.

**Example from notes: $\mathcal R(A)$ is a subspace.** Let $b_1\in\mathcal R(A)$, so $\exists x: b_1=Ax$. Is $2b_1\in\mathcal R(A)$? Yes: $A(2x)=2Ax=2b_1$. 🔍 In general $A(\alpha x)=\alpha b_1$ and $A(x+x')=b_1+b_2$, so both closure rules hold.

**Orthogonality picture (the boxes in the notes).**
$$\mathbb{R}^n=\mathcal N(A)\ \oplus\ \mathcal R(A^T),\qquad \mathbb{R}^m=\mathcal N(A^T)\ \oplus\ \mathcal R(A).$$
Two vectors $v,z$ are **orthogonal** iff $v^Tz=0$. The null space is orthogonal to the row space:

🔍 *Proof (missing in notes).* Let $x\in\mathcal N(A)$, so $Ax=0$, i.e. every row satisfies $a_i^Tx=0$. A row-space vector is $r=A^Ty=\sum_i y_ia_i$. Then
$$r^Tx=\Big(\sum_i y_ia_i\Big)^Tx=\sum_i y_i\,(a_i^Tx)=\sum_i y_i\cdot0=0.$$
This is exactly the picture in the notes: $A$ with rows $a_1^T,\dots,a_m^T$ times $x$ equals $0$ means $x$ is perpendicular to every row.

#### 3.6 Solving $Ax=b$: zero, one, or infinitely many solutions

A linear system has either **0, 1 or $\infty$** solutions — never exactly 2.

**Structure of the solution set.** Suppose $x_p$ is one particular solution ($Ax_p=b$) and $q\in\mathcal N(A)$ ($Aq=0$). Then
$$A(x_p+q)=Ax_p+Aq=b+0=b.$$
So $x_p+q$ is also a solution. Hence:
- If $\dim\mathcal N(A)=0$ (only $q=0$): **0 or 1** solution.
- If $\dim\mathcal N(A)\neq 0$: **0 or $\infty$** solutions (a whole line/plane of $q$'s can be added).

🔍 Conversely every solution has this form: if $Ax=b$ and $Ax_p=b$, then $A(x-x_p)=0$, so $x-x_p\in\mathcal N(A)$.

**Example 1 (row operations).** Solve
$$\begin{pmatrix}1&4&3\\7&10&6\end{pmatrix}\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}-2\\1\end{pmatrix}.$$
*Step 1:* eliminate the 7 under the first pivot: $R_2\leftarrow R_2-7R_1$.
- Row 2 of $A$: $(7-7\cdot1,\ 10-7\cdot4,\ 6-7\cdot3)=(0,-18,-15)$.
- Right side: $1-7(-2)=1+14=15$.

*Step 2:* make the second pivot 1: $R_2\leftarrow R_2\times(-\tfrac1{18})$.
- $(0,1,\tfrac{15}{18})$, right side $15\cdot(-\tfrac1{18})=-\tfrac{15}{18}=-\tfrac56$.

Now the system reads
$$x+4y+3z=-2,\qquad y+\tfrac{15}{18}z=-\tfrac56 .$$
🔍 *Finishing it (not done in the notes).* There are 3 unknowns but only 2 equations, so one variable is **free**. Let $z=t$ (any real). Then
$$y=-\tfrac56-\tfrac56t,\qquad x=-2-4y-3z=-2+\tfrac{10}{3}+\tfrac{10}{3}t-3t=\tfrac43+\tfrac13t.$$
So
$$\begin{pmatrix}x\\y\\z\end{pmatrix}=\underbrace{\begin{pmatrix}4/3\\-5/6\\0\end{pmatrix}}_{x_p}+t\underbrace{\begin{pmatrix}1/3\\-5/6\\1\end{pmatrix}}_{q\in\mathcal N(A)} ,$$
infinitely many solutions, exactly the "$x_p+q$" structure. (Check $q$: $1/3+4(-5/6)+3=1/3-10/3+3=0$ ✓.)

**Example 2 (inconsistent).**
$$A=\begin{pmatrix}1&4&3\\2&8&6\end{pmatrix},\quad b=\begin{pmatrix}-2\\15\end{pmatrix}.$$
Row 2 of $A$ is exactly $2\times$ row 1. Do $R_2\leftarrow R_2-2R_1$ on the augmented matrix $[A\mid b]$:
$$\left(\begin{array}{ccc|c}1&4&3&-2\\0&0&0&15-2(-2)=19\end{array}\right).$$
The second row says $0x+0y+0z=19$, i.e. $0=19$: impossible. **Inconsistent, number of solutions $=0$.** Geometrically $b\notin\mathcal R(A)$.

#### 3.7 Rank and Rank–Nullity

**Rank$(A)$** — three equivalent descriptions:
- number of non-zero rows in the **row-reduced echelon form** RREF$(A)$;
- number of linearly independent rows of $A$;
- number of linearly independent columns of $A$ (row rank = column rank, always).

Since you cannot have more independent rows than rows, or more independent columns than columns: for $A\in\mathbb{R}^{m\times n}$,
$$\operatorname{rank}(A)\le\min(m,n).$$
"Full rank" means equality.

**Rank–Nullity theorem:**
$$\operatorname{rank}(A)+\underbrace{\operatorname{nullity}(A)}_{=\dim\mathcal N(A)}=n\ \ (\text{number of columns}).$$
🔍 *Intuition.* In RREF, each non-zero row has a pivot variable; the remaining variables are free, and each free variable gives one independent null-space direction. So (#pivots) + (#free) = (#variables) = $n$. In Example 1: rank $2$, $n=3$, nullity $1$ (one free variable $t$) ✓.

#### 3.8 Eigenvalues and eigenvectors

For a **square** $A\in\mathbb{R}^{n\times n}$, a non-zero vector $v$ is an eigenvector with eigenvalue $\lambda$ if
$$Av=\lambda v,\qquad v\neq0 .$$
Meaning: $A$ does not rotate $v$, it only stretches it by $\lambda$.

🔍 *How to find them (steps skipped in notes).*
1. $Av=\lambda v\iff Av-\lambda Iv=0\iff (A-\lambda I)v=0$.
2. So $v\in\mathcal N(A-\lambda I)$ and $v\neq0$, so this null space is non-trivial, so $A-\lambda I$ is **not full rank** (not invertible).
3. A square matrix is not invertible iff its determinant is $0$: 
$$\boxed{\ |A-\lambda I|=0\ }\quad(\text{characteristic equation}).$$
Solve this polynomial for $\lambda$, then solve $(A-\lambda I)v=0$ for $v$.

The **eigenspace** of $\lambda$ is $\mathcal N(A-\lambda I)$ (all eigenvectors for $\lambda$, plus $0$). E.g. if $v_1,v_2$ both solve it, eigenspace $=\operatorname{span}\{v_1,v_2\}$.

**Note:** any multiple of an eigenvector is again an eigenvector (e.g. $100v$) — $A(100v)=100Av=\lambda(100v)$.

**Theorem (real symmetric matrices).** If $A=A^T$ is real, eigenvectors from *different* eigenvalues are orthogonal.

*Proof, every step.* Let $Av=\lambda_1v$, $Az=\lambda_2z$, $\lambda_1\neq\lambda_2$. Compute the number $\langle Av,z\rangle=(Av)^Tz$ in two ways.
- Way 1: $(Av)^Tz=(\lambda_1v)^Tz=\lambda_1v^Tz$.
- Way 2: $(Av)^Tz=v^TA^Tz$ (transpose rule $(Av)^T=v^TA^T$) $=v^TAz$ (since $A^T=A$) $=v^T(\lambda_2 z)=\lambda_2v^Tz$.

Equate: $\lambda_1v^Tz=\lambda_2v^Tz\Rightarrow(\lambda_1-\lambda_2)v^Tz=0$. Since $\lambda_1-\lambda_2\neq0$, we must have $v^Tz=0$. ∎

So $\mathbb{R}^n$ splits into mutually orthogonal eigenspaces $W_1,W_2,\dots$ (the stacked boxes in the notes). Also (fact): a real symmetric matrix has all **real** eigenvalues.

#### 3.9 Positive semidefinite (PSD) and positive definite (PD) matrices

For a **symmetric** $A\in\mathbb{R}^{n\times n}$, the expression
$$v^TAv=\sum_{i=1}^n\sum_{j=1}^n A_{ij}v_iv_j$$
is called a **quadratic form** (a single number; $(1\times n)(n\times n)(n\times 1)=1\times1$).

**PSD, written $A\succeq0$** — three equivalent conditions:
1. all eigenvalues $\ge0$;
2. $v^TAv\ge0$ for all $v\in\mathbb{R}^n$;
3. $A=BB^T$ for some $B\in\mathbb{R}^{n\times k}$.

🔍 *Why 3 ⇒ 2:* $v^TBB^Tv=(B^Tv)^T(B^Tv)=\|B^Tv\|_2^2\ge0$.
🔍 *Why 2 ⇒ 1:* if $Av=\lambda v$, $v\ne0$, then $0\le v^TAv=\lambda v^Tv=\lambda\|v\|^2$, and $\|v\|^2>0$, so $\lambda\ge0$.

**PD, written $A\succ0$:** $v^TAv>0$ for all $v\neq0$ ⟺ all eigenvalues $>0$.
- **Negative semidefinite:** $-A\succeq0$. **Negative definite:** $-A\succ0$.
- **Indefinite:** none of the above (some directions give positive, some negative $v^TAv$).

Why we care: later, a point is a minimum of $f$ when the Hessian (matrix of second derivatives) is PSD/PD — this is the matrix version of "second derivative $\ge 0$".

#### 3.10 Sylvester's criterion (testing definiteness with determinants)

- A **principal minor** = determinant of the submatrix obtained by keeping the *same* set of row and column indices (e.g. rows {1,3} & columns {1,3}).
- A **leading principal minor** = keep the first $k$ rows and first $k$ columns, $k=1,\dots,n$.

**Rules:**
- $A\succ0$ (PD) $\iff$ **all leading principal minors $>0$.**
- $A\succeq0$ (PSD) $\iff$ **all principal minors $\ge0$** (⚠️ for PSD you must check *all* principal minors, not only leading ones).

**Example 1 (PSD):** $A=\begin{pmatrix}1&2&-1\\2&4&-2\\-1&-2&1\end{pmatrix}$.
- 1×1 principal minors (diagonal): $1,4,1\ge0$.
- 2×2: rows/cols {1,2}: $\begin{vmatrix}1&2\\2&4\end{vmatrix}=1\cdot4-2\cdot2=0$; {2,3}: $\begin{vmatrix}4&-2\\-2&1\end{vmatrix}=4-4=0$; {1,3}: $\begin{vmatrix}1&-1\\-1&1\end{vmatrix}=1-1=0$.
- 🔍 3×3: $|A|=0$ (row 2 = 2×row 1).
All $\ge0$ ⇒ PSD (but not PD, since some are $0$). Indeed $A=uu^T$ with $u=(1,2,-1)^T$.

**Example 2 (PD):** $A=\begin{pmatrix}4&2&0\\2&3&1\\0&1&2\end{pmatrix}$. Leading principal minors:
- $k=1$: $4>0$.
- $k=2$: $\begin{vmatrix}4&2\\2&3\end{vmatrix}=12-4=8>0$.
- $k=3$: 🔍 expand along row 1: $|A|=4\begin{vmatrix}3&1\\1&2\end{vmatrix}-2\begin{vmatrix}2&1\\0&2\end{vmatrix}+0=4(6-1)-2(4-0)=20-8=12>0$.

All $>0$ ⇒ $A$ is PD. (The notes cross out $\begin{vmatrix}3&1\\1&2\end{vmatrix}$ and $\begin{vmatrix}4&0\\0&2\end{vmatrix}$ — these are principal minors but **not leading**, so they are not needed for the PD test.)

#### 3.11 Vector norms

A norm measures the "length" of a vector. For $x=(x_1,\dots,x_n)^T$:
- **2-norm (Euclidean):** $\|x\|_2=\sqrt{x^Tx}=\sqrt{x_1^2+\dots+x_n^2}$.
- **1-norm:** $\|x\|_1=|x_1|+|x_2|+\dots+|x_n|$.
- **$p$-norm:** $\|x\|_p=\big(|x_1|^p+\dots+|x_n|^p\big)^{1/p}$, $p\ge1$.
- **$\infty$-norm:** $\|x\|_\infty=\max_{1\le i\le n}|x_i|$ (the limit of $\|x\|_p$ as $p\to\infty$).

**Unit balls in 2-D** (the pictures):
- $\{x:\|x\|_2\le1\}$, i.e. $x^Tx\le 1$, i.e. $x_1^2+x_2^2\le1$: a **disc**.
- $\{x:\|x\|_1\le1\}$, i.e. $|x_1|+|x_2|\le1$: a **diamond** with corners $(\pm1,0),(0,\pm1)$.
- $\{x:\|x\|_\infty\le1\}$: the **square** $[-1,1]^2$.
They are nested: diamond ⊂ disc ⊂ square.

#### 3.12 Inner products and the $A$-norm

1. Standard inner product: $\langle x,y\rangle=x^Ty$.
2. If $A\succ0$, then $\langle x,y\rangle_A=x^TAy$ is also an inner product, and it gives the **$A$-norm** $\|x\|_A=\sqrt{x^TAx}$. (Notes write $\|x\|_A = x^TAx\ge0$; strictly, the norm is the square root.) PD-ness is exactly what guarantees $x^TAx>0$ for $x\ne0$.

#### 3.13 Cauchy–Schwarz inequality

Geometrically, $x^Ty=\|x\|_2\|y\|_2\cos\theta$, where $\theta$ = angle between $x$ and $y$. Since $-1\le\cos\theta\le1$:
$$-\|x\|_2\|y\|_2\ \le\ x^Ty\ \le\ \|x\|_2\|y\|_2,\qquad\text{i.e. } |x^Ty|\le\|x\|_2\|y\|_2 .$$
- Upper bound reached when $x,y$ point the **same** way ($\theta=0$).
- Lower bound reached when they point in **opposite** directions ($\theta=\pi$).

🔍 *An algebraic proof (no angles).* For any real $t$, $0\le\|x-ty\|_2^2=x^Tx-2t\,x^Ty+t^2y^Ty$. This is a quadratic in $t$ that is never negative, so its discriminant is $\le0$: $(2x^Ty)^2-4(y^Ty)(x^Tx)\le0$, i.e. $(x^Ty)^2\le\|x\|^2\|y\|^2$. ∎

This inequality is used again and again (e.g. to show that $-\nabla f$ is the steepest-descent direction).

---

### 4. Worked examples (extra practice)

1. **Is $\{(1,2)^T,(2,4)^T\}$ independent?** $(2,4)=2(1,2)$ ⇒ dependent; span is a line, dimension 1.
2. **Eigenvalues of $A=\begin{pmatrix}2&1\\1&2\end{pmatrix}$.** $|A-\lambda I|=(2-\lambda)^2-1=0\Rightarrow 2-\lambda=\pm1\Rightarrow\lambda=1,3$. Eigenvectors: $\lambda=3$: $(A-3I)v=\begin{pmatrix}-1&1\\1&-1\end{pmatrix}v=0\Rightarrow v=(1,1)^T$. $\lambda=1$: $v=(1,-1)^T$. Dot product $1-1=0$ — orthogonal, as the theorem says. Both eigenvalues $>0$ ⇒ PD; Sylvester agrees: $2>0$, $|A|=3>0$.
3. **Norms of $x=(3,-4)$:** $\|x\|_2=5$, $\|x\|_1=7$, $\|x\|_\infty=4$.

---

### 5. Summary

✅ Independence: $\sum c_iv_i=0\Rightarrow$ all $c_i=0$.
✅ Subspace = closed under scaling and adding; always contains $0$.
✅ Basis = smallest spanning independent set; dimension = its size.
✅ $\mathcal N(A)\perp$ row space; $\mathcal N(A^T)\perp\mathcal R(A)$.
✅ $Ax=b$ has 0, 1 or ∞ solutions; all solutions $=x_p+\mathcal N(A)$.
✅ rank + nullity $=n$.
✅ Eigen: $|A-\lambda I|=0$; symmetric ⇒ real eigenvalues, orthogonal eigenvectors.
✅ PSD ⟺ eigenvalues $\ge0$ ⟺ $v^TAv\ge0$; PD test by leading principal minors $>0$.
✅ Norms $\|\cdot\|_1,\|\cdot\|_2,\|\cdot\|_\infty$; Cauchy–Schwarz $|x^Ty|\le\|x\|\|y\|$.

<!-- LECTURE END: la -->

---

<!-- LECTURE START: calculus-aug7 -->
<a id="calculus-aug7"></a>

## Lecture 4: Calculus Review: limits, partial derivatives, gradient, contours, chain rule, Jacobian, directional derivative

*Source file:* `Calculus-aug7.pdf`

### 1. What this lecture is about

Optimization means finding where a function is smallest (or largest). To do that we need to know how a function **changes** when we move its input. Change is measured by derivatives. This lecture reviews:

1. Limits (what "approaches" means precisely, the $\varepsilon$–$N$ definition)
2. Partial derivatives and the **gradient** $\nabla f$
3. Contours (level sets), and why the gradient is perpendicular to them
4. Rules of differentiation: sum, product, **chain rule** (scalar and vector versions)
5. Vector-valued functions and the **Jacobian**
6. Two gradients you will use all the time: $\nabla(c^Tx)=c$ and $\nabla(x^TQx)=2Qx$
7. **Directional derivative** $\nabla f(x)^Tv$

---

### 2. Background you need first

- A **function** $f:\mathbb{R}^n\to\mathbb{R}$ takes a vector $x\in\mathbb{R}^n$ and returns one number. Example: $f(x)=x_1^2+x_2^2$ takes $(x_1,x_2)$ and returns a number.
- Ordinary derivative (1 variable): $f'(x)=\lim_{h\to0}\frac{f(x+h)-f(x)}{h}$, the slope of the tangent line.
- Basic 1-D rules: $(x^k)'=kx^{k-1}$, $(\sin x)'=\cos x$, $(e^x)'=e^x$, $(cf)'=cf'$, $(f+g)'=f'+g'$.
- The standard basis vector $e_i\in\mathbb{R}^n$ has a $1$ in position $i$ and $0$ elsewhere: $e_i=(0,\dots,0,1,0,\dots,0)^T$. So $x+\alpha e_i$ means "change only coordinate $i$ by $\alpha$".

---

### 3. Main content with full derivations

#### 3.1 Limits

**(1) $\lim_{x\to a}f(x)=L$**: as $x$ gets close to $a$, $f(x)$ gets close to $L$. The function does *not* have to be defined at $a$.

*Example:* $\displaystyle\lim_{x\to1}\frac{x^3-1}{x-1}=3$.
🔍 *Missing steps.* Plugging in $x=1$ gives $\frac{0}{0}$, which is undefined. Factor the top: $x^3-1=(x-1)(x^2+x+1)$ (check by multiplying out: $x^3+x^2+x-x^2-x-1=x^3-1$ ✓). For $x\neq1$ we can cancel:
$$\frac{x^3-1}{x-1}=x^2+x+1\ \xrightarrow{x\to1}\ 1+1+1=3.$$

**(2) $\lim_{x\to\infty}f(x)=L$** — precise definition:
$$\forall\varepsilon>0,\ \exists N>0\ \text{such that}\ x>N\Rightarrow|f(x)-L|<\varepsilon.$$
In words: however small a tolerance $\varepsilon$ you choose, eventually (beyond some $N$) $f(x)$ stays within $\varepsilon$ of $L$.

*Example:* $\lim_{x\to\infty}\frac1x=0$. Take $\varepsilon=0.3$. We need $N$ with $\left|\frac1x-0\right|<0.3$ for all $x>N$.
- For $x>0$: $\frac1x<0.3\iff x>\frac1{0.3}\approx3.33$.
- So $N=\frac1{0.3}$ works. 🔍 In general, for any $\varepsilon$ choose $N=\frac1\varepsilon$; this proves the limit.

**(3) $\lim_{x\to\infty}f(x)=\infty$** (the function grows without bound):
$$\forall N_1>0,\ \exists N_2>0\ \text{such that}\ x>N_2\Rightarrow f(x)>N_1.$$
*Example:* $\lim_{x\to\infty}x^2=\infty$. With $N_1=10000$: we need $x>N_2\Rightarrow x^2>10000$. Since $x^2>10000\iff x>100$ (for $x>0$), choose $N_2=100$. 🔍 In general $N_2=\sqrt{N_1}$.

#### 3.2 Partial derivatives

For $f:\mathbb{R}^n\to\mathbb{R}$, the **partial derivative with respect to $x_i$** is the ordinary derivative when only $x_i$ moves and all other coordinates are frozen:
$$\frac{\partial f}{\partial x_i}(x)=\lim_{\alpha\to0}\frac{f(x+\alpha e_i)-f(x)}{\alpha}.$$
(Picture in notes: from $x$, step a small amount $\alpha$ in the direction $e_i$ and see how $f$ changes.)

🔍 *How to compute in practice:* treat every other variable as a constant and differentiate as usual. E.g. $f=x_1^2x_2+3x_2$: $\frac{\partial f}{\partial x_1}=2x_1x_2$ (the $3x_2$ is a constant w.r.t. $x_1$, derivative 0), $\frac{\partial f}{\partial x_2}=x_1^2+3$.

#### 3.3 The gradient

Stack all partial derivatives into a column vector:
$$\nabla f(x)=\begin{pmatrix}\frac{\partial f}{\partial x_1}\\\frac{\partial f}{\partial x_2}\\\vdots\\\frac{\partial f}{\partial x_n}\end{pmatrix}\in\mathbb{R}^n .$$
✅ **The gradient points in the direction of maximum increase of $f$.** (Proved in §3.8 using Cauchy–Schwarz.) So $-\nabla f$ points in the direction of steepest *decrease* — the idea behind gradient descent.

#### 3.4 Contours (level sets)

The **contour** of $f$ at level $k$:
$$\{x: f(x)=k\}.$$
*Example:* $f(x)=x_1^2+x_2^2$. The contour $\{x:x_1^2+x_2^2=4\}$ is the circle of radius $2$. Different $k$ give concentric circles (radius $\sqrt k$).

The gradient: $\nabla f(x)=\begin{pmatrix}2x_1\\2x_2\end{pmatrix}=2x$. At a point $x_0$ on the circle, $\nabla f(x_0)=2x_0$ points **straight outward** from the centre. That is perpendicular to the circle's tangent line (the right angle drawn in the notes).

✅ **General fact: $\nabla f(x_0)$ is perpendicular (normal) to the contour through $x_0$.**
🔍 *Why.* Move along the contour on a curve $x(t)$ with $x(0)=x_0$. On the contour $f(x(t))=k$ is constant, so its derivative is $0$. By the chain rule (§3.6), $\frac{d}{dt}f(x(t))=\nabla f(x(t))^Tx'(t)=0$. At $t=0$: $\nabla f(x_0)^Tx'(0)=0$. Here $x'(0)$ is the tangent direction of the contour, and we just showed its dot product with the gradient is 0, i.e. they are perpendicular.

#### 3.5 Rules of differentiation

**(1) Sum rule.** $f,g:\mathbb{R}^n\to\mathbb{R}$: $\nabla(f+g)=\nabla f+\nabla g$ (do it coordinate by coordinate).

**(2) Product rule.** $h(x)=f(x)\,g(x)$.
- 1-D: $h'(x)=f'(x)g(x)+f(x)g'(x)$.
- $f,g:\mathbb{R}^n\to\mathbb{R}$: $\nabla h(x)=g(x)\nabla f(x)+f(x)\nabla g(x)$.

#### 3.6 The chain rule

**Scalar version.** $g:\mathbb{R}\to\mathbb{R}$, $f:\mathbb{R}\to\mathbb{R}$, $h(x)=(f\circ g)(x)=f(g(x))$:
$$h'(x)=\frac{df}{dg}\cdot\frac{dg}{dx}=f'(g(x))\,g'(x).$$
*Example:* $g(x)=\sin x$, $f(u)=u^2$, $h(x)=(\sin x)^2$.
🔍 $f'(u)=2u$, so $f'(g(x))=2\sin x$; $g'(x)=\cos x$. Hence $h'(x)=2\sin x\cos x$.

**Version $g:\mathbb{R}^n\to\mathbb{R}$, $f:\mathbb{R}\to\mathbb{R}$.** $h(x)=f(g(x))$:
$$\nabla h(x)=f'(g(x))\,\nabla g(x).$$
*Example:* $g(x)=x_1^2+\dots+x_n^2$, $f(u)=u^3$, so $h(x)=(x_1^2+\dots+x_n^2)^3$.
- $f'(u)=3u^2\Rightarrow f'(g(x))=3(x_1^2+\dots+x_n^2)^2$.
- $\nabla g(x)=(2x_1,2x_2,\dots,2x_n)^T$.
- $\nabla h(x)=3(x_1^2+\dots+x_n^2)^2\begin{pmatrix}2x_1\\\vdots\\2x_n\end{pmatrix}$.

**Neural-network example (homework in notes).** A single neuron takes inputs $a_1,\dots,a_d$, weights $w_1,\dots,w_d$, computes $w^Ta$ and outputs $z=\sigma(w^Ta)$ where $\sigma(u)=\frac1{1+e^{-u}}$ is the **sigmoid**. HW: find $\frac{\partial z}{\partial w}$.
🔍 *Solution.*
- Derivative of sigmoid: $\sigma(u)=(1+e^{-u})^{-1}$, so $\sigma'(u)=-(1+e^{-u})^{-2}\cdot(-e^{-u})=\frac{e^{-u}}{(1+e^{-u})^2}$. Rewrite: $\frac{1}{1+e^{-u}}\cdot\frac{e^{-u}}{1+e^{-u}}=\sigma(u)\big(1-\sigma(u)\big)$ (because $1-\sigma=\frac{e^{-u}}{1+e^{-u}}$).
- Inner function $g(w)=w^Ta=\sum_iw_ia_i$, $\nabla_wg=a$.
- Chain rule: $\frac{\partial z}{\partial w}=\sigma(w^Ta)\big(1-\sigma(w^Ta)\big)\,a$.

#### 3.7 Vector-valued functions and the Jacobian

$f:\mathbb{R}^n\to\mathbb{R}^m$ means $f(x)=\big(f_1(x),f_2(x),\dots,f_m(x)\big)^T$, i.e. $m$ scalar functions. Its derivative is the $m\times n$ **Jacobian** matrix:
$$Df(i,j)=\frac{\partial f_i}{\partial x_j},\qquad Df=\begin{pmatrix}-\ \nabla f_1^T\ -\\-\ \nabla f_2^T\ -\\\vdots\\-\ \nabla f_m^T\ -\end{pmatrix}.$$
Row $i$ is the (transposed) gradient of the $i$-th component.

**Chain rule for vector functions.** $g:\mathbb{R}^n\to\mathbb{R}^m$, $f:\mathbb{R}^m\to\mathbb{R}^t$, $h=f\circ g:\mathbb{R}^n\to\mathbb{R}^t$:
$$\boxed{Dh(x)=Df(g(x))\;Dg(x)}\qquad (t\times n)=(t\times m)(m\times n).$$
🔍 Note: for a scalar $f$ ($m=1$), $Df=\nabla f^T$ is a row vector. So the gradient is the transpose of the Jacobian.

#### 3.8 Two gradients you must know

**(a) Linear function** $f(x)=c^Tx=c_1x_1+c_2x_2+\dots+c_nx_n$.
$\frac{\partial f}{\partial x_i}=c_i$ (all other terms are constants), so
$$\nabla(c^Tx)=\begin{pmatrix}c_1\\\vdots\\c_n\end{pmatrix}=c.$$

**(b) Quadratic function** $f(x)=x^TQx+c^Tx+d$ with $Q$ **symmetric**:
$$\nabla f(x)=2Qx+c.$$
🔍 *Full derivation (skipped in notes).* Write out the quadratic form:
$$x^TQx=(x_1\cdots x_n)\begin{pmatrix}Q_{11}&\cdots&Q_{1n}\\\vdots&&\vdots\\Q_{n1}&\cdots&Q_{nn}\end{pmatrix}\begin{pmatrix}x_1\\\vdots\\x_n\end{pmatrix}=\sum_{i=1}^n\sum_{j=1}^nQ_{ij}x_ix_j .$$
Differentiate with respect to one variable $x_k$. A term $Q_{ij}x_ix_j$ depends on $x_k$ only if $i=k$ or $j=k$:
- terms with $i=k$: $\sum_jQ_{kj}x_kx_j$, derivative $\sum_jQ_{kj}x_j=(Qx)_k$;
- terms with $j=k$: $\sum_iQ_{ik}x_ix_k$, derivative $\sum_iQ_{ik}x_i=(Q^Tx)_k$;
- (the term $i=j=k$, $Q_{kk}x_k^2$, has derivative $2Q_{kk}x_k$ and is correctly counted once in each sum.)

So $\frac{\partial}{\partial x_k}x^TQx=(Qx)_k+(Q^Tx)_k$, i.e. $\nabla(x^TQx)=(Q+Q^T)x$. With $Q$ symmetric, $Q^T=Q$, so this is $2Qx$. Adding $\nabla(c^Tx)=c$ and $\nabla d=0$ gives $2Qx+c$. ∎

*Sanity check:* $n=1$: $f=qx^2+cx+d$, $f'=2qx+c$ ✓.

#### 3.9 Directional derivative

Question: **how fast does $f$ change at $x_0$ if we move in a direction $v$** (with $\|v\|=1$)?
$$D_vf(x_0)=\lim_{\alpha\to0}\frac{f(x_0+\alpha v)-f(x_0)}{\alpha}=\nabla f(x_0)^Tv .$$
🔍 *Why it equals $\nabla f^Tv$.* Define the 1-D function $\phi(\alpha)=f(x_0+\alpha v)$. The limit above is exactly $\phi'(0)$. Inner function $g(\alpha)=x_0+\alpha v$ has derivative $v$ (each coordinate $x_{0,i}+\alpha v_i$ has derivative $v_i$). By the chain rule, $\phi'(\alpha)=\nabla f(x_0+\alpha v)^Tv$, and at $\alpha=0$ this is $\nabla f(x_0)^Tv$.

🔍 *Why the gradient is the direction of steepest increase.* Among all unit vectors $v$, maximise $\nabla f(x_0)^Tv$. By Cauchy–Schwarz,
$$\nabla f(x_0)^Tv\le\|\nabla f(x_0)\|\,\|v\|=\|\nabla f(x_0)\|,$$
with equality exactly when $v$ points the same way as $\nabla f(x_0)$, i.e. $v=\nabla f(x_0)/\|\nabla f(x_0)\|$. Likewise the minimum, $-\|\nabla f(x_0)\|$, is reached at $v=-\nabla f/\|\nabla f\|$: **steepest descent**.

---

### 4. Worked examples

1. $f(x)=3x_1^2+2x_1x_2+x_2^2$. This is $x^TQx$ with $Q=\begin{pmatrix}3&1\\1&1\end{pmatrix}$ (the cross term $2x_1x_2$ is split as $Q_{12}=Q_{21}=1$). Then $\nabla f=2Qx=\begin{pmatrix}6x_1+2x_2\\2x_1+2x_2\end{pmatrix}$. Check directly: $\partial_1f=6x_1+2x_2$ ✓, $\partial_2f=2x_1+2x_2$ ✓.
2. Directional derivative of $f=x_1^2+x_2^2$ at $x_0=(1,1)$ along $v=(1,0)$: $\nabla f(x_0)=(2,2)$, so $D_vf=2$. Along $v=\frac1{\sqrt2}(1,-1)$ (tangent to the contour): $\frac{2-2}{\sqrt2}=0$. Moving along a contour does not change $f$ ✓.

---

### 5. Summary

✅ $\lim_{x\to\infty}f=L$: $\forall\varepsilon\,\exists N$: $x>N\Rightarrow|f(x)-L|<\varepsilon$.
✅ $\partial f/\partial x_i$ = derivative with only $x_i$ moving; $\nabla f$ = vector of these.
✅ $\nabla f$ is perpendicular to contours and points uphill fastest.
✅ Chain rule: $Dh=Df\,Dg$; scalar: $(f\circ g)'=f'(g)g'$.
✅ $\nabla(c^Tx)=c$; $\nabla(x^TQx)=2Qx$ for symmetric $Q$.
✅ Directional derivative $=\nabla f(x)^Tv$.

<!-- LECTURE END: calculus-aug7 -->

---

<!-- LECTURE START: aug10 -->
<a id="aug10"></a>

## Lecture 5: Calculus Review (Part 2): Hessian, continuity, Big-O / small-o, Taylor series

*Source file:* `Aug10.pdf`

> 📌 The first half of `Aug10.pdf` repeats the Aug 7 calculus notes: limits, gradient, chain rule, Jacobian, directional derivative. Those are explained in the block **"Calculus Review"**. This block covers only what is **new on Aug 10**.

### 1. What this lecture is about

1. The **Hessian** (matrix of second derivatives), the multivariable "second derivative"
2. Hessians of linear and quadratic functions
3. **Continuity** with the $\varepsilon$–$\delta$ definition
4. **Big-O** and **small-o** notation: comparing how fast functions grow or shrink
5. **Taylor series**: approximating a function by a polynomial, in 1-D and in $n$-D

These tools let us write $f(x_0+h)\approx f(x_0)+\nabla f(x_0)^Th+\tfrac12h^T\nabla^2f(x_0)h$. That one approximation is the basis of optimality conditions, gradient descent and Newton's method.

---

### 2. Background you need first

- Gradient $\nabla f(x)\in\mathbb{R}^n$, Jacobian, $\nabla(c^Tx)=c$, $\nabla(x^TQx)=2Qx$ (from the Calculus Review block).
- Factorial: $m!=1\cdot2\cdots m$, e.g. $3!=6$.
- $|\cdot|$ = absolute value; $\|\cdot\|$ = vector norm (length).

---

### 3. Main content with full derivations

#### 3.1 The Hessian

Let $f:\mathbb{R}^n\to\mathbb{R}$ be **twice differentiable**. The gradient $\nabla f(x)$ is itself a function $\mathbb{R}^n\to\mathbb{R}^n$. Its Jacobian is the **Hessian**:
$$\nabla^2f(x)=\begin{pmatrix}
\frac{\partial^2f}{\partial x_1^2}&\frac{\partial^2f}{\partial x_2\partial x_1}&\cdots&\frac{\partial^2f}{\partial x_n\partial x_1}\\
\frac{\partial^2f}{\partial x_1\partial x_2}&\frac{\partial^2f}{\partial x_2^2}&\cdots&\frac{\partial^2f}{\partial x_n\partial x_2}\\
\vdots&&\ddots&\vdots\\
\frac{\partial^2f}{\partial x_1\partial x_n}&\frac{\partial^2f}{\partial x_2\partial x_n}&\cdots&\frac{\partial^2f}{\partial x_n^2}
\end{pmatrix}\in\mathbb{R}^{n\times n}.$$
Entry $(i,j)$ = "differentiate w.r.t. $x_i$, then w.r.t. $x_j$". For twice *continuously* differentiable $f$ the order does not matter: $\frac{\partial^2f}{\partial x_i\partial x_j}=\frac{\partial^2f}{\partial x_j\partial x_i}$. That is why the notes mark the Hessian **symmetric**.

**Example 1: linear** $f(x)=c^Tx=\sum_ic_ix_i$.
$\nabla f(x)=c$, a constant vector. Differentiating a constant gives zero, so
$$\nabla^2f(x)=0_{n\times n}.$$

**Example 2: quadratic** $f(x)=x^TQx+c^Tx+d$, $Q$ symmetric.
$\nabla f(x)=2Qx+c$. 🔍 The $i$-th entry is $2\sum_jQ_{ij}x_j+c_i$; its derivative w.r.t. $x_j$ is $2Q_{ij}$. So
$$\nabla^2f(x)=2Q\quad(\text{constant, does not depend on }x).$$
**HW in notes: what if $Q$ is not symmetric?** 🔍 Answer: from the derivation in the Calculus Review, $\nabla(x^TQx)=(Q+Q^T)x$, so $\nabla f=(Q+Q^T)x+c$ and $\nabla^2f=Q+Q^T$. (Note $x^TQx=x^T\frac{Q+Q^T}{2}x$, so only the symmetric part of $Q$ matters anyway.)

#### 3.2 Continuity

$f$ is **continuous at $x_0$** if
$$\lim_{x\to x_0}f(x)=f(x_0).$$
Precisely ($\varepsilon$–$\delta$):
$$\forall\varepsilon>0,\ \exists\delta>0\ \text{such that}\ |x-x_0|<\delta\Rightarrow|f(x)-f(x_0)|<\varepsilon.$$
Meaning: to keep the output within $\varepsilon$ of $f(x_0)$, it is enough to keep the input within $\delta$ of $x_0$.

**Example from notes (discontinuous):**
$$f(x)=\begin{cases}\dfrac{1}{x-5}&x\neq5\\[4pt]2&x=5\end{cases}$$
(The notes write "$x-5$" with a graph of a hyperbola blowing up at $x=5$; the graph is that of $\frac1{x-5}$.) Take $\varepsilon=0.1$, so we want outputs inside $(2-0.1,\ 2+0.1)$. However small we take $\delta$ (e.g. $0.001$, giving inputs in $(4.999,5.001)$, or $\delta=10^{-4}$), points near 5 give $|f(x)|$ huge, nowhere near 2. So no $\delta$ works and $f$ is **not continuous at 5**.

#### 3.3 Big-O notation

**Asymptotic version ($x\to\infty$).** $f=O(g)$, read "$f$ is big-O of $g$", if
$$\exists C>0,\ N>0:\quad |f(x)|\le C|g(x)|\quad\forall x>N .$$
Meaning: eventually $f$ grows **no faster** than a constant times $g$. (Graph: after the point $N$, the curve of $Cg$ stays above $f$.)
*Example:* $n^2=O(n^3)$ as $n\to\infty$. 🔍 Take $C=1,N=1$: for $n>1$, $n^2\le n^3$ ✓.

**Local version (around $x\to a$).** $f=O(g)$ around $a$ if
$$\exists C>0,\ \delta>0:\quad |f(x)|\le C|g(x)|\quad\forall x:\ \|x-a\|<\delta .$$
*Example (around $n\to0$):* $n^3=O(n^2)$. Check: we need $C,\delta$ with $|n|^3\le C n^2$ for $|n|<\delta$. Take $C=1,\delta=1$: if $|n|<1$, then $|n|^3=|n|\cdot n^2\le n^2$ ✓.
⚠️ Near $0$, *higher* powers are *smaller* ($0.1^3=0.001<0.1^2=0.01$). Near $\infty$ it is the opposite.

#### 3.4 Small-o notation

**Asymptotic.** $f=o(g)$ ("$f$ is small-o of $g$", $f$ is *negligible* compared to $g$) if
$$\forall C>0\ \exists N>0:\ |f(x)|<C|g(x)|\ \forall x>N\quad\Longleftrightarrow\quad\lim_{x\to\infty}\frac{|f(x)|}{|g(x)|}=0 .$$
Difference from big-O: big-O needs the bound for **some** $C$; small-o needs it for **every** $C$, however tiny.

**Around $x\to a$:**
$$\lim_{x\to a}\frac{|f(x)|}{|g(x)|}=0,\quad\text{i.e. }\forall\varepsilon>0\ \exists\delta>0:\ \frac{|f(x)|}{|g(x)|}<\varepsilon\text{ whenever }\|x-a\|<\delta .$$

*Examples:* $n^2=o(n^3)$ around $n\to\infty$ (since $n^2/n^3=1/n\to0$); $n^3=o(n^2)$ around $n\to0$ (since $n^3/n^2=n\to0$).

**Exercises worked in the notes:**
1. Is $n^2+n^3=o(n^3)$ around $n\to0$?
   $\dfrac{|n^2+n^3|}{|n^3|}=\Big|\dfrac1n+1\Big|\to\infty\neq0$ as $n\to0$. ✗ **No.** (HW: the case $n\to\infty$. 🔍 There the ratio $\frac1n+1\to1\ne0$, so it is still **not** $o(n^3)$, although it *is* $O(n^3)$.)
2. Is $n^2+n^3=O(n^2)$ around $n\to0$? Need $C,\delta$ with $n^2+n^3\le Cn^2$ for $|n|<\delta$. Take $C=2$: need $n^2+n^3\le2n^2\iff n^3\le n^2$, true when $|n|<1$. ✓ **Yes** (with $\delta=1$; the notes take $|n|<\delta/2$, any small $\delta$ works).
3. But $n^2+n^3\neq o(n^2)$ around $0$: $\frac{n^2+n^3}{n^2}=1+n\to1\neq0$.

#### 3.5 Taylor series in one variable

For $f:\mathbb{R}\to\mathbb{R}$ smooth enough, near a point $x_0$:
$$f(x_0+h)=f(x_0)+hf'(x_0)+\frac{h^2}{2}f''(x_0)+\dots+\frac{h^m}{m!}f^{(m)}(x_0)+\dots$$
Truncations give approximations of increasing order:
- **0th order:** $f(x_0+h)\approx f(x_0)$.
- **1st order (linear / tangent line):** $f(x_0+h)\approx f(x_0)+hf'(x_0)$.
- **2nd order (quadratic):** $f(x_0+h)\approx f(x_0)+hf'(x_0)+\frac12h^2f''(x_0)$.

The error of the 2nd-order approximation is made of $h^3,h^4,\dots$ terms, which are $o(h^2)$. In general
$$\lim_{h\to0}\frac{h^{m+1}}{h^m}=\lim_{h\to0}h=0\quad\Rightarrow\quad h^{m+1}=o(h^m).$$
So we write exactly: $f(x_0+h)=f(x_0)+hf'(x_0)+\frac12h^2f''(x_0)+o(h^2)$.

🔍 *Where the coefficients come from.* Suppose $f(x_0+h)=a_0+a_1h+a_2h^2+a_3h^3+\dots$. Put $h=0$: $a_0=f(x_0)$. Differentiate once in $h$ and put $h=0$: $a_1=f'(x_0)$. Differentiate twice: $2a_2=f''(x_0)\Rightarrow a_2=\frac{f''(x_0)}{2}$. In general $m!\,a_m=f^{(m)}(x_0)$.

#### 3.6 Taylor series in $n$ variables

For $f:\mathbb{R}^n\to\mathbb{R}$, $x_0,h\in\mathbb{R}^n$:
$$\boxed{f(x_0+h)=f(x_0)+h^T\nabla f(x_0)+\tfrac12h^T\nabla^2f(x_0)\,h+\dots}$$
🔍 *Derivation from the 1-D version.* Let $\phi(t)=f(x_0+th)$, a function of one variable $t$. Then $f(x_0+h)=\phi(1)\approx\phi(0)+\phi'(0)+\frac12\phi''(0)$ (1-D Taylor with step 1).
- $\phi(0)=f(x_0)$.
- $\phi'(t)=\nabla f(x_0+th)^Th$ (directional derivative / chain rule), so $\phi'(0)=h^T\nabla f(x_0)$.
- $\phi'(t)=\sum_ih_i\frac{\partial f}{\partial x_i}(x_0+th)$; differentiate again: $\phi''(t)=\sum_i\sum_jh_ih_j\frac{\partial^2f}{\partial x_i\partial x_j}(x_0+th)=h^T\nabla^2f(x_0+th)h$, so $\phi''(0)=h^T\nabla^2f(x_0)h$.

Plug in to get the boxed formula. ∎

#### 3.7 Worked example from the notes: $f(x)=\ln(1+x)$ around $x_0=1$

- $f(1)=\ln 2$.
- $f'(x)=\dfrac{1}{1+x}\Rightarrow f'(1)=\dfrac12$.
- $f''(x)=-\dfrac{1}{(1+x)^2}\Rightarrow f''(1)=-\dfrac14$.

So
- **1st order:** $f(1+h)\approx\ln2+\dfrac h2$.
- **2nd order:** $f(1+h)\approx\ln2+\dfrac h2+\dfrac{h^2}{2}\cdot\Big(-\dfrac14\Big)=\ln2+\dfrac h2-\dfrac{h^2}{8}$.

🔍 *Numerical check* with $h=0.1$: true $\ln2.1\approx0.74194$. 1st order: $0.69315+0.05=0.74315$ (error $\approx1.2\times10^{-3}$). 2nd order: $0.74315-0.00125=0.74190$ (error $\approx4\times10^{-5}$). Each extra term makes the approximation much better for small $h$.

---

### 4. Worked example (multivariable Taylor)

$f(x)=x_1^2+x_1x_2+2x_2^2$ around $x_0=(1,1)$.
- $f(x_0)=1+1+2=4$.
- $\nabla f=(2x_1+x_2,\ x_1+4x_2)^T\Rightarrow\nabla f(x_0)=(3,5)^T$.
- $\nabla^2f=\begin{pmatrix}2&1\\1&4\end{pmatrix}$ (constant).
- $f(x_0+h)=4+3h_1+5h_2+\frac12(2h_1^2+2h_1h_2+4h_2^2)=4+3h_1+5h_2+h_1^2+h_1h_2+2h_2^2$. For a quadratic $f$ the 2nd-order Taylor expansion is **exact**: there are no higher terms.

---

### 5. Summary

✅ Hessian $\nabla^2f$ = matrix of second partials, symmetric. $\nabla^2(c^Tx)=0$, $\nabla^2(x^TQx)=2Q$.
✅ Continuity: $\lim_{x\to x_0}f(x)=f(x_0)$ ($\varepsilon$–$\delta$).
✅ $f=O(g)$: $|f|\le C|g|$ for some $C$. $f=o(g)$: $|f|/|g|\to0$.
✅ Near $0$: $h^{m+1}=o(h^m)$.
✅ Taylor: $f(x_0+h)=f(x_0)+h^T\nabla f(x_0)+\frac12h^T\nabla^2f(x_0)h+o(\|h\|^2)$.

<!-- LECTURE END: aug10 -->

---

<!-- LECTURE START: unconstrainedopt -->
<a id="unconstrainedopt"></a>

## Lecture 6: Unconstrained Optimization: problem setup, FONC & SONC (partial notes)

*Source file:* `UnconstrainedOpt.pdf`

> 📌 **This file is the first half of the Aug 17 notes (`Aug17.pdf`).** It stops in the middle of SONC. The full explanation, with complete proofs of FONC, SONC and SOSC and the worked HW example, is in the block **"Optimization Problems & Optimality Conditions (FONC, SONC, SOSC)"**. This block is only a quick checklist.

### Checklist of what this part covers

1. **Problem form:** $\min f(x)$ s.t. $x\in\Omega$. $f$ = objective; $\Omega$ = feasible set; $x^*=\arg\min$ = optimal solution.
2. $\max f=-\min(-f)$.
3. **Existence cases:** unique optimum, multiple optima, no optimum ($\min e^{-x}$), unbounded ($\min -x\to-\infty$).
4. **Global min:** $f(x^*)\le f(x)\ \forall x\in\Omega$. **Local min:** $\exists R>0$ with $f(x^*)\le f(x)$ whenever $\|x-x^*\|<R$.
5. **Examples:** least squares $\min_\beta\frac1n\sum(y^{(i)}-\beta^Tx^{(i)})^2$; logistic regression $\max_\beta$ log-likelihood with $\sigma(z)=1/(1+e^{-z})$.
6. **FONC:** local min ⇒ $\nabla f(x^*)=0$. Proof idea: $g(\alpha)=f(x^*+\alpha h)$ gives $\nabla f(x^*)^Th\ge0$ for $h$ and for $-h$, so it is $=0$. Then take $h=e_i$.
7. Saddle point example: $f(x)=x^3$ at $0$ (FONC holds, not a min).
8. **SONC:** local min ⇒ $\nabla f(x^*)=0$ and $\nabla^2f(x^*)\succeq0$, via Taylor expansion and the $o(\|v\|^2)$ argument.

✅ See the full block for every derivation step.

<!-- LECTURE END: unconstrainedopt -->

---

<!-- LECTURE START: topology -->
<a id="topology"></a>

## Lecture 7: Basic Topology: open/closed/bounded/compact sets, sup & inf, Weierstrass theorem

*Source file:* `Topology..pdf`

### 1. What this lecture is about

Before asking "where is the minimum?", we must ask "**does a minimum even exist?**". For example, $\min x$ over $0<x<1$ has no answer: you can get closer and closer to 0 but never reach it. Topology gives the language to say when a minimum is guaranteed:

1. Convergent sequences
2. Balls (neighbourhoods)
3. Open sets, limit points, closed sets, interior
4. Bounded functions and bounded sets; compact sets
5. Least upper bound / greatest lower bound, supremum vs maximum, infimum vs minimum
6. **Weierstrass theorem:** a continuous function on a closed and bounded set always attains its min and max.

---

### 2. Background you need first

- $\|x\|_2=\sqrt{x_1^2+\dots+x_n^2}$ is the distance from $x$ to the origin, so $\|x-x_0\|_2$ is the distance between $x$ and $x_0$.
- Interval notation: $(a,b)=\{a<x<b\}$ (open, endpoints excluded), $[a,b]=\{a\le x\le b\}$ (closed, endpoints included), $[a,b)$ = half-open.
- $\Omega$ (capital omega) will denote a set in $\mathbb{R}^n$.

---

### 3. Main content with full explanations

#### 3.1 Convergent sequences

A sequence $x_1,x_2,\dots$ **converges to $x$** if $\lim_{n\to\infty}x_n=x$, i.e. 🔍 for every $\varepsilon>0$ there is $N$ such that $\|x_n-x\|<\varepsilon$ for all $n>N$ (eventually all terms stay within any tolerance of $x$). Example: $x_n=\frac1n\to0$.

#### 3.2 Ball / sphere / neighbourhood

The (open) **ball** around $x_0\in\mathbb{R}^n$ of radius $r>0$:
$$B(x_0,r)=\{x:\ \|x-x_0\|_2<r\}.$$
In 2-D it is a disc without its boundary circle. It is called a **neighbourhood** of $x_0$.

#### 3.3 Open sets

$\Omega\subseteq\mathbb{R}^n$ is **open** $\iff$ around **every** point $x_0\in\Omega$ there is a ball of *some* radius $r>0$ completely contained in $\Omega$.

Intuition: an open set does not contain any of its "edge" points. From any point you can wiggle a little in every direction and stay inside.

**Examples:**
- $\mathbb{R}^n$ itself is open (any ball fits).
- $\Omega=\{x\in\mathbb{R}:0<x<1\}=(0,1)$ is **open**. 🔍 For a point $x_0\in(0,1)$, take $r=\min(x_0,1-x_0)$ (distance to the nearer end); then $(x_0-r,x_0+r)\subseteq(0,1)$.
- $\{0\le x<1\}$ is **not open** ✗: around $x_0=0$ every ball $(-r,r)$ contains negative numbers, which are not in the set.
- $\Omega=\{x\in\mathbb{R}^2:\|x\|_1<1\}$ (the inside of the diamond $|x_1|+|x_2|<1$) is **open**. Every interior point has a small disc around it inside the diamond (the little circles in the drawing).

#### 3.4 Limit points

$x$ is a **limit point** of $\Omega$ if **every** $\varepsilon$-neighbourhood of $x$ contains a point $x_0\in\Omega$ with $x_0\neq x$.
Equivalent: there is a sequence $\{x_n\}$ in $\Omega$ (with $x_n\neq x$) that converges to $x$.

Note: a limit point need **not** itself belong to $\Omega$.

#### 3.5 Closed sets

$\Omega$ is **closed** $\iff$ **all limit points of $\Omega$ belong to $\Omega$**. (You cannot "escape" a closed set by taking limits.)

**Examples:**
1. $\Omega=[0,1]=\{0\le x\le1\}$: **closed**. Every $x\in\Omega$ is a limit point, and nothing outside is (e.g. $-0.5$ has a ball around it that misses $[0,1]$ entirely).
2. $\Omega=(0,1)$: **not closed**. Its limit points are $[0,1]$. 🔍 $0$ is a limit point (the sequence $\frac1n\in\Omega$ converges to 0) but $0\notin\Omega$.
3. $\Omega=\{x\in\mathbb{R}^2:\|x\|_2<1,\ x\neq0\}$ (punctured open disc): **not closed**. Limit points $=\{x:\|x\|_2\le1\}$, which include the boundary circle and the origin, neither in $\Omega$.
4. $\Omega=\{x\in\mathbb{R}:x\text{ is an integer}\}$: **closed**. It has **no** limit points at all (each integer is isolated: a ball of radius $0.5$ contains no other integer). So "all limit points belong to $\Omega$" is true trivially.

⚠️ "Not open" does **not** mean "closed". $[0,1)$ is neither open nor closed. $\mathbb{R}^n$ and $\emptyset$ are both open and closed.

🔍 *Useful fact:* $\Omega$ is closed iff its complement $\mathbb{R}^n\setminus\Omega$ is open.

#### 3.6 Interior

The **interior** of $\Omega$ is the **largest open set contained in $\Omega$** (all points that have a whole ball inside $\Omega$). E.g. interior of $[0,1]$ is $(0,1)$.

#### 3.7 Bounded functions and bounded sets

A function $f:\mathbb{R}^n\to\mathbb{R}$ is
- **bounded above** if $\exists M:\ f(x)<M\ \forall x$ ($M$ = an upper bound);
- **bounded below** if $\exists M:\ f(x)>M\ \forall x$ ($M$ = a lower bound).

A set $\Omega\subseteq\mathbb{R}^n$ is **bounded** if $\exists M:\ \|x\|<M\ \forall x\in\Omega$ (the set fits inside some big ball).

Examples: $\mathbb{R}^n$ — not bounded ✗. $\{x\in\mathbb{R}^2:x_1>0,x_2>0\}$ (first quadrant) — not bounded ✗. $\{x\in\mathbb{R}^2:\|x\|<1\}$ — bounded ✓.

#### 3.8 Compact sets

In $\mathbb{R}^n$: **compact = closed and bounded.**
E.g. $[0,1]$, the closed disc $\{\|x\|\le1\}$. Not compact: $(0,1)$ (not closed), $[0,\infty)$ (not bounded).

#### 3.9 Least upper bound and greatest lower bound

- $M$ is the **least upper bound (l.u.b.)** of a set if it is an upper bound and for any upper bound $M_1$, $M_1\ge M$ (it is the smallest upper bound).
- $M$ is the **greatest lower bound (g.l.b.)** if it is a lower bound and for any lower bound $M_1$, $M_1\le M$.

**Supremum** $=$ l.u.b. **Infimum** $=$ g.l.b.

#### 3.10 Supremum vs maximum, infimum vs minimum

$$\sup\ x\quad\text{s.t.}\quad x<5 .$$
The feasible set is $(-\infty,5)$. Upper bounds: every number $\ge5$. The least one is $5$, so **$\sup=5$**. But $5$ is not in the set, so there is **no maximum**. 🔍 Any candidate $x<5$ is beaten by $\frac{x+5}{2}$, which is also $<5$ and larger.

Similarly $\sup x$ s.t. $0<x<5$ is $5$ (no max), and $\inf x$ over it is $0$ (no min).

✅ The maximum exists exactly when the supremum is **achieved** (belongs to the set). The same holds for minimum and infimum.

#### 3.11 Weierstrass theorem (existence of optimal solutions)

> Let $f:\Omega\to\mathbb{R}$ with $\Omega\subseteq\mathbb{R}^n$. If **$f$ is continuous** and **$\Omega$ is closed and bounded** (compact), then the infimum and supremum of $f$ over $\Omega$ are **achieved**. That is, the problems
> $$\min_{x\in\Omega}f(x)\qquad\text{and}\qquad\max_{x\in\Omega}f(x)$$
> both have solutions.

🔍 *Why each condition is needed (counter-examples):*
- **Not closed:** $\min x$ over $(0,1)$: $\inf=0$, never reached.
- **Not bounded:** $\min e^{-x}$ over $[0,\infty)$: $\inf=0$, never reached ($e^{-x}>0$ always).
- **Not continuous:** $f(x)=x$ for $x\in(0,1]$, $f(0)=1$ on $[0,1]$: values approach $0$ but never equal 0.

⚠️ The theorem is **sufficient, not necessary**. A minimum can exist without these conditions (e.g. $\min x^2$ over $\mathbb{R}$ has min $0$ at $x=0$, though $\mathbb{R}$ is unbounded).

---

### 4. Worked examples

1. Is $\{x\in\mathbb{R}^2: x_1^2+x_2^2\le4\}$ compact? Closed (contains its boundary circle) and bounded (inside a ball of radius 3) ⇒ **compact**. So, e.g., $\min (x_1+x_2)$ over it has a solution (by Weierstrass).
2. $\{x\in\mathbb{R}: x\ge0\}$: closed but not bounded ⇒ not compact.
3. $\sup\{1-\frac1n:n=1,2,\dots\}=1$, not attained ⇒ no max.

---

### 5. Summary

✅ Open: every point has a ball inside. Closed: contains all its limit points.
✅ Bounded: fits in a ball. Compact (in $\mathbb{R}^n$) = closed + bounded.
✅ sup/inf always exist (possibly $\pm\infty$); max/min exist only if attained.
✅ **Weierstrass:** continuous $f$ + compact $\Omega$ ⇒ min and max exist.

<!-- LECTURE END: topology -->

---

<!-- LECTURE START: aug17 -->
<a id="aug17"></a>

## Lecture 8: Optimization Problems & Optimality Conditions (FONC, SONC, SOSC)

*Source file:* `Aug17.pdf`

### 1. What this lecture is about

1. The standard form of an optimization problem: objective, feasible set, optimal solution
2. Turning max into min
3. When an optimum exists: unique, multiple, none, unbounded
4. Global vs local minima
5. Two machine-learning examples: least squares (linear regression) and logistic regression
6. **Optimality conditions** for unconstrained problems, with full proofs:
   - **FONC** (first-order necessary): $\nabla f(x^*)=0$
   - **SONC** (second-order necessary): $\nabla^2f(x^*)\succeq0$
   - **SOSC** (second-order sufficient): $\nabla f(x^*)=0$ and $\nabla^2f(x^*)\succ0$ ⇒ local min
7. A key fact used in the proof: $\min_{\|v\|=1}v^TAv=\lambda_{\min}(A)$

---

### 2. Background you need first

- Gradient, Hessian, Taylor expansion $f(x^*+h)=f(x^*)+\nabla f(x^*)^Th+\frac12h^T\nabla^2f(x^*)h+o(\|h\|^2)$ (Calculus blocks).
- PSD/PD matrices and eigenvalues (Linear Algebra block): $A\succeq0\iff v^TAv\ge0\ \forall v$.
- Small-o: $g(v)=o(\|v\|^2)$ means $\frac{g(v)}{\|v\|^2}\to0$ as $v\to0$.

---

### 3. Main content with full derivations

#### 3.1 The general problem

$$\min_x\ f(x)\quad\text{s.t.}\quad x\in\Omega$$
- $f$: **objective function** (what we want small: cost, error, …).
- $\Omega$: **constraint set / feasible set / feasible region**.
- "s.t." = "subject to".
- Any $x\in\Omega$ is a **feasible point**.
- $x^*=\arg\min_{x\in\Omega}f(x)$ is an **optimal solution**: the feasible point with the smallest $f$. ("arg min" = the *point* where the min happens; "min" = the *value* $f(x^*)$.)

**Max ⇄ min.**
$$\max_{x\in\Omega}f(x)\;=\;-\min_{x\in\Omega}\big(-f(x)\big).$$
🔍 *Why:* flipping the graph upside down (multiplying by $-1$) turns the highest peak of $f$ into the lowest valley of $-f$ at the same $x$. Then negate the value back. So it is enough to study minimisation.

#### 3.2 Does an optimum exist?

- The feasible set must be **non-empty** (otherwise there is nothing to choose).
- **Unique optimum:** e.g. a single valley.
- **Multiple optima:** several valleys of the same depth.
- **No optimum though bounded below:** $\min e^{-x}$ over $\mathbb{R}$. The values approach $0$ but never reach it (inf $=0$, no min).
- **Unbounded:** $\min -x$ over $\mathbb{R}$. You can make $-x\to-\infty$; there is a sequence $x_n$ (e.g. $x_n=n$) with $f(x_n)\to-\infty$.

#### 3.3 Unconstrained optimization; global vs local minima

Unconstrained: $\Omega=\mathbb{R}^n$:
$$\min f(x)\quad\text{s.t.}\quad x\in\mathbb{R}^n .$$

- **Global minimum** $x^*$: $f(x^*)\le f(x)$ for **all** $x\in\Omega$.
- **Local minimum** $x^*$: $\exists R>0$ such that $\|x-x^*\|<R\Rightarrow f(x^*)\le f(x)$. Best only in a small neighbourhood.

(The wavy curve in the notes has several local minima, one global min and a global max.) Every global min is a local min, not conversely. With constraints, minima can also sit on the **boundary or corners** of $\Omega$ (the pentagon picture).

#### 3.4 Examples from machine learning

**(1) Least squares / linear regression.** Data $\mathcal D=\{(x^{(i)},y^{(i)})\}_{i=1}^n$, $x^{(i)}\in\mathbb{R}^d$ (features), $y^{(i)}\in\mathbb{R}$ (label). Model: $y=\beta^Tx+\varepsilon$ (linear + noise). Choose $\beta$ to minimise the mean squared error:
$$\min_{\beta\in\mathbb{R}^d}\ \frac1n\sum_{i=1}^n\big(y^{(i)}-\beta^Tx^{(i)}\big)^2 .$$
(Solved in the Linear Regression block.)

**(2) Logistic regression** (labels $y_i\in\{0,1\}$). Model: $y_i=1$ with probability $\sigma(\beta^Tx^i)$ and $0$ with probability $1-\sigma(\beta^Tx^i)$, with $\sigma(z)=\frac{1}{1+e^{-z}}$. Choose $\beta$ by
$$\max_\beta\ \text{log-likelihood}(\beta).$$
🔍 *Written out:* the probability of observing $y_i$ is $\sigma_i^{y_i}(1-\sigma_i)^{1-y_i}$ with $\sigma_i=\sigma(\beta^Tx^i)$; assuming independence, multiply over $i$ and take $\log$:
$$\ell(\beta)=\sum_{i=1}^n\Big[y_i\log\sigma(\beta^Tx^i)+(1-y_i)\log\big(1-\sigma(\beta^Tx^i)\big)\Big].$$

---

#### 3.5 FONC: first-order necessary condition

> **Theorem (FONC).** Let $f:\mathbb{R}^n\to\mathbb{R}$ be differentiable and $x^*$ a local minimum. Then
> $$\nabla f(x^*)=0 .$$

**Proof, every step.**
1. Since $x^*$ is a local min, $\exists R>0$: $f(x^*)\le f(x)$ for all $x$ with $\|x-x^*\|<R$.
2. Take **any** direction $h\in\mathbb{R}^n$ and look along the line $x^*+\alpha h$, $\alpha\in\mathbb{R}$. Define the 1-D function
$$g(\alpha)=f(x^*+\alpha h),\qquad g(0)=f(x^*).$$
3. By the chain rule (directional derivative): $g'(\alpha)=\nabla f(x^*+\alpha h)^Th$, so $g'(0)=\nabla f(x^*)^Th$.
4. By definition of the derivative,
$$g'(0)=\lim_{v\to0^+}\frac{g(v)-g(0)}{v}=\lim_{v\to0^+}\frac{f(x^*+vh)-f(x^*)}{v}.$$
For small enough $v>0$ the point $x^*+vh$ is inside the ball of radius $R$, so the numerator is $\ge0$ (step 1). The denominator $v>0$. So the fraction is $\ge0$, hence its limit is $\ge0$:
$$\nabla f(x^*)^Th\ge0 .$$
5. This holds for *every* $h$; in particular for $-h$: $\nabla f(x^*)^T(-h)\ge0$, i.e. $\nabla f(x^*)^Th\le0$.
6. Both $\ge0$ and $\le0$ ⇒ $\nabla f(x^*)^Th=0$ for every $h$.
7. Choose $h=e_1=(1,0,\dots,0)^T$: $\nabla f(x^*)^Te_1=\frac{\partial f}{\partial x_1}(x^*)=0$. Similarly $h=e_i$ gives $\frac{\partial f}{\partial x_i}(x^*)=0$ for each $i$. So every component of $\nabla f(x^*)$ is zero. ∎

A point with $\nabla f(x^*)=0$ is called a **critical (stationary) point**.

⚠️ **Necessary, not sufficient.** $f(x)=x^3$: $f'(0)=0$ but $0$ is not a min (to the left $f<0$). Such points are **saddle points / inflection points**. Maxima also satisfy $\nabla f=0$.

#### 3.6 SONC: second-order necessary condition

> **Theorem (SONC).** Let $f$ be twice differentiable and $x^*$ a local minimum. Then
> 1. $\nabla f(x^*)=0$, and
> 2. $\nabla^2f(x^*)\succeq0$ (positive semidefinite: $v^T\nabla^2f(x^*)v\ge0\ \forall v$).

**Proof.**
1. Part 1 is FONC.
2. Taylor expansion around $x^*$ with step $v$ (write $v=\alpha h$):
$$f(x^*+v)=f(x^*)+\underbrace{\nabla f(x^*)^Tv}_{=0\text{ by FONC}}+\tfrac12v^T\nabla^2f(x^*)v+g(v),\qquad g(v)=o(\|v\|^2).$$
3. For $\|v\|$ small enough, local minimality gives $f(x^*+v)-f(x^*)\ge0$, so
$$\tfrac12v^T\nabla^2f(x^*)v+g(v)\ge0 .$$
4. Divide by $\|v\|^2>0$:
$$\frac{\frac12v^T\nabla^2f(x^*)v}{\|v\|^2}\ \ge\ -\frac{g(v)}{\|v\|^2}.$$
5. 🔍 The left side depends only on the *direction* $u=v/\|v\|$: it equals $\frac12u^T\nabla^2f(x^*)u$ (scaling $v$ by $t$ multiplies top and bottom by $t^2$). So fix a unit direction $u$ and let $v=tu$, $t\to0^+$.
6. Since $g(v)=o(\|v\|^2)$: $\forall\varepsilon>0\ \exists\delta>0$ with $\frac{|g(v)|}{\|v\|^2}<\varepsilon$ whenever $\|v\|<\delta$. So the right side is $>-\varepsilon$.
7. Hence $\frac12u^T\nabla^2f(x^*)u>-\varepsilon$ for **every** $\varepsilon>0$. A fixed number bigger than every negative number must be $\ge0$. So $u^T\nabla^2f(x^*)u\ge0$ for every unit $u$, hence (scale) for every $v\in\mathbb{R}^n$. ∎

**Examples.**
1. $f(x)=x_1^2+0.5x_2^2+3x_2$ on $\mathbb{R}^2$.
   - $\nabla f=\begin{pmatrix}2x_1\\x_2+3\end{pmatrix}=0\Rightarrow x_1=0,\ x_2=-3$. Critical point $(0,-3)$.
   - $\nabla^2f=\begin{pmatrix}2&0\\0&1\end{pmatrix}$. Eigenvalues (diagonal entries) $2,1>0$ ⇒ PD ⇒ SONC holds (and, by SOSC below, it is a local min).
2. $f(x)=x_1^2-0.5x_2^2+3x_2$.
   - $\nabla f=\begin{pmatrix}2x_1\\-x_2+3\end{pmatrix}=0\Rightarrow x_1=0,\ x_2=3$. Critical point $(0,3)$.
   - $\nabla^2f=\begin{pmatrix}2&0\\0&-1\end{pmatrix}$: eigenvalues $2$ and $-1$. **Not PSD** ($\not\succeq0$), so by SONC it is **not a local min**. Also $-\nabla^2f$ is not PSD (eigenvalue $-2$), so not a local max either: it is a **saddle point** (up along $x_1$, down along $x_2$).

#### 3.7 A key lemma: $\min_{\|v\|=1}v^TAv=\lambda_1$

Let $A$ be symmetric with eigenvalues $\lambda_1\le\lambda_2\le\dots\le\lambda_n$ (roots of $|A-\lambda I|=0$). Then
$$\boxed{\min_{\|v\|=1}v^TAv=\lambda_1\ (\text{smallest eigenvalue})}\quad\Rightarrow\quad v^TAv\ge\lambda_1\|v\|^2\ \ \forall v .$$

🔍 *Proof (sketched in notes, completed here).*
1. For a symmetric $A$, **the eigenvectors span $\mathbb{R}^n$** and can be chosen orthonormal $v_1,\dots,v_n$ ($v_i^Tv_j=0$ for $i\ne j$, $\|v_i\|=1$). For distinct eigenvalues orthogonality was proved in the LA block. If an eigenvalue repeats, its eigenspace $\mathcal N(A-\lambda_iI)$ has dimension $>1$ (e.g. $w_1,w_2$ with $Aw_1=\lambda_iw_1$, $Aw_2=\lambda_iw_2$), and we pick an orthonormal basis inside it.
2. Any $z\in\mathbb{R}^n$ can be written $z=\alpha_1v_1+\dots+\alpha_nv_n$. If $\|z\|=1$, then $\sum\alpha_i^2=1$ (orthonormality: $z^Tz=\sum_{i,j}\alpha_i\alpha_jv_i^Tv_j=\sum_i\alpha_i^2$).
3. $Az=\sum\alpha_i\lambda_iv_i$, so $z^TAz=\sum_{i,j}\alpha_j\alpha_i\lambda_iv_j^Tv_i=\sum_i\lambda_i\alpha_i^2$.
4. Since $\lambda_i\ge\lambda_1$: $z^TAz\ge\lambda_1\sum_i\alpha_i^2=\lambda_1$. Equality at $z=v_1$. ∎

#### 3.8 SOSC: second-order sufficient condition

> **Theorem (SOSC).** Let $f$ be twice differentiable. If $x^*$ satisfies
> $$\nabla f(x^*)=0\quad\text{and}\quad\nabla^2f(x^*)\succ0\ (\text{positive definite}),$$
> then $x^*$ is a (strict) **local minimizer**.

**Proof.**
1. Take any $y\in\mathbb{R}^n$ (small step). Taylor:
$$f(x^*+y)=f(x^*)+\underbrace{\nabla f(x^*)^Ty}_{=0}+\tfrac12y^T\nabla^2f(x^*)y+g(y),\qquad g(y)=o(\|y\|^2).$$
So $f(x^*+y)-f(x^*)=\frac12y^T\nabla^2f(x^*)y+g(y)$. We must show this is $\ge0$ for small $y$.
2. Let $\lambda_1>0$ be the smallest eigenvalue of $\nabla^2f(x^*)$ (positive because PD). By the lemma, $y^T\nabla^2f(x^*)y\ge\lambda_1\|y\|^2$. Hence
$$f(x^*+y)-f(x^*)\ \ge\ \frac{\lambda_1}{2}\|y\|^2+g(y)=\|y\|^2\Big(\frac{\lambda_1}{2}+\frac{g(y)}{\|y\|^2}\Big).$$
3. Because $\lim_{\|y\|\to0}\frac{g(y)}{\|y\|^2}=0$: for every $\varepsilon>0$ there is $\delta>0$ with $\frac{|g(y)|}{\|y\|^2}<\varepsilon$ whenever $\|y\|<\delta$ (equivalently $|g(y)|<\varepsilon\|y\|^2$).
4. **Choose $\varepsilon=\frac{\lambda_1}{2}$.** Then for all $y$ with $0<\|y\|<\delta$: $\frac{g(y)}{\|y\|^2}>-\frac{\lambda_1}{2}$, so the bracket $\frac{\lambda_1}{2}+\frac{g(y)}{\|y\|^2}>0$.
5. Therefore $f(x^*+y)-f(x^*)>0$ for all $0<\|y\|<\delta$: $x^*$ is a strict local minimizer. ∎

⚠️ SOSC is **sufficient but not necessary**: $f(x)=x^4$ has a min at 0 but $f''(0)=0$ (not PD). There SOSC is silent.

#### 3.9 Summary of the logic

| Condition | Statement | Type |
|---|---|---|
| FONC | local min ⇒ $\nabla f=0$ | necessary |
| SONC | local min ⇒ $\nabla f=0$, $\nabla^2f\succeq0$ | necessary |
| SOSC | $\nabla f=0$, $\nabla^2f\succ0$ ⇒ strict local min | sufficient |

Recipe: (1) solve $\nabla f=0$ to get candidates; (2) compute $\nabla^2f$ at each; PD ⇒ local min, ND ⇒ local max, indefinite ⇒ saddle, semidefinite ⇒ inconclusive (look closer).

---

### 4. Worked example: HW from the notes

$f(x)=2x_1^2+x_2^2-2x_1x_2+2x_1^3+x_1^4$.

**Step 1: gradient.**
$$\frac{\partial f}{\partial x_1}=4x_1-2x_2+6x_1^2+4x_1^3,\qquad\frac{\partial f}{\partial x_2}=2x_2-2x_1 .$$
**Step 2: critical points.** From the second equation, $x_2=x_1$. Substitute into the first:
$$4x_1-2x_1+6x_1^2+4x_1^3=2x_1+6x_1^2+4x_1^3=2x_1(2x_1^2+3x_1+1)=2x_1(2x_1+1)(x_1+1)=0.$$
So $x_1\in\{0,-\tfrac12,-1\}$, giving critical points $(0,0)$, $(-\tfrac12,-\tfrac12)$, $(-1,-1)$.

**Step 3: Hessian.**
$$\nabla^2f=\begin{pmatrix}4+12x_1+12x_1^2&-2\\-2&2\end{pmatrix}.$$
- $(0,0)$: $\begin{pmatrix}4&-2\\-2&2\end{pmatrix}$. Leading minors $4>0$, $8-4=4>0$ ⇒ PD ⇒ **strict local min**, $f=0$.
- $(-\frac12,-\frac12)$: top-left $4-6+3=1$; $\begin{pmatrix}1&-2\\-2&2\end{pmatrix}$, det $2-4=-2<0$ ⇒ indefinite ⇒ **saddle point**.
- $(-1,-1)$: top-left $4-12+12=4$; same matrix as at the origin ⇒ PD ⇒ **strict local min**. $f(-1,-1)=2+1-2-2+1=0$.

Both local minima have value $0$. 🔍 In fact $f=(x_1-x_2)^2+x_1^2(x_1+1)^2\ge0$ (expand to check), so both are **global** minima.

---

### 5. Summary

✅ $\max f=-\min(-f)$. An optimum may be unique, multiple, non-existent or unbounded.
✅ Local min: best in a neighbourhood; global: best everywhere.
✅ FONC $\nabla f(x^*)=0$, proved by looking along every line through $x^*$.
✅ SONC adds $\nabla^2f(x^*)\succeq0$; SOSC: $\nabla f=0$ & $\nabla^2f\succ0$ ⇒ local min.
✅ $v^TAv\ge\lambda_{\min}\|v\|^2$ for symmetric $A$.

<!-- LECTURE END: aug17 -->

---

<!-- LECTURE START: linearreg -->
<a id="linearreg"></a>

## Lecture 9: Linear Regression via Optimality Conditions (Normal Equations)

*Source file:* `LinearReg.pdf`

### 1. What this lecture is about

This lecture applies the optimality conditions (FONC/SONC/SOSC) to a real problem: **fitting a line/hyperplane to data by least squares**. We:

1. Recap FONC, SONC, SOSC
2. Write linear regression in matrix form $\min_\beta\|X\beta-z\|^2$
3. Expand, differentiate, and solve $\nabla f(\beta)=0$, giving the **normal equations** $\beta^*=(X^TX)^{-1}X^Tz$
4. Prove that $X^TX$ is invertible when the columns of $X$ are independent
5. Show the Hessian $2X^TX$ is PD, so $\beta^*$ is a (in fact global) minimum

---

### 2. Background you need first

- $\nabla(c^Tx)=c$, $\nabla(x^TQx)=2Qx$ for symmetric $Q$; Hessian of $x^TQx$ is $2Q$ (Calculus blocks).
- $\|v\|^2=v^Tv$; $(AB)^T=B^TA^T$; a number equals its own transpose.
- Null space $\mathcal N(A)=\{q:Aq=0\}$; a square matrix is invertible ⟺ $\det\neq0$ ⟺ full rank ⟺ $\mathcal N=\{0\}$ (LA block).

---

### 3. Main content with full derivations

#### 3.1 Recap of optimality conditions

For $\min f(x)$, $x\in\mathbb{R}^n$, $f:\mathbb{R}^n\to\mathbb{R}$:
1. $f$ differentiable, $x^*$ local min ⇒ $\nabla f(x^*)=0$ (FONC).
2. $f$ twice differentiable, $x^*$ local min ⇒ $\nabla f(x^*)=0$ and $\nabla^2f(x^*)\succeq0$ (SONC).
3. $f$ twice differentiable, $\nabla f(x^*)=0$ and $\nabla^2f(x^*)\succ0$ ⇒ $x^*$ local min (SOSC).

$\{x:\nabla f(x)=0\}$ = **critical points**.

#### 3.2 The regression problem

Data: pairs $(x,y)$ with $x\in\mathbb{R}^d$ (features), $y\in\mathbb{R}$ (target); dataset $\mathcal D=\{(x^{(i)},y^{(i)})\}_{i=1}^n$. Model $y\approx\beta^Tx=\beta_1x_1+\dots+\beta_dx_d$. Least squares:
$$\min_\beta\sum_{i=1}^n\big(y^{(i)}-\beta^Tx^{(i)}\big)^2 .$$

**Matrix form.** Stack the data:
$$X=\begin{pmatrix}-\ x^{(1)T}\ -\\-\ x^{(2)T}\ -\\\vdots\\-\ x^{(n)T}\ -\end{pmatrix}\in\mathbb{R}^{n\times d},\qquad\beta=\begin{pmatrix}\beta_1\\\vdots\\\beta_d\end{pmatrix},\qquad z=\begin{pmatrix}y^{(1)}\\\vdots\\y^{(n)}\end{pmatrix}.$$
🔍 Then $(X\beta)_i=x^{(i)T}\beta=\beta^Tx^{(i)}$ (the $i$-th prediction), so $X\beta-z$ is the vector of errors, and
$$\sum_i\big(y^{(i)}-\beta^Tx^{(i)}\big)^2=\|X\beta-z\|^2=(X\beta-z)^T(X\beta-z)=:f(\beta).$$

**Assumption:** the columns of $X$ are **linearly independent** (typical when $n\gg d$, many more data points than features).

#### 3.3 Expanding $f(\beta)$

$$f(\beta)=(\beta^TX^T-z^T)(X\beta-z)\qquad[\text{since }(X\beta-z)^T=\beta^TX^T-z^T]$$
Multiply out (four terms):
$$=\beta^TX^TX\beta-\beta^TX^Tz-z^TX\beta+z^Tz .$$
🔍 The two middle terms are equal: $z^TX\beta$ is a $1\times1$ number, so it equals its transpose $(z^TX\beta)^T=\beta^TX^Tz$. Hence
$$f(\beta)=\beta^T(X^TX)\beta-2(X^Tz)^T\beta+z^Tz .$$
This is a quadratic $\beta^TQ\beta+c^T\beta+d$ with $Q=X^TX$ (symmetric, since $(X^TX)^T=X^TX$), $c=-2X^Tz$, $d=z^Tz$.

#### 3.4 FONC → normal equations

Using $\nabla(\beta^TQ\beta)=2Q\beta$ and $\nabla(c^T\beta)=c$:
$$\nabla f(\beta)=2X^TX\beta-2X^Tz=0\quad\Longrightarrow\quad X^TX\beta=X^Tz\quad(\textbf{normal equations}).$$
If $X^TX$ is invertible:
$$\boxed{\beta^*=(X^TX)^{-1}X^Tz}$$

#### 3.5 Is $X^TX$ invertible?

$\det(X^TX)\neq0\iff X^TX$ full rank $\iff\mathcal N(X^TX)=\{0\}$. Let's prove the last statement.

Take any $q\in\mathcal N(X^TX)$, so $X^TXq=0_{d\times1}$.
1. Multiply on the left by $q^T$: $q^TX^TXq=0$.
2. Recognise $q^TX^TXq=(Xq)^T(Xq)=\|Xq\|^2$. So $\|Xq\|^2=0$.
3. A vector with zero length is the zero vector: $Xq=0$, so $q\in\mathcal N(X)$.
4. But the columns of $X$ are independent. 🔍 $Xq=\sum_jq_j(\text{column }j)=0$ forces all $q_j=0$. So $q=0$.

Therefore $\mathcal N(X^TX)=\{0\}$ and $X^TX$ is invertible ✓.

#### 3.6 SOSC: is it a minimum?

$$\nabla^2f(\beta)=2X^TX .$$
For any $v$: $v^TX^TXv=\|Xv\|^2\ge0$, so $X^TX\succeq0$. Moreover, by §3.5, $\|Xv\|^2=0$ only when $v=0$, so $v^TX^TXv>0$ for $v\neq0$: $X^TX\succ0$ (positive definite).

FONC holds at $\beta^*$ and the Hessian is PD ⇒ by SOSC, $\beta^*$ is a **local min**.

✅ And in fact a **global** min. 🔍 *Preview of why:* the Hessian is PD *everywhere* (it does not depend on $\beta$), which makes $f$ (strictly) **convex**, and for convex functions every local min is global. (Proved in the Convexity blocks.)

---

### 4. Worked example

Fit $y=\beta x$ (one feature, no intercept) to points $(1,2),(2,3),(3,7)$.
- $X=\begin{pmatrix}1\\2\\3\end{pmatrix}$, $z=\begin{pmatrix}2\\3\\7\end{pmatrix}$.
- $X^TX=1+4+9=14$, $X^Tz=2+6+21=29$.
- $\beta^*=29/14\approx2.07$.

With an intercept ($y=\beta_0+\beta_1x$), add a column of ones: $X=\begin{pmatrix}1&1\\1&2\\1&3\end{pmatrix}$, $X^TX=\begin{pmatrix}3&6\\6&14\end{pmatrix}$, $X^Tz=\begin{pmatrix}12\\29\end{pmatrix}$.
- $\det=42-36=6$, $(X^TX)^{-1}=\frac16\begin{pmatrix}14&-6\\-6&3\end{pmatrix}$.
- $\beta^*=\frac16\begin{pmatrix}14\cdot12-6\cdot29\\-6\cdot12+3\cdot29\end{pmatrix}=\frac16\begin{pmatrix}-6\\15\end{pmatrix}=\begin{pmatrix}-1\\2.5\end{pmatrix}$.

So the best line is $y=-1+2.5x$.

---

### 5. Summary

✅ Least squares: $\min_\beta\|X\beta-z\|^2$, where $X$ has one data point per row.
✅ $\nabla f=2X^TX\beta-2X^Tz$; FONC gives the normal equations $X^TX\beta=X^Tz$.
✅ Columns of $X$ independent ⇒ $X^TX\succ0$ ⇒ invertible.
✅ $\beta^*=(X^TX)^{-1}X^Tz$ is the (global) minimiser.

<!-- LECTURE END: linearreg -->

---

<!-- LECTURE START: aug21 -->
<a id="aug21"></a>

## Lecture 10: One-Dimensional Optimization: bisection search (Aug 21 class, partial)

*Source file:* `Aug21.pdf`

> 📌 **This file is the first part of `OneDim.pdf`.** It covers only bisection search. The full explanation is in the block **"One-Dimensional Optimization: bisection, Newton's method, root finding, secant method"**. That block also covers bisection for constrained problems, Newton's method, root finding and the secant method. This block is a short checklist.

### Checklist

1. **Assumptions:** $x^*\in[a_0,b_0]$; $f$ continuous, differentiable and **unimodal** (one valley); $f'$ changes sign once.
2. **Bisection:** $x_t=\frac{a_t+b_t}{2}$. If $|f'(x_t)|<\varepsilon$, return $x_t$. If $f'(x_t)>0$ (going up), keep the left half $[a_t,x_t]$. Otherwise keep the right half $[x_t,b_t]$.
3. **Interval length:** $I_t=\frac{I_0}{2^t}$. Requiring $\frac{I_0}{2^t}<\varepsilon$ gives $T\ge\log_2\frac{I_0}{\varepsilon}$ iterations.

✅ See the full One-Dimensional Optimization block for the derivations and worked examples.

<!-- LECTURE END: aug21 -->

---

<!-- LECTURE START: onedim -->
<a id="onedim"></a>

## Lecture 11: One-Dimensional Optimization: bisection, Newton's method, root finding, secant method

*Source file:* `OneDim.pdf`

### 1. What this lecture is about

Many algorithms reduce to minimising a function of **one** variable (e.g. choosing the step size $\alpha$ in a line search). Here we learn 1-D methods:

1. Assumptions: unimodal $f$, known bracket $[a_0,b_0]$
2. **Bisection search** on $f'$, and how many iterations it needs ($T\ge\log_2\frac{I_0}{\varepsilon}$)
3. Bisection used for **constrained** problems (via the epigraph trick)
4. **Newton's method** for minimisation (quadratic approximation), with derivation
5. Newton's method for **root finding** (linear approximation), and an example where it fails ($x^{1/3}$)
6. **Secant method** (Newton without second derivatives)

---

### 2. Background you need first

- $f'(x)=0$ at a minimum (FONC in 1-D).
- Taylor: $f(x_k+d)\approx f(x_k)+f'(x_k)d+\frac12f''(x_k)d^2$.
- $\log_2$: $2^t<y\iff t<\log_2y$.

---

### 3. Main content with full derivations

#### 3.1 Setting

$\min f(x)$, $f:\mathbb{R}\to\mathbb{R}$. Suppose we know the minimiser lies in an interval: $x^*\in[a_0,b_0]$, and
- $f$ is continuous and differentiable;
- $f$ is **unimodal** on $[a_0,b_0]$: it has a single valley (goes down then up). Wiggly functions with several dips (crossed out in the notes) are excluded.
- Equivalently, $f'(x)$ **changes sign only once** in $[a_0,b_0]$: negative left of $x^*$, positive right of it.

#### 3.2 Bisection search

```
Initialize a_0, b_0, t = 0
loop {
    x_t = (a_t + b_t)/2              # midpoint
    if |f'(x_t)| < ε : return x_t    # nearly flat ⇒ done
    if f'(x_t) > 0 :                 # f increasing at x_t ⇒ minimum is to the LEFT
        b_{t+1} = x_t ; a_{t+1} = a_t
    else :                           # f decreasing ⇒ minimum is to the RIGHT
        a_{t+1} = x_t ; b_{t+1} = b_t
    t = t + 1
}
```
🔍 *Why the rule is right.* If $f'(x_t)>0$, the function is going **up** at $x_t$; since $f$ has only one valley, the bottom must be to the left, so throw away $[x_t,b_t]$. If $f'(x_t)<0$, $f$ is going down, so the bottom is to the right.

**How fast does the interval shrink?** Let $I_t=b_t-a_t$ (interval length).
$$I_0=b_0-a_0,\quad I_1=\frac{I_0}{2},\quad I_2=\frac{I_1}{2}=\frac{I_0}{2^2},\quad\dots,\quad I_t=\frac{I_0}{2^t}.$$
To guarantee the interval (and hence the error in $x^*$) is below $\varepsilon$:
$$\frac{I_0}{2^t}<\varepsilon\iff2^t>\frac{I_0}{\varepsilon}\iff t>\log_2\frac{I_0}{\varepsilon}.$$
So the number of iterations needed is $T\ge\log_2\frac{I_0}{\varepsilon}$.
*Example:* $I_0=1$, $\varepsilon=10^{-6}$: $T\ge\log_2 10^6\approx19.9$, so **20 iterations**.

#### 3.3 Bisection for constrained optimization (epigraph trick)

$$\min_x f(x)\ \text{s.t. } g(x)\le0\qquad\Longleftrightarrow\qquad\min_{x,t}\ t\ \text{s.t. } f(x)\le t,\ g(x)\le0 .$$
🔍 *Why equivalent:* for fixed $x$, the smallest allowed $t$ is $t=f(x)$, so minimising $t$ is the same as minimising $f(x)$.

**Bisection on $t$:** for a candidate level $t$, ask the *feasibility* question "is there $x$ with $f(x)\le t$ and $g(x)\le0$?" (the horizontal lines in the sketch).
- If **yes**, the optimal value is $\le t$: lower the upper bound.
- If **no** (e.g. level $t'$ below the bottom of the curve), the optimal value is $>t'$: raise the lower bound.

Halving the interval of $t$ values each time pins down the optimal value.

#### 3.4 Newton's method for minimisation

Idea: at the current point $x_k$, approximate $f$ by a parabola $g(x)=ax^2+bx+c$ that matches $f$'s value, slope and curvature at $x_k$. Then jump to the bottom of that parabola.

**Matching conditions** ($g'(x)=2ax+b$, $g''(x)=2a$):
- curvature: $2a=f''(x_k)\Rightarrow a=\frac12f''(x_k)$;
- slope: $2ax_k+b=f'(x_k)\Rightarrow b=f'(x_k)-2ax_k=f'(x_k)-f''(x_k)x_k$;
- value: $ax_k^2+bx_k+c=f(x_k)\Rightarrow c=f(x_k)-ax_k^2-bx_k$. (We don't need $c$: it does not affect where the minimum is.)

**Equivalent and simpler: Taylor form.** Writing $x=x_k+d$:
$$g(x)=f(x_k)+f'(x_k)(x-x_k)+\tfrac12f''(x_k)(x-x_k)^2 .$$
Minimise: set $g'(x)=f'(x_k)+f''(x_k)(x-x_k)=0$, so
$$\boxed{x_{k+1}=x_k-\frac{f'(x_k)}{f''(x_k)}}$$
(The bottom of the parabola is a minimum only if $f''(x_k)>0$.)

**Algorithm.**
```
Initialize x_0, k = 0
loop: x_{k+1} = x_k − f'(x_k)/f''(x_k);  k = k+1
until |f'(x_{k+1})| < ε
```
⚠️ The sketches show what can go wrong: where $f''<0$ (near a hump) the parabola opens downward and Newton moves toward a **maximum**; where $f$ is very flat (a plateau), $f''\approx0$ and the step is huge.

#### 3.5 Newton's method for root finding

Problem: find $x$ with $f(x)=0$. Approximate $f$ at $x_k$ by its tangent line $g(x)=ax+b$:
- $g'(x)=a=f'(x_k)$;
- $ax_k+b=f(x_k)\Rightarrow b=f(x_k)-f'(x_k)x_k$.

So $g(x)=f'(x_k)x+f(x_k)-f'(x_k)x_k=f(x_k)+f'(x_k)(x-x_k)$. Set $g(x)=0$:
$$\boxed{x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}}$$
(The next point is where the tangent crosses the axis.)

🔍 **Connection:** minimising $f$ = finding a root of $f'$. Apply root-finding Newton to $f'$ and you get exactly $x_{k+1}=x_k-f'(x_k)/f''(x_k)$.

**Failure example: $f(x)=x^{1/3}$** (root at $0$).
- $f'(x)=\frac13x^{-2/3}$.
- $x_{k+1}=x_k-\dfrac{x_k^{1/3}}{\frac13x_k^{-2/3}}=x_k-3x_k^{1/3}x_k^{2/3}=x_k-3x_k=-2x_k$.

So $|x_{k+1}|=2|x_k|$: each step **doubles** the distance from the root and flips sides ($1\to-2\to4\to-8\dots$). Newton **diverges**, however close you start. (The graph is very steep at 0, so tangents overshoot.)

#### 3.6 Secant method

Newton needs $f''(x_k)$, which may be expensive or unavailable. Replace it by a finite difference of the last two derivatives:
$$f''(x_k)=\lim_{h\to0}\frac{f'(x_k+h)-f'(x_k)}{h}\ \approx\ \frac{f'(x_k)-f'(x_{k-1})}{x_k-x_{k-1}} .$$
Plugging in:
$$x_{k+1}=x_k-f'(x_k)\,\frac{x_k-x_{k-1}}{f'(x_k)-f'(x_{k-1})}.$$
It needs two starting points $x_0,x_1$. It converges a bit slower than Newton but each step is cheaper.

---

### 4. Worked examples

**Bisection.** $f(x)=(x-2)^2$ on $[0,3]$, $f'(x)=2(x-2)$.
- $t=0$: $x_0=1.5$, $f'=-1<0$ ⇒ $[1.5,3]$.
- $t=1$: $x_1=2.25$, $f'=0.5>0$ ⇒ $[1.5,2.25]$.
- $t=2$: $x_2=1.875$, $f'=-0.25$ ⇒ $[1.875,2.25]$. … closing in on 2.

**Newton (min).** $f(x)=x^4-3x$, $f'=4x^3-3$, $f''=12x^2$, $x_0=1$.
- $x_1=1-\frac{1}{12}=0.9167$.
- $x_2=0.9167-\frac{4(0.770)-3}{12(0.840)}=0.9167-\frac{0.0807}{10.08}\approx0.9087$.
- True min: $x^*=(3/4)^{1/3}\approx0.9086$. Converges very fast (quadratically).

---

### 5. Summary

✅ Bisection: halve $[a,b]$ using the sign of $f'$ at the midpoint; $I_t=I_0/2^t$; need $t\ge\log_2(I_0/\varepsilon)$.
✅ Constrained problems → $\min t$ s.t. $f\le t$; bisect on $t$.
✅ Newton (min): $x_{k+1}=x_k-f'(x_k)/f''(x_k)$, from the quadratic Taylor model.
✅ Newton (root): $x_{k+1}=x_k-f(x_k)/f'(x_k)$; can diverge (e.g. $x^{1/3}$).
✅ Secant: replace $f''$ by $\frac{f'(x_k)-f'(x_{k-1})}{x_k-x_{k-1}}$.

<!-- LECTURE END: onedim -->

---

<!-- LECTURE START: apr28-uncon -->
<a id="apr28-uncon"></a>

## Lecture 12: Descent Methods: descent directions, gradient & Newton directions, step sizes, stopping rules

*Source file:* `Apr28-uncon.pdf`

### 1. What this lecture is about

For most functions we cannot solve $\nabla f(x)=0$ by hand. Instead we use **iterative algorithms**: start at a guess $x_0$, then repeatedly improve it, $x_0\to x_1\to x_2\to\cdots$. This lecture gives the general template and its ingredients:

1. The update rule $x_{k+1}=x_k+\alpha_kd_k$ (step size × direction)
2. What a **descent direction** is, and the test $\nabla f(x_k)^Td<0$ (with proof)
3. Choices of direction: steepest descent, scaled gradient $-B_k\nabla f$, Newton, modified Newton, diagonally scaled, quasi-Newton
4. Choices of step size: exact line search, diminishing, constant, inexact line search
5. **Stopping criteria**

---

### 2. Background you need first

- Directional derivative $\nabla f(x)^Td$ and Cauchy–Schwarz (Calculus & LA blocks).
- First-order Taylor: $f(x+\alpha d)=f(x)+\alpha\nabla f(x)^Td+o(\|\alpha d\|)$.
- Positive definite $B\succ0$: $v^TBv>0$ for $v\ne0$.
- $\mathbb{R}_+$ = positive real numbers.

---

### 3. Main content with full derivations

#### 3.1 The general iteration

Problem: $\min f(x)$ s.t. $x\in\mathbb{R}^n$. Generate
$$x_{k+1}=x_k+\alpha_kd_k,\qquad\alpha_k\in\mathbb{R}_+\ (\text{step size}),\quad d_k\in\mathbb{R}^n\ (\text{direction}).$$
Picture: from $x_k$, pick an arrow $d_k$ and walk a distance controlled by $\alpha_k$.

#### 3.2 Descent direction

We want $d_k$ to be a **descent direction**: there is an $r>0$ such that
$$f(x_k+\beta d_k)<f(x_k)\qquad\forall\beta\in(0,r).$$
I.e. moving a *small enough* amount along $d_k$ strictly lowers $f$. (A big step might overshoot, as in the 1-D sketch where the function goes down and then up again.)

#### 3.3 Characterisation: $\nabla f(x_k)^Td<0$ ⇒ descent direction

> **Claim.** Any $d$ with $\nabla f(x_k)^Td<0$ is a descent direction for $f$ at $x_k$.

**Proof, step by step.**
1. First-order Taylor expansion:
$$f(x+\alpha d)=f(x)+\alpha\underbrace{\nabla f(x)^Td}_{<0}+g(\alpha d),\qquad g(\alpha d)=o(\|\alpha d\|).$$
2. Meaning of $o(\|\alpha d\|)$: $\lim_{\|\alpha d\|\to0}\frac{g(\alpha d)}{\|\alpha d\|}=0$. Since $\|\alpha d\|=|\alpha|\,\|d\|$ and $\|d\|$ is a fixed positive number, this is the same as $\lim_{\alpha\to0}\frac{g(\alpha d)}{|\alpha|}=0$.
3. So: $\forall\varepsilon>0\ \exists\delta>0$ such that $g(\alpha d)<\varepsilon|\alpha|$ for all $\alpha\in(0,\delta)$.
4. 🔍 *Missing finish.* Let $c=-\nabla f(x)^Td>0$. Choose $\varepsilon=c/2$. Then for $\alpha\in(0,\delta)$:
$$f(x+\alpha d)-f(x)=-\alpha c+g(\alpha d)<-\alpha c+\tfrac c2\alpha=-\tfrac c2\alpha<0 .$$
So $f$ strictly decreases for every $\alpha\in(0,\delta)$; take $r=\delta$. ∎

**Interpretation.** At $x_k$ the rate of change of $f$ along $d$ is the directional derivative
$$\phi(\alpha)=f(x_k+\alpha d),\quad\phi'(\alpha)=\nabla f(x_k+\alpha d)^Td,\quad\phi'(0)=\nabla f(x_k)^Td .$$
A negative slope at $\alpha=0$ means going downhill.

**Which direction decreases fastest?** By Cauchy–Schwarz,
$$-\|\nabla f(x_k)\|\,\|d\|\ \le\ \nabla f(x_k)^Td\ \le\ \|\nabla f(x_k)\|\,\|d\| .$$
The upper bound is reached when $d=\nabla f(x_k)$ (steepest ascent), the lower bound when $d=-\nabla f(x_k)$ (**steepest descent**).

#### 3.4 Gradient methods (choices of $d_k$)

**(1) Gradient descent / steepest descent:** $d_k=-\nabla f(x_k)$.
Check: $\nabla f(x_k)^Td_k=-\|\nabla f(x_k)\|^2<0$ whenever $\nabla f(x_k)\neq0$ ✓.

**(2) Scaled gradient:** pick $B_k\in\mathbb{R}^{n\times n}$, $B_k\succ0$, and set $d_k=-B_k\nabla f(x_k)$.
Check: $\nabla f(x_k)^Td_k=-\underbrace{\nabla f(x_k)^TB_k\nabla f(x_k)}_{>0\text{ since }B_k\succ0}<0$ ✓. So **any PD scaling gives a descent direction**.

**(3) Newton's method:** $d_k=-\big(\nabla^2f(x_k)\big)^{-1}\nabla f(x_k)$.
This is case (2) with $B_k=(\nabla^2f(x_k))^{-1}$, so it is a descent direction **if $\nabla^2f(x_k)$ is positive definite** (the inverse of a PD matrix is PD). 🔍 *Where it comes from:* minimise the 2nd-order Taylor model $m(d)=f(x_k)+\nabla f(x_k)^Td+\frac12d^T\nabla^2f(x_k)d$; setting $\nabla m(d)=\nabla f(x_k)+\nabla^2f(x_k)d=0$ gives exactly this $d$.

**(4) Modified Newton:** $d_k=-\big(\nabla^2f(x_0)\big)^{-1}\nabla f(x_k)$. Reuse the Hessian at the starting point for every step; cheaper, since the inverse is computed only once.

**(5) Diagonally scaled steepest descent:**
$$d_k=-\begin{pmatrix}\frac{\partial^2f}{\partial x_1^2}&&0\\&\ddots&\\0&&\frac{\partial^2f}{\partial x_n^2}\end{pmatrix}^{-1}_{x=x_k}\nabla f(x_k).$$
Keep only the diagonal of the Hessian (easy to invert: just $1/$entry). Good when the variables have very different scales.

**(6) Quasi-Newton methods:** approximate $\big(\nabla^2f(x_k)\big)^{-1}$ by "some" matrix $T_k$ built from gradient information (e.g. BFGS), avoiding second derivatives.

#### 3.5 How to fix $\alpha_k$ (step size)

1. **Exact line search:** $\alpha_k=\arg\min_{\alpha\ge0}f(x_k+\alpha d_k)$. The best step along the line (a 1-D minimisation each iteration).
2. **Diminishing step sizes:** $\alpha_k\to0$ but $\sum_{k=1}^\infty\alpha_k=\infty$, e.g. $\alpha_k=\frac1k$. 🔍 Why the sum must diverge: the total distance you can travel is at most $\sum\alpha_k\|d_k\|$. If $\sum\alpha_k<\infty$ you might stop before reaching the minimum.
3. **Constant step size:** $\alpha_k=s$ for all $k$. Simple, but $s$ must be small enough or the method may diverge.
4. **Inexact line search:** don't find the exact minimum along the line, just accept an $\alpha$ that gives "enough" decrease. The sketch shows $\phi(\alpha)=f(x_k+\alpha d_k)$ and lines through $\phi(0)$ with slopes $\phi'(0)$ and a flatter fraction of it; accept $\alpha$ where $\phi$ lies below the flatter line (the **Armijo** idea; see the Backtracking Line Search block).

(The notes also ponder "$\alpha_k=|\nabla f(x_k)^Td_k|$?". It is a possible heuristic, but not a standard rule.)

#### 3.6 When to stop

Choose a small tolerance $\varepsilon>0$ and stop when one of these holds:
1. $\|\nabla f(x_k)\|<\varepsilon$ (nearly stationary: FONC almost satisfied);
2. $|f(x_k)-f(x_{k+1})|<\varepsilon$ (function value stopped changing);
3. $\|x_k-x_{k+1}\|<\varepsilon$ (iterates stopped moving);
4. $\dfrac{|f(x_k)-f(x_{k+1})|}{\max(1,|f(x_k)|)}<\varepsilon$ (**relative** change; scale-free);
5. $\dfrac{\|x_k-x_{k+1}\|}{\max(1,\|x_k\|)}<\varepsilon$ (relative step).

🔍 The $\max(1,\cdot)$ avoids dividing by numbers near zero while still giving a relative test when values are large.

---

### 4. Worked example

$f(x)=x_1^2+10x_2^2$, $x_0=(1,1)$, $\nabla f=(2x_1,20x_2)$, so $\nabla f(x_0)=(2,20)$.
- Steepest descent $d_0=(-2,-20)$; $\nabla f^Td_0=-4-400=-404<0$ ✓.
- Newton: $\nabla^2f=\operatorname{diag}(2,20)$, $d_0=-(\tfrac22,\tfrac{20}{20})=(-1,-1)$. With $\alpha=1$: $x_1=(0,0)$, the exact minimum in one step (for quadratics Newton is exact).
- Here the diagonally-scaled direction equals Newton's, because the Hessian is diagonal.

---

### 5. Summary

✅ Iterate $x_{k+1}=x_k+\alpha_kd_k$.
✅ $\nabla f(x_k)^Td_k<0$ ⇒ descent direction (Taylor + small-o argument).
✅ $d=-B\nabla f$ with $B\succ0$ always descends; $B=I$: gradient descent; $B=(\nabla^2f)^{-1}$: Newton.
✅ Step sizes: exact, diminishing ($\sum\alpha_k=\infty$), constant, inexact (Armijo).
✅ Stop on small gradient / small change (absolute or relative).

<!-- LECTURE END: apr28-uncon -->

---

<!-- LECTURE START: graddescent-exactls -->
<a id="graddescent-exactls"></a>

## Lecture 13: Steepest Descent with Exact Line Search on a Quadratic

*Source file:* `Graddescent-exactls.pdf`

### 1. What this lecture is about

We take the simplest important function, a **quadratic** $f(x)=\frac12x^TQx+b^Tx+c$, and:
1. find its exact minimiser $x^*=-Q^{-1}b$;
2. run **steepest descent** $x_{k+1}=x_k-\alpha_k\nabla f(x_k)$;
3. derive the **exact line-search step** in closed form: $\alpha_k=\dfrac{\nabla f_k^T\nabla f_k}{\nabla f_k^TQ\nabla f_k}$;
4. work a numerical example.

---

### 2. Background you need first

- $\nabla(\frac12x^TQx)=Qx$ for symmetric $Q$ (because $\nabla(x^TQx)=2Qx$, times $\frac12$). $\nabla(b^Tx)=b$.
- Minimum of a 1-D parabola $A\alpha^2+B\alpha+C$ with $A>0$: at $\alpha=-\frac{B}{2A}$ (set derivative $2A\alpha+B=0$).
- $Q\succ0$ ⇒ $v^TQv>0$ for $v\neq0$, and $Q$ is invertible.

---

### 3. Main content with full derivations

#### 3.1 The problem and its exact solution

$$\min_{x\in\mathbb{R}^n}f(x)=\tfrac12x^TQx+b^Tx+c,\qquad Q\in\mathbb{R}^{n\times n}\ \text{symmetric},\ Q\succ0,\ b\in\mathbb{R}^n,\ c\in\mathbb{R}.$$
- Gradient: $\nabla f(x)=Qx+b$. Hessian: $Q\succ0$.
- FONC: $Qx+b=0\Rightarrow x^*=-Q^{-1}b$. The Hessian is PD (everywhere), so this is the unique (global) minimiser.

(Why the $\frac12$? It cancels the 2 from differentiation, so the gradient is simply $Qx+b$.)

We study the iterative method anyway, because for big $n$ inverting $Q$ is expensive, and because quadratics show how gradient descent behaves near any minimum (where $f$ looks quadratic by Taylor).

#### 3.2 Steepest descent with exact line search

Write $\nabla f_k:=\nabla f(x_k)=Qx_k+b$.
$$x_{k+1}=x_k-\alpha_k\nabla f_k,\qquad\alpha_k=\arg\min_{\alpha\ge0}\ \phi(\alpha),\quad\phi(\alpha)=f(x_k-\alpha\nabla f_k).$$

#### 3.3 Deriving $\alpha_k$ in closed form

🔍 *Full expansion (the notes write "( ... )" for the messy terms).* Let $g=\nabla f_k$ to keep notation short.
$$\phi(\alpha)=\tfrac12(x_k-\alpha g)^TQ(x_k-\alpha g)+b^T(x_k-\alpha g)+c.$$
Expand the quadratic part:
$$(x_k-\alpha g)^TQ(x_k-\alpha g)=x_k^TQx_k-2\alpha\,g^TQx_k+\alpha^2g^TQg$$
(the two cross terms $-\alpha x_k^TQg$ and $-\alpha g^TQx_k$ are equal because $Q$ is symmetric). So
$$\phi(\alpha)=\underbrace{\tfrac12g^TQg}_{A>0}\,\alpha^2\;+\;\underbrace{\big(-g^TQx_k-b^Tg\big)}_{B}\,\alpha\;+\;\underbrace{\tfrac12x_k^TQx_k+b^Tx_k+c}_{C=f(x_k)} .$$
Simplify $B$: $-g^TQx_k-g^Tb=-g^T(Qx_k+b)=-g^Tg$.

So $\phi(\alpha)=\frac12(g^TQg)\alpha^2-(g^Tg)\alpha+f(x_k)$, an upward parabola in $\alpha$ (since $g^TQg>0$). Minimise: $\phi'(\alpha)=(g^TQg)\alpha-g^Tg=0$ gives
$$\boxed{\alpha_k=\frac{\nabla f_k^T\nabla f_k}{\nabla f_k^TQ\,\nabla f_k}}$$
which is positive ✓.

🔍 **A property worth knowing:** with exact line search, consecutive gradients are **orthogonal**, $\nabla f_{k+1}^T\nabla f_k=0$. Reason: at the best $\alpha$, $\phi'(\alpha_k)=0$, and $\phi'(\alpha)=-\nabla f(x_k-\alpha g)^Tg$, so $\nabla f_{k+1}^Tg=0$. That is why steepest descent **zig-zags** at right angles.

#### 3.4 Worked example from the notes

$$f(x)=5x_1^2+x_2^2+4x_1x_2-6x_1-4x_2+15 .$$
**Step 1: write it in the form $\frac12x^TQx+b^Tx+c$.**
- $\frac12x^TQx=\frac12(Q_{11}x_1^2+2Q_{12}x_1x_2+Q_{22}x_2^2)$. Match: $\frac12Q_{11}=5\Rightarrow Q_{11}=10$; $Q_{12}=4$; $\frac12Q_{22}=1\Rightarrow Q_{22}=2$.
$$Q=\begin{pmatrix}10&4\\4&2\end{pmatrix},\quad b=\begin{pmatrix}-6\\-4\end{pmatrix},\quad c=15 .$$
- Is $Q\succ0$? Leading minors: $10>0$, $\det=20-16=4>0$ ✓.

**Step 2: exact minimiser** (🔍 not in the notes). $Q^{-1}=\frac14\begin{pmatrix}2&-4\\-4&10\end{pmatrix}$.
$$x^*=-Q^{-1}b=-\tfrac14\begin{pmatrix}2(-6)-4(-4)\\-4(-6)+10(-4)\end{pmatrix}=-\tfrac14\begin{pmatrix}4\\-16\end{pmatrix}=\begin{pmatrix}-1\\4\end{pmatrix},$$
with $f(x^*)=5+16-16+6-16+15=10$. Check: $\nabla f(x^*)=Qx^*+b=(-10+16-6,\ -4+8-4)=(0,0)$ ✓.

**Step 3: one iteration of steepest descent from $x_0=(0,0)$** (🔍 not in the notes).
- $\nabla f_0=Qx_0+b=(-6,-4)$.
- $\nabla f_0^T\nabla f_0=36+16=52$.
- $Q\nabla f_0=(10(-6)+4(-4),\ 4(-6)+2(-4))=(-76,-32)$; $\nabla f_0^TQ\nabla f_0=(-6)(-76)+(-4)(-32)=456+128=584$.
- $\alpha_0=52/584\approx0.0890$.
- $x_1=x_0-\alpha_0\nabla f_0=(0.534,\ 0.356)$.
- $\nabla f_1=Qx_1+b\approx(5.342+1.425-6,\ 2.137+0.712-4)=(0.767,\ -1.151)$.
- Orthogonality check: $(-6)(0.767)+(-4)(-1.151)\approx-4.60+4.60=0$ ✓.

Continuing, the iterates zig-zag toward $(-1,4)$. The ratio of eigenvalues of $Q$ ($\approx11.7/0.34\approx34$) is large, so progress is slow: the **condition number** effect (see Convergence Results).

---

### 4. Summary

✅ Quadratic $\frac12x^TQx+b^Tx+c$, $Q\succ0$: $\nabla f=Qx+b$, $x^*=-Q^{-1}b$.
✅ Exact line search along $-\nabla f_k$: $\alpha_k=\frac{\nabla f_k^T\nabla f_k}{\nabla f_k^TQ\nabla f_k}$.
✅ Successive gradients are orthogonal ⇒ zig-zag path.

<!-- LECTURE END: graddescent-exactls -->

---

<!-- LECTURE START: rates-of-convergence -->
<a id="rates-of-convergence"></a>

## Lecture 14: Rates of Convergence: definitions only (partial notes)

*Source file:* `Rates of Convergence..pdf`

> 📌 **This file is the first few lines of the Sep 2 notes (`Sep2.pdf`).** The full lecture, with the limit-ratio test, three worked examples ($1/k$, $a^k$, $a^{2^k}$) and the proof of why the ratio test works, is in the block **"Rates of Convergence (full lecture)"**. This block is a quick checklist.

### Checklist

- **Definition:** $\{x_k\}\to x^*$ has order $p$ and factor $\gamma>0$ if $\exists N$ such that $\|x_{k+1}-x^*\|\le\gamma\|x_k-x^*\|^p$ for all $k>N$.
- $p=1,\ \gamma<1$: **linear**. $p=2$: **quadratic**. $p>1$: **superlinear**. $p=1,\ \gamma=1$: **sublinear**.
- ✅ Desirable: **large $p$, small $\gamma$**. Squaring an error of $0.1$ gives $0.01$, which is why higher orders are so much faster.

<!-- LECTURE END: rates-of-convergence -->

---

<!-- LECTURE START: sep2 -->
<a id="sep2"></a>

## Lecture 15: Rates of Convergence (full lecture): order $p$, factor $\gamma$, the limit-ratio test, examples

*Source file:* `Sep2.pdf`

### 1. What this lecture is about

Two algorithms can both converge to $x^*$ while one needs 10 iterations and the other a million. **Order of convergence** measures *how fast* the error shrinks per iteration. We:

1. Define convergence of order $p$ with factor $\gamma$
2. Name the cases: linear, superlinear, quadratic, sublinear
3. Learn the practical test: compute $\lim_{k\to\infty}\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|^p}$
4. Work three examples: $\frac1k$ (sublinear), $a^k$ (linear), $a^{2^k}$ (quadratic)
5. Prove why the limit of the ratio gives the definition

---

### 2. Background you need first

- Error at iteration $k$: $e_k=\|x_k-x^*\|$ (distance from the answer).
- Limits: $\lim_{k\to\infty}\frac{k}{k+1}=1$; for $0<a<1$, $a^k\to0$ and $a^{-k}\to\infty$.
- $\varepsilon$–$N$ definition of a limit (Calculus block): $\lim_{k\to\infty}g_k=\gamma\iff\forall\varepsilon>0\ \exists N:\ |g_k-\gamma|<\varepsilon\ \forall k>N$.
- Exponent rules: $\frac{a^m}{a^n}=a^{m-n}$, $(a^m)^p=a^{mp}$.

---

### 3. Main content with full derivations

#### 3.1 Definition

Let $\{x_k\}\to x^*$. Convergence is **of order $p$ with factor $\gamma>0$** if
$$\exists N\ \text{such that}\quad\|x_{k+1}-x^*\|\le\gamma\,\|x_k-x^*\|^p\qquad\forall k>N .$$

| Case | Name |
|---|---|
| $p=1,\ \gamma<1$ | **linear** convergence |
| $p>1$ | **superlinear** convergence |
| $p=2$ | **quadratic** convergence |
| $p=1,\ \gamma=1$ | **sublinear** convergence |

Note: $\|x_{k+1}-x^*\|\le\|x_k-x^*\|$ (i.e. $p=1$, $\gamma=1$) is the weakest statement; "always" true for a method that never gets worse. ✅ **Desirable: large $p$, small $\gamma$.** Side note in the notes: a squared error $0.1\to0.01$ shows how squaring speeds things up.

#### 3.2 The practical test

Compute
$$\boxed{\ \lim_{k\to\infty}\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|^p}\ }$$
for candidate values of $p$. The order is the **largest $p$ for which this limit is finite** (and, for $p=1$, the limit is the rate $\gamma$).

#### 3.3 Example 1: $\|x_k-x^*\|=\frac1k$

- $p=1$: $\displaystyle\lim_{k\to\infty}\frac{1/(k+1)}{1/k}=\lim_{k\to\infty}\frac{k}{k+1}=1=\gamma$.
- $p=2$: $\displaystyle\lim_{k\to\infty}\frac{1/(k+1)}{1/k^2}=\lim_{k\to\infty}\frac{k^2}{k+1}=\infty$.

So $p=1$ with $\gamma=1$: **sublinear**. (🔍 $\frac{k}{k+1}=\frac{1}{1+1/k}\to1$; $\frac{k^2}{k+1}\approx k\to\infty$.)

#### 3.4 Example 2: $\|x_k-x^*\|=a^k$, $0<a<1$ (e.g. $0.5^k$)

- $p=1$: $\displaystyle\lim\frac{a^{k+1}}{a^k}=a<1$ ⇒ **linear convergence** with $\gamma=a$.
- $p=2$: $\displaystyle\lim\frac{a^{k+1}}{(a^k)^2}=\lim\frac{a^{k+1}}{a^{2k}}=\lim a^{1-k}=\lim\frac{a}{a^k}=\infty$ (since $a^k\to0$).

So linear, not quadratic.

#### 3.5 Example 3: $\|x_k-x^*\|=a^{2^k}$, $0<a<1$

E.g. $a=0.1$: errors $(0.1)^2,(0.1)^4,(0.1)^8,\dots$. Each is the square of the previous one.

General $p$:
$$\lim_{k\to\infty}\frac{a^{2^{k+1}}}{\big(a^{2^k}\big)^p}=\lim_{k\to\infty}\frac{a^{2^k\cdot2}}{a^{2^k\cdot p}}=\lim_{k\to\infty}a^{2^k(2-p)} .$$
🔍 (Used $2^{k+1}=2\cdot2^k$ and $(a^m)^p=a^{mp}$.)
- **$p=2$:** exponent $2^k\cdot0=0$, so $a^0=1$. A finite non-zero limit ⇒ **quadratic convergence** with $\gamma=1$.
- **$p=1$:** $\lim a^{2^k}=0$. The ratio goes to 0, so it is "better than linear" (superlinear).
- 🔍 **$p>2$:** exponent $2^k(2-p)\to-\infty$, so $a^{(\text{very negative})}\to\infty$. Not of order $p>2$.

So the order is exactly 2: **only $p=2$ gives a finite non-zero limit.**

#### 3.6 Why do we look at the limit of the ratio?

Suppose
$$\lim_{k\to\infty}\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|}=\gamma .$$
By the $\varepsilon$–$N$ definition: $\forall\varepsilon>0\ \exists N>0$ such that
$$\left|\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|}-\gamma\right|<\varepsilon\qquad\forall k>N .$$
Let $g=\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|}-\gamma$ (the quantity inside $|\cdot|$). Two cases:
- If $g\ge0$: $|g|=g<\varepsilon$ ⇒ $\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|}<\gamma+\varepsilon$ ⇒ $\|x_{k+1}-x^*\|<(\gamma+\varepsilon)\|x_k-x^*\|$.
- If $g<0$: then the ratio is $<\gamma<\gamma+\varepsilon$ directly. (The notes do this as: $-\varepsilon<0<\gamma-\text{ratio}<\varepsilon$, which rearranges to $(-\gamma-\varepsilon)\|x_k-x^*\|<-\|x_{k+1}-x^*\|$, i.e. the same inequality.)

So effectively:
$$\forall\varepsilon>0\ \exists N>0:\quad\|x_{k+1}-x^*\|<(\gamma+\varepsilon)\|x_k-x^*\|\quad\forall k>N .$$
This is the definition with factor $\gamma+\varepsilon$ for any tiny $\varepsilon$. So the limit value $\gamma$ is (essentially) the convergence factor. If $\gamma<1$, pick $\varepsilon$ with $\gamma+\varepsilon<1$ to get genuine linear convergence. The same argument works with $\|x_k-x^*\|^p$ in the denominator.

---

### 4. More intuition

🔍 Starting from error 0.1:
- Linear, $\gamma=0.5$: $0.1,0.05,0.025,\dots$ (about 23 steps to reach $10^{-8}$).
- Quadratic: $0.1,10^{-2},10^{-4},10^{-8}$ (3 steps).
- Sublinear ($1/k$): about $10^8$ steps.

Unrolling linear: $e_k\le\gamma^ke_0$ (geometric decay). Unrolling quadratic ($\gamma=1$): $e_k\le e_0^{2^k}$.

---

### 5. Summary

✅ Order $p$, factor $\gamma$: $\|x_{k+1}-x^*\|\le\gamma\|x_k-x^*\|^p$ eventually.
✅ Test: the limit of $\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|^p}$. The order is the largest $p$ with a finite limit.
✅ $1/k$: sublinear. $a^k$: linear ($\gamma=a$). $a^{2^k}$: quadratic.
✅ A limit $\gamma$ of the ratio ⇒ bound with factor $\gamma+\varepsilon$ for every $\varepsilon>0$.

<!-- LECTURE END: sep2 -->

---

<!-- LECTURE START: convergenceresults -->
<a id="convergenceresults"></a>

## Lecture 16: Convergence Results for Gradient Descent and Newton's Method

*Source file:* `ConvergenceResults.pdf`

### 1. What this lecture is about

We now know the algorithms (gradient descent, Newton) and how to measure speed (order $p$, factor $\gamma$). This lecture states **how fast they actually converge**:

1. Steepest descent + exact line search on a quadratic: **linear** convergence, with rate $\frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\max}+\lambda_{\min}}=\frac{\kappa-1}{\kappa+1}$ ($\kappa$ = condition number)
2. Gradient descent with a **fixed** step $\alpha$ on a quadratic: converges for every start **iff** $0<\alpha<\frac{2}{\lambda_{\max}}$ (proved fully below)
3. Newton's method: **quadratic** convergence near $x^*$ (local result)

---

### 2. Background you need first

- Quadratic $f(x)=\frac12x^TQx+b^Tx+c$, $Q\succ0$; $\nabla f=Qx+b$; $x^*=-Q^{-1}b$ (Exact line search block).
- Eigen-decomposition of symmetric $Q$: orthonormal eigenvectors $v_i$, eigenvalues $\lambda_i>0$. For any vector $e=\sum_i c_iv_i$: $Qe=\sum_i\lambda_ic_iv_i$ and $\|e\|^2=\sum c_i^2$.
- Rates of convergence (Sep 2 block).
- $C^3$ = three times continuously differentiable.

---

### 3. Main content with full derivations

#### 3.1 Result 1: steepest descent with exact line search (quadratic $f$)

> For $f(x)=\frac12x^TQx+b^Tx+c$ with $Q\succ0$, steepest descent with exact line search converges **linearly** to the global minimum from **any** $x_0$. With $\lambda_{\max},\lambda_{\min}>0$ the largest and smallest eigenvalues of $Q$:
> $$\lim_{k\to\infty}\frac{\|x_{k+1}-x^*\|}{\|x_k-x^*\|}\ \le\ \frac{\lambda_{\max}-\lambda_{\min}}{\lambda_{\max}+\lambda_{\min}}=\frac{\kappa-1}{\kappa+1},\qquad\kappa=\frac{\lambda_{\max}}{\lambda_{\min}} .$$
> (Strictly, the classical bound, Kantorovich's inequality, is stated in the $Q$-norm / for $f(x_k)-f^*$. The meaning is the same.)

🔍 *Dividing top and bottom by $\lambda_{\min}$* gives the second form: $\frac{\lambda_{\max}/\lambda_{\min}-1}{\lambda_{\max}/\lambda_{\min}+1}$. Here $\kappa=\lambda_{\max}/\lambda_{\min}\ge1$ is the **condition number** of $Q$.

**Interpretation.**
- $\kappa=1$ (all eigenvalues equal, circular contours): rate $0$, so the minimum is reached in **one** step.
- $\kappa$ large (long thin elliptical contours): rate $\approx1-\frac2\kappa$, close to 1, so convergence is **very slow** (the zig-zag). E.g. $\kappa=100$ gives rate $\approx0.98$: about 115 iterations per extra digit of accuracy.

#### 3.2 Result 2: fixed step size

> Gradient descent $x_{k+1}=x_k-\alpha\nabla f(x_k)$ with a fixed $\alpha$ on the quadratic above converges to $x^*$ from every $x_0$ **if and only if**
> $$0<\alpha<\frac{2}{\lambda_{\max}} .$$

🔍 *Full proof (not in notes).*
1. Error $e_k=x_k-x^*$. Since $Qx^*+b=0$, we have $\nabla f(x_k)=Qx_k+b=Qx_k-Qx^*=Qe_k$.
2. Then $e_{k+1}=x_{k+1}-x^*=x_k-\alpha Qe_k-x^*=(I-\alpha Q)e_k$, and so $e_k=(I-\alpha Q)^ke_0$.
3. Expand $e_0=\sum_ic_iv_i$ in eigenvectors of $Q$. Since $(I-\alpha Q)v_i=(1-\alpha\lambda_i)v_i$:
$$e_k=\sum_i(1-\alpha\lambda_i)^kc_iv_i .$$
4. $e_k\to0$ for **every** $e_0$ (every choice of $c_i$) $\iff|1-\alpha\lambda_i|<1$ for all $i$.
5. $|1-\alpha\lambda_i|<1\iff-1<1-\alpha\lambda_i<1\iff0<\alpha\lambda_i<2\iff0<\alpha<\frac2{\lambda_i}$.
6. This must hold for all $i$; the tightest is the largest eigenvalue: $0<\alpha<\frac{2}{\lambda_{\max}}$. ∎

The rate is $\max_i|1-\alpha\lambda_i|$ — linear. 🔍 The best fixed step is $\alpha=\frac{2}{\lambda_{\max}+\lambda_{\min}}$, which again gives rate $\frac{\kappa-1}{\kappa+1}$.

#### 3.3 Result 3: Newton's method converges quadratically (locally)

> If $f\in C^3$ and $\nabla^2f(x^*)$ is **invertible** (at a point with $\nabla f(x^*)=0$), then there is a ball around $x^*$ such that, starting inside it, Newton's method
> $$x_{k+1}=x_k-\big(\nabla^2f(x_k)\big)^{-1}\nabla f(x_k)$$
> converges to $x^*$ **quadratically**: $\|x_{k+1}-x^*\|\le\gamma\|x_k-x^*\|^2$.

🔍 *Why (1-D sketch).* Let $e_k=x_k-x^*$. Taylor-expand $f'$ around $x_k$, evaluated at $x^*$:
$$0=f'(x^*)=f'(x_k)-f''(x_k)e_k+\tfrac12f'''(\xi)e_k^2 .$$
Divide by $f''(x_k)$ and rearrange: $e_k-\frac{f'(x_k)}{f''(x_k)}=\frac{f'''(\xi)}{2f''(x_k)}e_k^2$. The left side is exactly $x_{k+1}-x^*=e_{k+1}$. So
$$|e_{k+1}|=\Big|\frac{f'''(\xi)}{2f''(x_k)}\Big|\,e_k^2\le\gamma\,e_k^2 ,$$
where $\gamma$ is bounded near $x^*$ because $f'''$ is continuous ($C^3$) and $f''(x^*)\ne0$ (invertibility). Hence quadratic convergence.

⚠️ This is a **local** result. Far from $x^*$, Newton can diverge (recall $x^{1/3}$) or head to a maximum. In practice one adds a line search ("damped Newton").

---

### 4. Worked example

$Q=\begin{pmatrix}10&4\\4&2\end{pmatrix}$ from the exact-line-search example. Eigenvalues: $\lambda^2-12\lambda+4=0\Rightarrow\lambda=6\pm\sqrt{32}$, i.e. $\lambda_{\max}\approx11.66$, $\lambda_{\min}\approx0.343$.
- $\kappa\approx34$; rate $\frac{\kappa-1}{\kappa+1}\approx0.943$, so slow.
- Fixed step must satisfy $\alpha<2/11.66\approx0.1715$. With $\alpha=0.2$ the iterates blow up.

---

### 5. Summary

✅ GD + exact line search on a quadratic: linear, rate $\le\frac{\kappa-1}{\kappa+1}$; bad conditioning means slow.
✅ GD fixed step: converges iff $0<\alpha<2/\lambda_{\max}$ (error $e_{k+1}=(I-\alpha Q)e_k$).
✅ Newton: quadratic convergence in a neighbourhood of $x^*$ if $f\in C^3$ and $\nabla^2f(x^*)$ is invertible.

<!-- LECTURE END: convergenceresults -->

---

<!-- LECTURE START: backtracking-line-search -->
<a id="backtracking-line-search"></a>

## Lecture 17: Backtracking Line Search (Armijo condition)

*Source file:* `Backtracking line search.pdf`

### 1. What this lecture is about

Exact line search ($\alpha_k=\arg\min_\alpha f(x_k+\alpha d_k)$) is usually too expensive: it is a whole 1-D optimization at every step. **Backtracking** finds a step that is *good enough* very cheaply: start with a big step and shrink it until the function has decreased "sufficiently". This lecture explains:

1. The 1-D function $\phi(\alpha)=f(x_k+\alpha d_k)$ and its slope at 0
2. The tangent line $g(\alpha)$ and the flatter **Armijo line** $g_1(\alpha)=\phi(0)+c\,\alpha\,\phi'(0)$, $c\in(0,0.5)$
3. The backtracking algorithm and why it always stops

---

### 2. Background you need first

- Descent direction: $\nabla f(x_k)^Td_k<0$ (Descent Methods block).
- Directional derivative via chain rule: $\frac{d}{d\alpha}f(x_k+\alpha d_k)=\nabla f(x_k+\alpha d_k)^Td_k$.
- First-order Taylor: $\phi(\alpha)=\phi(0)+\alpha\phi'(0)+o(\alpha)$.

---

### 3. Main content with full derivations

#### 3.1 The 1-D slice

Fix the current point $x_k$ and a descent direction $d_k$ ($\nabla f(x_k)^Td_k<0$). Define
$$\phi(\alpha)=f(x_k+\alpha d_k),\qquad\phi(0)=f(x_k).$$
- $\phi'(\alpha)=\nabla f(x_k+\alpha d_k)^Td_k$ (chain rule),
- $\phi'(0)=\nabla f(x_k)^Td_k<0$: the slice starts going **down**.

#### 3.2 Two lines through $(0,\phi(0))$

- **Tangent line:** $g(\alpha)=\phi(0)+\alpha\,\phi'(0)$, slope $\phi'(0)$ (e.g. $-1$).
- **Armijo (sufficient-decrease) line:** a *flatter* line
$$g_1(\alpha)=\phi(0)+c\,\alpha\,\phi'(0),\qquad c\in(0,0.5).$$
With $c=0.4$ and $\phi'(0)=-1$, its slope is $c\phi'(0)=-0.4$.

**Accept $\alpha$ if** the curve is below the Armijo line:
$$\boxed{f(x_k+\alpha d_k)\le f(x_k)+c\,\alpha\,\nabla f(x_k)^Td_k}\quad(\textbf{Armijo condition}).$$
The hatched region in the sketch (between curve and line) shows the acceptable $\alpha$'s. The condition says: the actual decrease is at least a fraction $c$ of the decrease predicted by the tangent.

🔍 *Why acceptable steps always exist near 0.* By Taylor, $\phi(\alpha)=\phi(0)+\alpha\phi'(0)+o(\alpha)$. So
$$\phi(\alpha)-g_1(\alpha)=(1-c)\,\alpha\,\phi'(0)+o(\alpha).$$
Since $(1-c)\phi'(0)<0$, for small $\alpha>0$ the negative term $(1-c)\alpha\phi'(0)$ dominates $o(\alpha)$ (same argument as in the descent-direction proof). So $\phi(\alpha)<g_1(\alpha)$ for all sufficiently small $\alpha$. The tangent itself ($c=1$) would not work: for a convex-looking $\phi$, the curve lies above its tangent.

#### 3.3 The algorithm

```
Initialize α_0 (e.g. 1), choose c ∈ (0, 0.5), ρ ∈ (0,1) (e.g. 0.5);  j = 0
while  f(x_k + α_j d_k) > f(x_k) + c α_j ∇f(x_k)ᵀ d_k :     # not enough decrease
       α_{j+1} = ρ α_j ;   j = j + 1                         # shrink ("backtrack")
return α_j
```
(The notes use the same letter $k$ for this inner counter; here it is $j$ to avoid confusion with the outer iteration.)

🔍 *It terminates.* After $j$ shrinks, $\alpha_j=\rho^j\alpha_0\to0$, and by §3.2 the condition holds once $\alpha_j$ is small enough. So the loop ends after finitely many steps.

🔍 *Why not accept tiny steps immediately?* Starting from a large $\alpha_0$ and shrinking only as needed gives the **largest** step of the form $\rho^j\alpha_0$ that passes. That prevents the method from crawling with needlessly small steps.

---

### 4. Worked example

$f(x)=x^2$, $x_k=1$, $d_k=-f'(1)=-2$, $c=0.4$, $\rho=0.5$, $\alpha_0=1$.
- $\phi(\alpha)=(1-2\alpha)^2$, $\phi(0)=1$, $\phi'(0)=f'(1)d_k=2(-2)=-4$.
- Armijo line: $1-1.6\alpha$.
- $\alpha=1$: $\phi=(1-2)^2=1$; line $=-0.6$. $1>-0.6$ ⇒ reject.
- $\alpha=0.5$: $\phi=0$; line $=0.2$. $0\le0.2$ ⇒ **accept**. (Lucky: this is the exact minimiser.)

---

### 5. Summary

✅ $\phi(\alpha)=f(x_k+\alpha d_k)$, $\phi'(0)=\nabla f(x_k)^Td_k<0$.
✅ Armijo: accept $\alpha$ if $f(x_k+\alpha d_k)\le f(x_k)+c\alpha\nabla f(x_k)^Td_k$, $c\in(0,0.5)$.
✅ Backtracking: start big, multiply by $\rho$ until Armijo holds; always terminates.

<!-- LECTURE END: backtracking-line-search -->

---

<!-- LECTURE START: convexitysep7 -->
<a id="convexitysep7"></a>

## Lecture 18: Convexity (Part 1): convex sets, convex hull, convex functions, why convex optimization is "easy"

*Source file:* `ConvexitySep7.pdf`

### 1. What this lecture is about

In general optimization a local minimum need not be the global one (think of a landscape with many valleys). **Convexity** is the property that removes this difficulty: for convex problems, **every local minimum is a global minimum**. This lecture:

1. Defines convex sets, with examples (hyperplane, half-space, ball, PSD matrices)
2. Defines convex combinations and the convex hull
3. Defines convex functions (the "chord lies above the graph" picture)
4. Proves two key theorems for convex optimization problems:
   - every local optimum is global;
   - the set of global optima is convex
5. Lists operations that preserve convexity of sets (intersection, Minkowski sum); union does **not**
6. Shows that linear programs are convex problems

---

### 2. Background you need first

- Line segment between $x$ and $y$: all points $\lambda x+(1-\lambda)y$ with $\lambda\in[0,1]$. ($\lambda=1$ gives $x$, $\lambda=0$ gives $y$, $\lambda=\frac12$ the midpoint.)
- **Triangle inequality:** $\|u+v\|\le\|u\|+\|v\|$. **Homogeneity:** $\|\lambda u\|=|\lambda|\,\|u\|$.
- Local vs global minimum (Aug 17 block).

---

### 3. Main content with full derivations

#### 3.1 Convex sets

$S\subseteq\mathbb{R}^n$ is **convex** if
$$x,y\in S\ \Rightarrow\ \lambda x+(1-\lambda)y\in S\qquad\forall\lambda\in[0,1].$$
In words: the straight segment joining any two points of $S$ stays inside $S$. A disc is convex. A "bow-tie"/star shape is not: some segments leave the set.

**Example 1: hyperplane** $S=\{x\in\mathbb{R}^n:a^Tx=b\}$ (a line in 2-D, a plane in 3-D).
Take $x_1,x_2\in S$, so $a^Tx_1=b$ and $a^Tx_2=b$. Let $z=\lambda x_1+(1-\lambda)x_2$:
$$a^Tz=\lambda a^Tx_1+(1-\lambda)a^Tx_2=\lambda b+(1-\lambda)b=b .$$
So $z\in S$ ✓ convex.

**Example 2: half-space** $S=\{x:a^Tx\le b\}$ (one side of a hyperplane).
🔍 *Proof (not written):* $a^Tz=\lambda a^Tx_1+(1-\lambda)a^Tx_2\le\lambda b+(1-\lambda)b=b$, using $\lambda\ge0$ and $1-\lambda\ge0$ so the inequalities are preserved when multiplied.

**Example 3: ball** $S=\{x:\|x-x_0\|\le r\}$.
Take $x_1,x_2\in S$: $\|x_1-x_0\|\le r$, $\|x_2-x_0\|\le r$. Let $z=\lambda x_1+(1-\lambda)x_2$.
$$\|z-x_0\|=\|\lambda x_1+(1-\lambda)x_2-\lambda x_0-(1-\lambda)x_0\|$$
(🔍 we wrote $x_0=\lambda x_0+(1-\lambda)x_0$, a clever way of writing $x_0$)
$$=\|\lambda(x_1-x_0)+(1-\lambda)(x_2-x_0)\|\ \le\ \|\lambda(x_1-x_0)\|+\|(1-\lambda)(x_2-x_0)\|\quad(\text{triangle ineq.})$$
$$=\lambda\|x_1-x_0\|+(1-\lambda)\|x_2-x_0\|\ \le\ \lambda r+(1-\lambda)r=r .$$
So $z\in S$ ✓.

**Example 4 (HW): PSD cone** $S=\{X\in\mathbb{R}^{n\times n}:X\text{ symmetric},\ X\succeq0\}$.
🔍 *Solution:* If $X,Y\succeq0$, then for any $v$: $v^T(\lambda X+(1-\lambda)Y)v=\lambda\,v^TXv+(1-\lambda)\,v^TYv\ge0$ (both terms $\ge0$). A combination of symmetric matrices is symmetric. So it is convex ✓.

#### 3.2 Convex combinations and convex hull

A **convex combination** of $x_1,\dots,x_k$:
$$\sum_{i=1}^k\lambda_ix_i\quad\text{with}\quad\sum_{i=1}^k\lambda_i=1,\ \lambda_i\ge0 .$$
(A weighted average; for $k=2$ it is the segment.)

**Convex hull** of $S$: the **smallest convex set containing $S$**. Equivalently, all convex combinations of points of $S$. Picture: stretch a rubber band around the set. E.g. the convex hull of two separate squares is the hexagon-like shape enclosing both.

#### 3.3 Convex functions

$f:\Omega\to\mathbb{R}$, $\Omega\subseteq\mathbb{R}^n$, is **convex** if
1. $\operatorname{Dom}(f)=\Omega$ is a convex set, and
2. for all $x_1,x_2\in\Omega$ and $\lambda\in[0,1]$:
$$f\big(\lambda x_1+(1-\lambda)x_2\big)\ \le\ \lambda f(x_1)+(1-\lambda)f(x_2).$$

**Picture.** Draw the chord from $(x_1,f(x_1))$ to $(x_2,f(x_2))$; call it $g$. At $z=\lambda x_1+(1-\lambda)x_2$:
- LHS $=f(z)$, the height of the **curve**;
- RHS $=\lambda f(x_1)+(1-\lambda)f(x_2)$, the height of the **chord** above $z$.

E.g. $\lambda=0.2$: $z=0.2x_1+0.8x_2$ and the chord height is $0.2f(x_1)+0.8f(x_2)$. Convex means **the chord is always above (or on) the graph**: a "bowl" shape. A function with two dips is not convex (a chord across the hump goes *below* the graph).

#### 3.4 Convex optimization problems

$$(P)\qquad\min f(x)\ \text{ s.t. } x\in\Omega,\qquad f\text{ convex function},\ \Omega\text{ convex set}.$$

**Theorem 1. Every local optimum of (P) is a global optimum.**

*Proof by contradiction, every step.*
1. Let $x^*$ be a local optimum: $\exists r>0$ with $f(x^*)\le f(x)$ for all $x\in\Omega$, $\|x-x^*\|<r$.
2. Suppose it is **not** global: $\exists z\in\Omega$ with $f(z)<f(x^*)$.
3. Consider points on the segment from $x^*$ toward $z$:
$$y=(1-\lambda)x^*+\lambda z,\qquad\lambda\in(0,1].$$
$y\in\Omega$ because $\Omega$ is convex.
4. By convexity of $f$:
$$f(y)\le(1-\lambda)f(x^*)+\lambda f(z)<(1-\lambda)f(x^*)+\lambda f(x^*)=f(x^*).$$
(We replaced $f(z)$ by the larger $f(x^*)$; strict because $\lambda>0$.)
5. Now **choose $\lambda$ close to 0** so that $y$ is inside the ball $B(x^*,r)$. 🔍 $\|y-x^*\|=\lambda\|z-x^*\|$, so any $\lambda<\frac{r}{\|z-x^*\|}$ works.
6. Then $y$ is in the neighbourhood **and** $f(y)<f(x^*)$, contradicting step 1. ✗ So no such $z$ exists, and $x^*$ is global. ∎

**Theorem 2. The set $S^*$ of global optima of (P) is convex.**
(Picture: a flat-bottomed bowl, where the optima form an interval.)

*Proof.* Let $x_1,x_2\in S^*$, so $f(x_1)\le f(x)$ and $f(x_2)\le f(x)$ for all $x\in\Omega$. For $\lambda\in[0,1]$:
$$f(\lambda x_1+(1-\lambda)x_2)\le\lambda f(x_1)+(1-\lambda)f(x_2)\le\lambda f(x)+(1-\lambda)f(x)=f(x)\quad\forall x\in\Omega .$$
So $\lambda x_1+(1-\lambda)x_2$ is also a global optimum, i.e. in $S^*$. ∎
(Consequence: a convex problem has 0, 1 or infinitely many optimal solutions, never exactly two.)

#### 3.5 Operations that give convex sets

**(1) Intersection.** $S_1,S_2$ convex ⇒ $S_1\cap S_2$ convex.
🔍 *Proof:* $x,y\in S_1\cap S_2$ ⇒ the segment is in $S_1$ (convex) and in $S_2$ (convex) ⇒ in both.

*Application.* $H_1=\{x:a_1^Tx\le b_1\}$, $H_2=\{x:a_2^Tx\le b_2\}$ are half-spaces, so $H_1\cap H_2=\{x:a_i^Tx\le b_i,\ i=1,2\}$ is convex. Similarly any number of half-spaces. So the feasible set of a **Linear Program**
$$\text{LP:}\quad\min c^Tx\ \text{ s.t. } a_i^Tx\le b_i,\ i=1,\dots,K$$
is convex (a polyhedron), and $c^Tx$ is convex (linear). **LPs are convex problems.**

**Union is NOT convex in general.** $S_1\cup S_2$: two triangles touching at a point (bow-tie) is not convex. A segment between the two far corners leaves the set.

**(2) Minkowski sum.** $S=\{x:x=x_1+x_2,\ x_1\in S_1,\ x_2\in S_2\}$ is convex if $S_1,S_2$ are convex.
🔍 *Proof:* take $x=x_1+x_2$, $y=y_1+y_2$. Then $\lambda x+(1-\lambda)y=\underbrace{(\lambda x_1+(1-\lambda)y_1)}_{\in S_1}+\underbrace{(\lambda x_2+(1-\lambda)y_2)}_{\in S_2}\in S$.

Reference: Stephen Boyd & Lieven Vandenberghe, *Convex Optimization*.

---

### 4. Worked examples

1. Is $\{x\in\mathbb{R}^2:x_1^2+x_2^2\ge1\}$ (outside of the disc) convex? No: $(1,0)$ and $(-1,0)$ are in it, but their midpoint $(0,0)$ is not.
2. Is $\{x:x_1\ge0,\ x_2\ge0,\ x_1+x_2\le1\}$ convex? Yes: an intersection of three half-spaces (a triangle).
3. $f(x)=x^2$: $f(\frac{x+y}{2})\le\frac{f(x)+f(y)}2\iff\frac{(x+y)^2}{4}\le\frac{x^2+y^2}{2}\iff0\le\frac{(x-y)^2}{4}$ ✓ (the midpoint case of convexity).

---

### 5. Summary

✅ Convex set: contains the segment between any two of its points. Hyperplanes, half-spaces, balls and the PSD cone are convex.
✅ Convex hull = all convex combinations = smallest convex superset.
✅ Convex function: $f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$ (chord above graph).
✅ Convex problem: local optimum = global optimum; the optimal set is convex.
✅ Intersections and Minkowski sums of convex sets are convex; unions need not be. LPs are convex.

<!-- LECTURE END: convexitysep7 -->

---

<!-- LECTURE START: convexity-sep-9 -->
<a id="convexity-sep-9"></a>

## Lecture 19: Convexity (Part 2): Jensen's inequality, first- & second-order tests, examples, convexity-preserving operations, epigraph

*Source file:* `Convexity sep 9.pdf`

### 1. What this lecture is about

Checking convexity straight from the definition (chords) is tedious. This lecture gives practical tests and tools:

1. Concave functions; **Jensen's inequality** $f(E[X])\le E[f(X)]$
2. **First-order characterisation:** $f$ convex ⟺ $f(y)\ge f(x)+\nabla f(x)^T(y-x)$ (tangent lies below the graph), with a full proof of both directions
3. **Second-order characterisation:** $f$ convex ⟺ $\nabla^2f(x)\succeq0$ everywhere
4. Examples: $e^{ax}$, $\log x$, $\|x\|$, $\max(x_1,\dots,x_n)$, $c^Tx$
5. Operations preserving convexity: sum, max, $e^{g(x)}$, composition
6. **Epigraph:** $f$ convex ⟺ $\operatorname{epi}(f)$ is a convex set

---

### 2. Background you need first

- Convex sets and convex functions (Convexity Part 1).
- Directional derivative: $\lim_{\lambda\to0}\frac{f(y+\lambda d)-f(y)}{\lambda}=\nabla f(y)^Td$.
- Expected value of a discrete random variable: $E[X]=\sum_xx\,P(X=x)$.

---

### 3. Main content with full derivations

#### 3.1 Definitions recap

$f:\Omega\to\mathbb{R}$, $\Omega\subseteq\mathbb{R}^n$, is convex if (1) $\operatorname{Dom}(f)$ is convex and (2) $f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$ for all $x,y\in\operatorname{Dom}(f)$, $\lambda\in[0,1]$.

$f$ is **concave** if $-f$ is convex (chord lies *below* the graph: a "cap" shape).

#### 3.2 Jensen's inequality

For a random variable $X$ and a convex function $f$:
$$\boxed{f\big(E[X]\big)\le E\big[f(X)\big]}$$
Discrete case: $X=x_i$ with probability $p_i$, $i=1,\dots,K$, $\sum p_i=1$, $p_i\ge0$. Then $E[X]=\sum_ip_ix_i$ and $E[f(X)]=\sum_ip_if(x_i)$, so Jensen says
$$f\Big(\sum_{i=1}^Kp_ix_i\Big)\le\sum_{i=1}^Kp_if(x_i)\qquad\forall p:\ \textstyle\sum p_i=1,\ p_i\ge0 .$$
This is the convexity definition extended from 2 points to $K$ points (a convex combination).

🔍 *Proof by induction on $K$.* $K=2$ is the definition. Suppose it holds for $K-1$ points. If $p_K=1$ there is nothing to prove; otherwise write
$$\sum_{i=1}^Kp_ix_i=(1-p_K)\underbrace{\sum_{i=1}^{K-1}\frac{p_i}{1-p_K}x_i}_{=:u}+p_Kx_K .$$
The weights $\frac{p_i}{1-p_K}$ sum to 1. By the 2-point definition, $f(\cdot)\le(1-p_K)f(u)+p_Kf(x_K)$. By the induction hypothesis, $f(u)\le\sum_{i<K}\frac{p_i}{1-p_K}f(x_i)$. Combine to get $\sum_ip_if(x_i)$. ∎

#### 3.3 First-order characterisation

> **Theorem.** Let $f$ be differentiable. Then
> $$f\text{ convex}\iff f(y)\ge f(x)+\nabla f(x)^T(y-x)\qquad\forall x,y\in\operatorname{Dom}(f).$$

Picture: the right side is the **tangent line (plane)** at $x$. Convex ⟺ the graph lies **above every tangent**. A global under-estimator from purely local information!

**Proof (⇐): tangent condition ⇒ convex.**
1. Take any $x,y$, $\lambda\in[0,1]$, and set $z=\lambda x+(1-\lambda)y$ (also in the domain).
2. Apply the tangent inequality at the point $z$, once towards $y$ and once towards $x$:
$$f(y)\ge f(z)+\nabla f(z)^T(y-z)\quad\cdot(1-\lambda)$$
$$f(x)\ge f(z)+\nabla f(z)^T(x-z)\quad\cdot\lambda$$
3. Multiply as indicated (both multipliers $\ge0$, so the inequalities keep direction) and add:
$$(1-\lambda)f(y)+\lambda f(x)\ge f(z)+\nabla f(z)^T\big[(1-\lambda)y+\lambda x-z\big].$$
(🔍 the $f(z)$ terms add to $(1-\lambda+\lambda)f(z)=f(z)$; the gradient terms combine because $\nabla f(z)^T$ is common.)
4. But $(1-\lambda)y+\lambda x=z$, so the bracket is $z-z=0$. Hence
$$\lambda f(x)+(1-\lambda)f(y)\ge f(\lambda x+(1-\lambda)y),$$
which is convexity. ∎

**Proof (⇒): convex ⇒ tangent condition.**
1. Convexity: $f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$.
2. Rewrite the argument: $\lambda x+(1-\lambda)y=y+\lambda(x-y)$; and the right side: $\lambda f(x)+(1-\lambda)f(y)=f(y)+\lambda\big(f(x)-f(y)\big)$. So
$$f\big(y+\lambda(x-y)\big)\le f(y)+\lambda\big(f(x)-f(y)\big).$$
3. Subtract $f(y)$ and divide by $\lambda>0$:
$$\frac{f(y+\lambda(x-y))-f(y)}{\lambda}\le f(x)-f(y).$$
4. Let $\lambda\to0^+$. The left side becomes the directional derivative of $f$ at $y$ in direction $x-y$: $\nabla f(y)^T(x-y)$. So
$$\nabla f(y)^T(x-y)\le f(x)-f(y)\iff f(x)\ge f(y)+\nabla f(y)^T(x-y).$$
Swap the names $x\leftrightarrow y$ to match the theorem. ∎

🔍 **Big consequence:** if $\nabla f(x^*)=0$ and $f$ is convex, then $f(y)\ge f(x^*)+0$ for all $y$. So **for convex functions FONC is also sufficient for a global minimum**. (This is why the least-squares solution is global.)

#### 3.4 Second-order characterisation

> **Theorem.** Let $f$ be twice differentiable on an open convex domain. Then
> $$f\text{ convex}\iff\nabla^2f(x)\succeq0\quad\forall x\in\operatorname{Dom}(f).$$
> In 1-D: $f''(x)\ge0$ for all $x$ (the slope never decreases; the curve bends upward).

🔍 *Proof idea.* (⇐) Taylor with remainder: $f(y)=f(x)+\nabla f(x)^T(y-x)+\frac12(y-x)^T\nabla^2f(\xi)(y-x)$ for some $\xi$ on the segment. If the Hessian is PSD, the last term is $\ge0$, which gives the first-order condition. (⇒) If $v^T\nabla^2f(x)v<0$ for some $v$, then near $x$ along $v$ the function lies below its tangent, contradicting the first-order condition.

#### 3.5 Examples

1. $f(x)=e^{ax}$ is convex **for every real $a$**: $f'=ae^{ax}$, $f''=a^2e^{ax}\ge0$ for all $x$.
2. $f(x)=\log x$ ($x>0$): $f'=\frac1x$, $f''=-\frac1{x^2}<0$ ⇒ **concave**.
3. $f(x)=\|x\|$ (any norm) is convex. 🔍 By the triangle inequality and homogeneity: $\|\lambda x+(1-\lambda)y\|\le\lambda\|x\|+(1-\lambda)\|y\|$.
4. $f(x)=\max(x_1,\dots,x_n)$ is convex:
$$f(\lambda x+(1-\lambda)y)=\max_i\big[\lambda x_i+(1-\lambda)y_i\big]\le\lambda\max_ix_i+(1-\lambda)\max_iy_i .$$
🔍 For each $i$: $\lambda x_i\le\lambda\max_jx_j$ and $(1-\lambda)y_i\le(1-\lambda)\max_jy_j$; add, then take the max over $i$ of the left side.
5. $f(x)=c^Tx=\sum_ic_ix_i$: $c^T(\lambda x+(1-\lambda)y)=\lambda c^Tx+(1-\lambda)c^Ty$, **equality**. So linear functions are both **convex and concave**.

#### 3.6 Operations preserving convexity of functions

1. **Sum:** $f=f_1+f_2$ with $f_1,f_2$ convex (same domain) is convex. (Add the two defining inequalities.) Also $\alpha f$ for $\alpha\ge0$.
2. **Pointwise max:** $f(x)=\max(f_1(x),f_2(x))$ is convex. (The upper envelope of two bowls, as drawn.) 🔍 $f(\lambda x+(1-\lambda)y)=f_j(\lambda x+(1-\lambda)y)$ for some $j$; that is $\le\lambda f_j(x)+(1-\lambda)f_j(y)\le\lambda f(x)+(1-\lambda)f(y)$.
3. **(HW)** $f(x)=e^{g(x)}$ with $g$ convex (1-D, twice differentiable):
$$f'(x)=e^{g(x)}g'(x),\qquad f''(x)=e^{g(x)}g''(x)+\big(g'(x)\big)^2e^{g(x)}\ \ge0\ \forall x,$$
since $e^{g}>0$, $g''\ge0$, and a square is $\ge0$. ✓ Convex.
4. **(HW)** composition $f(x)=g(h(x))$. 🔍 Rule (1-D): $f''=g''(h)\,(h')^2+g'(h)\,h''$. So $f$ is convex if **$g$ is convex and non-decreasing and $h$ is convex** (both terms $\ge0$). Also if $g$ is convex and $h$ is affine ($h''=0$). Example 3 is the case $g=\exp$ (convex, increasing).

#### 3.7 Epigraph

For $f:\mathbb{R}^n\to\mathbb{R}$, the **epigraph** is the region on or above the graph:
$$\operatorname{epi}(f)=\Big\{\begin{pmatrix}x\\t\end{pmatrix}\in\mathbb{R}^{n+1}:\ x\in\mathbb{R}^n,\ t\in\mathbb{R},\ f(x)\le t\Big\}.$$
> $$\operatorname{epi}(f)\text{ is a convex set}\iff f\text{ is a convex function.}$$

(A wiggly function has a non-convex region above it; a bowl has a convex one.)

🔍 *Proof (⇒ direction for functions):* if $(x,t),(y,s)\in\operatorname{epi}f$ then $f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)\le\lambda t+(1-\lambda)s$, so the combination is in the epigraph. (⇐) Apply set-convexity to the points $(x,f(x))$ and $(y,f(y))$.

---

### 4. Worked examples

1. Is $f(x_1,x_2)=x_1^2+x_1x_2+x_2^2$ convex? Hessian $\begin{pmatrix}2&1\\1&2\end{pmatrix}$: minors $2>0$, $3>0$ ⇒ PD ⇒ convex.
2. $f(x)=\max(x^2,\ 2x+1)$: max of two convex functions ⇒ convex.
3. Jensen with $f=x^2$: $(E[X])^2\le E[X^2]$, i.e. $\operatorname{Var}(X)=E[X^2]-(E[X])^2\ge0$.

---

### 5. Summary

✅ Concave ⟺ $-f$ convex. Jensen: $f(\sum p_ix_i)\le\sum p_if(x_i)$.
✅ 1st order: convex ⟺ $f(y)\ge f(x)+\nabla f(x)^T(y-x)$ ⇒ $\nabla f(x^*)=0$ gives a global min.
✅ 2nd order: convex ⟺ $\nabla^2f\succeq0$ everywhere ($f''\ge0$ in 1-D).
✅ Convex: $e^{ax}$, norms, max, linear. Concave: $\log x$.
✅ Sums, maxima, $e^{\text{convex}}$ are convex. $f$ is convex ⟺ its epigraph is a convex set.

<!-- LECTURE END: convexity-sep-9 -->

---

<!-- LECTURE START: sep11 -->
<a id="sep11"></a>

## Lecture 20: Convex Hulls, Carathéodory's Theorem, and "every problem is convex"

*Source file:* `Sep11.pdf`

### 1. What this lecture is about

1. A precise formula for the convex hull $\operatorname{conv}(S)$
2. **Carathéodory's theorem:** in $\mathbb{R}^n$, every point of $\operatorname{conv}(S)$ is a convex combination of **at most $n+1$** points of $S$ (full proof)
3. Minimising a linear function over $S$ gives the same value as over $\operatorname{conv}(S)$
4. Using that (plus the epigraph trick), **every** optimization problem can be rewritten as a convex one. So "convex" does **not** automatically mean "easy".

---

### 2. Background you need first

- Convex combination $\sum\lambda_ix_i$, $\lambda_i\ge0$, $\sum\lambda_i=1$; convex hull (Convexity Part 1).
- Rank–Nullity: an $(n+1)\times m$ matrix with $m>n+1$ has rank $\le n+1<m$, so nullity $=m-\text{rank}>0$, i.e. a non-zero solution $y$ of $Ay=0$ exists (LA block).
- Epigraph reformulation $\min f\iff\min t$ s.t. $f(x)\le t$ (1-D optimization & Convexity Part 2).

---

### 3. Main content with full derivations

#### 3.1 Convex hull formula

$$\operatorname{conv}(S)=\Big\{x\ \Big|\ x=\sum_{i=1}^m\lambda_ix_i,\ x_i\in S,\ \lambda_i\ge0,\ \sum_{i=1}^m\lambda_i=1,\ \text{for some }m\Big\}.$$
Note "for some $m$": a priori you might need arbitrarily many points. Carathéodory says you never need more than $n+1$.

#### 3.2 Carathéodory's theorem

> Let $S\subseteq\mathbb{R}^n$. Every point in $\operatorname{conv}(S)$ can be written as a convex combination of **$n+1$ (or fewer) points** of $S$.

(In the plane, $n=2$: every point of the hull lies in a triangle with corners in $S$.)

**Proof, all steps.**
1. Take $x\in\operatorname{conv}(S)$: $x=\sum_{i=1}^m\lambda_ix_i$, $x_i\in S$, $\lambda_i\ge0$, $\sum\lambda_i=1$. Take this representation with the **smallest** possible $m$ (so all $\lambda_i>0$).
2. If $m\le n+1$: done.
3. Assume $m>n+1$. We show we can remove one point, which contradicts minimality (or, as an algorithm: repeat until $m\le n+1$).
4. **Find a "balanced" set of weights $y$.** We want $y_1,\dots,y_m\in\mathbb{R}$, not all zero, with
$$\sum_{i=1}^my_i=0\quad\text{and}\quad\sum_{i=1}^my_ix_i=0 .$$
In matrix form, stack a row of 1's on top of the points:
$$A=\begin{pmatrix}1&1&\cdots&1\\x_1&x_2&\cdots&x_m\end{pmatrix}\in\mathbb{R}^{(n+1)\times m},\qquad Ay=0 .$$
(The first row gives $\sum y_i=0$; the remaining $n$ rows give $\sum y_ix_i=0$.)
5. $\operatorname{rank}(A)\le n+1<m$, so by rank–nullity, $\text{nullity}=m-\operatorname{rank}(A)>0$: **a non-trivial null-space vector $y$ exists.**
6. Since $\sum y_i=0$ and $y\neq0$, **at least one $y_i>0$** (if all were $\le0$ with sum 0, all would be 0).
7. **Slide the weights.** For any $\varepsilon$:
$$\sum_{i=1}^m(\lambda_i-\varepsilon y_i)x_i=\underbrace{\sum_i\lambda_ix_i}_{=x}-\varepsilon\underbrace{\sum_iy_ix_i}_{=0}=x,\qquad\sum_i(\lambda_i-\varepsilon y_i)=1-\varepsilon\cdot0=1 .$$
So $\lambda_i-\varepsilon y_i$ are new weights giving the same $x$ and summing to 1. We just need them $\ge0$, and we want one of them to become exactly 0.
8. **Choose $\varepsilon$.** We need $\lambda_i-\varepsilon y_i\ge0$ for all $i$.
   - If $y_i\le0$: automatically true for $\varepsilon\ge0$.
   - If $y_i>0$: need $\varepsilon\le\frac{\lambda_i}{y_i}$.
   So take the largest allowed value:
$$\varepsilon=\min_{i:\ y_i>0}\frac{\lambda_i}{y_i}\ (>0).$$
9. For the index $i$ achieving this minimum, $\lambda_i-\varepsilon y_i=0$. **At least one weight becomes zero**, so that point can be dropped. Now $x$ is a convex combination of $m-1$ points of $S$, contradicting minimality of $m$. Hence $m\le n+1$. ∎

#### 3.3 Linear objective: optimizing over $S$ vs over $\operatorname{conv}(S)$

$$P1=\min_{x\in S}\sum_{i=1}^nc_ix_i=\min_{x\in S}c^Tx,\qquad P2=\min_{x\in\operatorname{conv}(S)}c^Tx .$$
$S$ may be non-convex. **Claim: $P1=P2$.**

- **$P2\le P1$:** $S\subseteq\operatorname{conv}(S)$, and minimising over a bigger set can only give a smaller (or equal) value.
- **$P2\ge P1$:** let $z=\arg\min_{x\in\operatorname{conv}(S)}c^Tx$. By Carathéodory, $z=\sum_{i=1}^{n+1}\lambda_iq_i$ with $q_i\in S$, $\lambda_i\ge0$, $\sum\lambda_i=1$. Then
$$P2=c^Tz=c^T\Big(\sum_i\lambda_iq_i\Big)=\sum_i\lambda_i\,(c^Tq_i)\ \ge\ \sum_i\lambda_i\,P1=P1,$$
since each $q_i\in S$ gives $c^Tq_i\ge\min_{x\in S}c^Tx=P1$, and $\sum\lambda_i=1$.

So $P1=P2$ ∎. **A linear function over a non-convex set has the same minimum as over its convex hull.**

#### 3.4 Every optimization problem can be made convex

Start with a general (possibly non-convex) problem
$$Q1:\quad\min f(x)\ \text{ s.t. } x\in S\subseteq\mathbb{R}^n .$$
**Step 1: epigraph reformulation.**
$$Q2:\quad\min_{x,t}\ t\ \text{ s.t. } f(x)\le t,\ x\in S .$$
The objective is now **linear** in the variable $\binom{t}{x}$: $t=(1,0,\dots,0)\binom tx$.
**Step 2: convexify the feasible set.** Let $\mathcal X=\{\binom tx: f(x)\le t,\ x\in S\}$. By §3.3 (linear objective!),
$$Q3:\quad\min t\ \text{ s.t. }\binom tx\in\operatorname{conv}(\mathcal X)$$
has the **same optimal value**. And $Q3$ is a convex problem (linear objective, convex feasible set).

(Sketch: the shaded region above the wiggly graph; the lowest horizontal level $t^*$ is the same for the region and for its convex hull.)

✅ **Takeaway:** every optimization problem can be reformulated as a convex optimization problem. So **convex optimization is not always easy!** The difficulty moves into *describing* $\operatorname{conv}(\mathcal X)$, which can be extremely hard. Convex problems are easy only when the convex set/function is given in a tractable form.

---

### 4. Worked example (Carathéodory in 1-D)

$n=1$: any point of $\operatorname{conv}(S)$ is a combination of at most 2 points. Take $S=\{0,1,3\}$ and $x=1.6=0.2\cdot0+0.4\cdot1+0.4\cdot3$ (3 points, weights $\lambda=(0.2,0.4,0.4)$).
- Null-space $y$: need $y_1+y_2+y_3=0$ and $0\cdot y_1+1\cdot y_2+3\cdot y_3=0$. One solution: $y=(2,-3,1)$.
- Positive $y$'s: $y_1=2$ (ratio $\lambda_1/y_1=0.1$) and $y_3=1$ (ratio $0.4$). So $\varepsilon=0.1$.
- New weights $\lambda-0.1y=(0.2-0.2,\ 0.4+0.3,\ 0.4-0.1)=(0,\ 0.7,\ 0.3)$.
- Check: $0.7\cdot1+0.3\cdot3=1.6=x$ ✓, and $0.7+0.3=1$ ✓. Only 2 points are used now.

---

### 5. Summary

✅ $\operatorname{conv}(S)$ = all finite convex combinations of points of $S$.
✅ Carathéodory: $n+1$ points suffice. Proof: a null-space vector of $\binom{1\cdots1}{x_1\cdots x_m}$ lets us slide weights until one hits 0.
✅ $\min_Sc^Tx=\min_{\operatorname{conv}S}c^Tx$.
✅ Epigraph + convex hull ⇒ any problem becomes convex. Convex ≠ easy.

<!-- LECTURE END: sep11 -->

---

<!-- LECTURE START: sephyperplane -->
<a id="sephyperplane"></a>

## Lecture 21: Separating Hyperplane Theorem (partial notes)

*Source file:* `SepHyperplane.pdf`

> 📌 **This file is the first part of `FL.pdf`.** The full explanation is in the block **"Separating Hyperplane Theorem and Farkas' Lemma"**. That block has the complete separating-hyperplane proof, Farkas' lemma with both directions proved, and the variants. This block is a checklist.

### Checklist

1. **Theorem:** for disjoint convex sets $C,D\subseteq\mathbb{R}^n$ there exist $a,b$ with $a^Tx\le b$ on $C$ and $a^Tx\ge b$ on $D$.
2. **Proof for closed, bounded sets:** take the closest pair $c\in C$, $d\in D$ (exists by Weierstrass). Set $a=d-c$ and $b=\frac{\|d\|^2-\|c\|^2}{2}$, so the wall passes through the midpoint $\frac{c+d}{2}$.
3. If some $\tilde c\in C$ had $a^T\tilde c>b$, then $\nabla f(c)^T(\tilde c-c)<-\|c-d\|^2<0$ for $f(x)=\|x-d\|^2$. So $\tilde c-c$ would be a descent direction inside $C$ (by convexity), contradicting that $c$ is the closest point.
4. **Point vs closed convex set:** if $x_0\notin C$, there are $a,b$ with $a^Tx>b$ on $C$ and $a^Tx_0\le b$.
5. **Farkas' lemma (statement):** exactly one of $\{x:Ax=b,\ x\ge0\}$ and $\{y:A^Ty\le0,\ b^Ty>0\}$ is non-empty. Geometrically: either $b$ lies in the cone spanned by the columns $a_i$, or a hyperplane separates them.

<!-- LECTURE END: sephyperplane -->

---

<!-- LECTURE START: oct7 -->
<a id="oct7"></a>

## Lecture 22: Linear Programming: standard form and applications (max flow, transportation, diet)

*Source file:* `Oct7.pdf`

### 1. What this lecture is about

We now specialise to the most important constrained problem class, **Linear Programs (LPs)**: linear objective, linear constraints. This lecture:
1. Gives the LP form and compares it with general constrained optimization
2. Models three classic applications as LPs: **max flow**, **transportation**, **diet**

---

### 2. Background you need first

- Vectors and matrix–vector products; $x\ge0$ means every component $x_i\ge0$.
- Graphs: nodes and directed edges $(i,j)$.
- Linear function: $c^Tx=\sum c_ix_i$ (no squares, no products of variables).

---

### 3. Main content

#### 3.1 The LP form

$$\min_{x\in\mathbb{R}^n}\ c^Tx\qquad\text{s.t.}\quad Ax=b,\ \ x\ge0,\qquad A\in\mathbb{R}^{m\times n}.$$
- Objective **linear** in $x$; constraint functions **linear** in $x$.

Compare the general constrained problem
$$\min f(x)\quad\text{s.t.}\quad g_i(x)\le0,\ \ h_j(x)=0 .$$
An LP is the special case where $f,g_i,h_j$ are all linear (affine).

#### 3.2 Application 1: Maximum flow

Send the maximum quantity of product from $s$ to $t$ in a network, subject to link capacities $c_{ij}$ and flow conservation.
- Variables: $x_{ij}$ = flow on edge $(i,j)$; $v$ = total flow value.
- $$\max\ v$$
$$\text{s.t.}\quad c_{ij}\ge x_{ij}\ge0\quad\forall(i,j)\in E,$$
$$\sum_{k:(k,i)\in E}x_{ki}-\sum_{k:(i,k)\in E}x_{ik}=\begin{cases}0&i\neq s,t\\ v&i=t\\-v&i=s\end{cases}$$
🔍 Read the second line as "flow **into** $i$ minus flow **out of** $i$". Intermediate nodes store nothing ($=0$). The sink receives net $v$, the source sends out net $v$ (in − out $=-v$).

#### 3.3 Application 2: Transportation problem

$m$ sources, $n$ destinations. Source $i$ has **supply** $a_i$; destination $j$ has **demand** $b_j$; sending one unit from $i$ to $j$ costs $c_{ij}$. Satisfy all demands at least cost.
- Variables: $x_{ij}$ = number of items sent from source $i$ to destination $j$.
$$\min\sum_{i=1}^m\sum_{j=1}^nc_{ij}x_{ij}$$
$$\text{s.t.}\quad\sum_{i=1}^mx_{ij}\ge b_j\ \ \forall j=1..n\quad(\text{each destination gets its demand}),$$
$$\sum_{j=1}^nx_{ij}\le a_i\ \ \forall i=1..m\quad(\text{a source can't ship more than it has}),$$
$$x_{ij}\ge0\ \ \forall i,j .$$
🔍 For this to be feasible we need total supply $\ge$ total demand: $\sum_ia_i\ge\sum_jb_j$.

#### 3.4 Application 3: Diet problem

$n$ foods (rice, dal, chapathi, …), $m$ nutrients (carbs, protein, vitamins, fat, …).
- One unit of food $i$ contains $a_{ij}$ units of nutrient $j$.
- Body requirement of nutrient $j$: $b_j$.
- Cost of food $i$ per unit: $c_i$.

Obtain the required nutrition at least overall cost. Variables: $x_i$ = quantity of food $i$ to buy.
$$\min\sum_{i=1}^nc_ix_i\quad\text{s.t.}\quad\sum_{i=1}^nx_ia_{ij}\ge b_j\ \ \forall j=1..m,\qquad x_i\ge0\ \ \forall i .$$
🔍 In matrix form with $M_{ji}=a_{ij}$: $\min c^Tx$ s.t. $Mx\ge b$, $x\ge0$. That is an LP with inequality constraints, convertible to $Ax=b$ form with surplus variables (next lectures).

---

### 4. Worked example (tiny diet)

Two foods: rice (₹2/unit, 4 carbs, 1 protein) and dal (₹5/unit, 1 carb, 4 protein). Need ≥8 carbs, ≥8 protein.
$$\min2x_1+5x_2\ \text{ s.t. }4x_1+x_2\ge8,\ x_1+4x_2\ge8,\ x\ge0 .$$
Corners of the feasible region: $(8,0)$ cost 16; $(0,8)$ cost 40; intersection $4x_1+x_2=8$, $x_1+4x_2=8$ gives $x_1=x_2=1.6$, cost $3.2+8=11.2$. Best: $(1.6,1.6)$, cost ₹11.2. (The optimum is at a corner. The next lectures explain why.)

---

### 5. Summary

✅ LP: $\min c^Tx$ s.t. $Ax=b,\ x\ge0$ (everything linear).
✅ Max flow, transportation and diet are all LPs: choose variables, write a linear cost and linear constraints.

<!-- LECTURE END: oct7 -->

---

<!-- LECTURE START: fl -->
<a id="fl"></a>

## Lecture 23: Separating Hyperplane Theorem and Farkas' Lemma (theorem of alternatives)

*Source file:* `FL.pdf`

### 1. What this lecture is about

Two foundational results that later give us **LP duality** and **KKT conditions**:

1. **Separating hyperplane theorem:** two disjoint convex sets can be separated by a flat "wall" $a^Tx=b$. Full proof for closed, bounded sets, using Weierstrass + descent directions.
2. Point-vs-closed-convex-set separation.
3. **Farkas' lemma:** for $Ax=b,\ x\ge0$, *exactly one* of two systems has a solution. Full proof via cones and separation.
4. Variants: converting $Ax\ge b$ and free variables into the standard form.

---

### 2. Background you need first

- Convex sets (Convexity Part 1); closed/bounded/compact; **Weierstrass:** continuous $f$ on a compact set attains its min (Topology).
- Descent direction: $\nabla f(c)^Td<0$ ⇒ $f$ decreases a little along $d$ (Descent Methods block).
- $\|x-d\|^2=(x-d)^T(x-d)$, gradient $2(x-d)$.
- Hyperplane $\{x:a^Tx=b\}$ splits space into the half-spaces $a^Tx\le b$ and $a^Tx\ge b$.

---

### 3. Main content with full derivations

#### 3.1 Statement

> **Separating Hyperplane Theorem.** Let $C,D\subseteq\mathbb{R}^n$ be convex with $C\cap D=\emptyset$. Then there exist $a\in\mathbb{R}^n$ ($a\neq0$) and $b\in\mathbb{R}$ with
> $$a^Tx\le b\ \ \forall x\in C\qquad\text{and}\qquad a^Tx\ge b\ \ \forall x\in D .$$

Picture: two blobs, a wall between them; $C$ on the side $a^Tx\le b$, $D$ on the side $a^Tx\ge b$.

We prove it when $C,D$ are additionally **closed and bounded**.

#### 3.2 Construction of $a$ and $b$

1. Consider $\inf\{\|x-y\|^2:\ x\in C,\ y\in D\}$, the squared distance between the sets.
2. The objective is continuous in $(x,y)$ and we optimise over a compact set ($C\times D$, closed and bounded). By **Weierstrass**, the infimum is attained at some $x^*=c\in C$, $y^*=d\in D$: the two **closest points**. Since the sets are disjoint, $c\neq d$.
3. Define
$$a=d-c,\qquad b=\frac{\|d\|^2-\|c\|^2}{2}.$$
🔍 *Where this $b$ comes from:* the wall should pass through the midpoint $\frac{c+d}{2}$ and be perpendicular to $d-c$. So $b=a^T\frac{c+d}{2}=\frac{(d-c)^T(d+c)}{2}=\frac{d^Td-c^Tc}{2}$ (the cross terms $d^Tc-c^Td$ cancel).

#### 3.3 Proof of the claim $a^Tx\le b$ for all $x\in C$

(The $D$ side is symmetric.)

1. Suppose not: $\exists\tilde c\in C$ with $a^T\tilde c>b$, i.e.
$$(d-c)^T\tilde c>\frac{\|d\|^2-\|c\|^2}{2}. \tag{*}$$
2. Since $c$ is the closest point of $C$ to $d$: $\inf_{x\in C}\|x-d\|^2=\|c-d\|^2$, i.e. $c$ minimises $f(x)=\|x-d\|^2=(x-d)^T(x-d)$ over $C$.
3. Points on the segment $c+\alpha(\tilde c-c)=(1-\alpha)c+\alpha\tilde c$, $\alpha\in[0,1]$, are in $C$ (**$C$ is convex**).
4. $\nabla f(x)=2(x-d)$, so $\nabla f(c)=2(c-d)$. Is $\nabla f(c)^T(\tilde c-c)<0$? (If yes, $\tilde c-c$ is a descent direction and we can get closer to $d$ while staying in $C$, a contradiction.)
5. Compute:
$$\nabla f(c)^T(\tilde c-c)=2(c-d)^T\tilde c-2c^Tc+2d^Tc .$$
By (*), $2(c-d)^T\tilde c=-2(d-c)^T\tilde c<-(\|d\|^2-\|c\|^2)=\|c\|^2-\|d\|^2$. So
$$\nabla f(c)^T(\tilde c-c)<\|c\|^2-\|d\|^2-2\|c\|^2+2d^Tc=-\|c\|^2-\|d\|^2+2d^Tc=-\big[(c-d)^T(c-d)\big]=-\|c-d\|^2<0 .$$
(🔍 $\|c-d\|^2=\|c\|^2-2c^Td+\|d\|^2$, so $-\|c\|^2-\|d\|^2+2d^Tc=-\|c-d\|^2$.)
6. So $\tilde c-c$ **is a descent direction for $f$ at $c$**: for small $\alpha>0$, $f(c+\alpha(\tilde c-c))<f(c)$, and this point is in $C$. This contradicts $c$ being the closest point. ✗ Hence $a^Tx\le b$ for all $x\in C$. ∎

#### 3.4 Separating a point from a closed convex set

> If $C\subseteq\mathbb{R}^n$ is closed and convex and $x_0\notin C$, then $\exists a\in\mathbb{R}^n,\ b\in\mathbb{R}$ with
> $$a^Tx>b\ \ \forall x\in C\qquad\text{and}\qquad a^Tx_0\le b .$$

(Apply the same closest-point argument with $D=\{x_0\}$; a single point is compact, and for closed $C$ the closest point exists. The separation can be made strict.)

#### 3.5 Farkas' lemma

Let $A\in\mathbb{R}^{m\times n}$ with columns $a_1,\dots,a_n\in\mathbb{R}^m$, and $b\in\mathbb{R}^m$.

> **Farkas' Lemma.** Exactly one of the following two systems is feasible (non-empty):
> 1. $\{x\in\mathbb{R}^n:\ Ax=b,\ x\ge0\}$ (i.e. $x_1,\dots,x_n\ge0$)
> 2. $\{y\in\mathbb{R}^m:\ a_i^Ty\le0\ \forall i=1,\dots,n,\ \ b^Ty>0\}=\{y:\ A^Ty\le0,\ b^Ty>0\}$

This is a **theorem of alternatives**.

**Geometric meaning.** $Ax=\sum_ix_ia_i$ with $x_i\ge0$ = a non-negative combination of the columns: the **cone** spanned by $a_1,\dots,a_n$ (the shaded wedge). System 1 says "$b$ is in this cone". System 2 says "there is a $y$ making an obtuse angle ($\le90°$ … i.e. $a_i^Ty\le0$) with every $a_i$ but an acute angle with $b$". That $y$ is the normal of a hyperplane through the origin separating $b$ from the cone. Either $b$ is inside the cone, or a separating hyperplane exists, never both. (If the $a_i$'s surround the origin, the cone is everything and system 1 is always feasible.)

**Proof.**

*(a) Both cannot hold.* Suppose (1) is feasible: $\exists x$, $Ax=b$, $x_i\ge0$. If (2) were also feasible with some $y$ ($A^Ty\le0$, $b^Ty>0$), then
$$b^Ty=(Ax)^Ty=x^TA^Ty=\sum_{i=1}^n\underbrace{x_i}_{\ge0}\underbrace{(A^Ty)_i}_{\le0}\le0,$$
contradicting $b^Ty>0$. ✗ So (2) is infeasible.

*(b) If (1) fails, (2) holds.* Suppose $\nexists x$ with $Ax=b,\ x\ge0$.
1. Let $K_1=\{v:\ v=\sum_{i=1}^na_ix_i,\ x_i\ge0\}$. Then "(1) infeasible" means exactly **$b\notin K_1$**.
2. **Cone:** a set $K$ with $v\in K\Rightarrow\alpha v\in K$ for all $\alpha\ge0$. $K_1$ is a cone, it is **convex** (HW: a non-negative combination of non-negative combinations is again one), and it is **closed** (HW; true for finitely generated cones).
3. Separate the point $b$ from the closed convex set $K_1$ (§3.4): $\exists y,r$ with
$$z^Ty\le r\ \ \forall z\in K_1,\qquad b^Ty>r .$$
4. **Pin down $r$.**
   - $0\in K_1$ (all $x_i=0$), so $0=0^Ty\le r$: **$r\ge0$**.
   - Can there be $z\in K_1$ with $z^Ty>0$? Then $\alpha z\in K_1$ (cone) and $(\alpha z)^Ty=\alpha z^Ty\to\infty$ as $\alpha\to\infty$, eventually exceeding $r$ (the notes test "$r=1500$?" and see that any positive bound gets beaten). ✗ So $z^Ty\le0$ for all $z\in K_1$: we may take $r=0$.
   - And $b^Ty>r\ge0$, so $b^Ty>0$.
5. Each column $a_i\in K_1$ (take $x_i=1$, others 0), so $a_i^Ty\le0$ for all $i$.
6. Thus $y$ satisfies $A^Ty\le0$, $b^Ty>0$: **(2) is feasible**. ∎

#### 3.6 Variants (reducing other systems to the Farkas form)

**V1: inequalities.** $\{x\in\mathbb{R}^n:\ Ax\ge b,\ x\ge0\}$. Introduce **surplus variables** $z=Ax-b\ge0$:
$$\iff\{x\in\mathbb{R}^n,z\in\mathbb{R}^m:\ Ax-z=b,\ x\ge0,\ z\ge0\},\quad\text{i.e.}\quad\begin{bmatrix}A&-I\end{bmatrix}\begin{bmatrix}x\\z\end{bmatrix}=b .$$
Now it has the form "equalities + non-negativity", and Farkas applies with matrix $[A\ -I]$.
🔍 Resulting alternative: either $Ax\ge b,x\ge0$ is feasible, or $\exists y$ with $A^Ty\le0$, $-y\le0$ (i.e. $y\ge0$), $b^Ty>0$.

**V2: free variables.** $\{x\in\mathbb{R}^n:\ Ax=b\}$ ($x$ unrestricted in sign). Write $x=x^+-x^-$ with $x^+\ge0$, $x^-\ge0$ (any real number is a difference of two non-negative ones):
$$\iff\{x^+,x^-\in\mathbb{R}^n:\ A(x^+-x^-)=b,\ x^+\ge0,\ x^-\ge0\}.$$
🔍 Farkas with matrix $[A\ -A]$: either $Ax=b$ is solvable, or $\exists y$ with $A^Ty\le0$ and $-A^Ty\le0$ (i.e. $A^Ty=0$) and $b^Ty>0$. (This is the familiar "$b\notin\mathcal R(A)$ ⟺ some $y\perp$ all columns has $b^Ty\neq0$".)

---

### 4. Worked example

$A=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ (columns $a_1=(1,0)$, $a_2=(1,1)$), $b=(-1,0)$.
- System 1: $x_1+x_2=-1$, $x_2=0$ ⇒ $x_1=-1<0$. Infeasible.
- System 2 must then be feasible: try $y=(-1,0)$: $a_1^Ty=-1\le0$, $a_2^Ty=-1\le0$, $b^Ty=1>0$ ✓.

---

### 5. Summary

✅ Disjoint convex $C,D$ ⇒ a hyperplane $a^Tx=b$ separates them. Proof: closest points $c,d$; $a=d-c$; $b=\frac{\|d\|^2-\|c\|^2}{2}$; a violating $\tilde c$ would give a descent direction for $\|x-d\|^2$.
✅ Point vs closed convex set: strict separation.
✅ Farkas: exactly one of $\{Ax=b,x\ge0\}$ and $\{A^Ty\le0,b^Ty>0\}$ is feasible.
✅ Use surplus variables ($Ax\ge b$) and $x=x^+-x^-$ (free variables) to reach the standard form.

<!-- LECTURE END: fl -->

---

<!-- LECTURE START: lp-bfs -->
<a id="lp-bfs"></a>

## Lecture 24: LP Solutions: standard form, slack/surplus/free variables, basic (feasible) solutions, Fundamental Theorem

*Source file:* `Lp-bfs.pdf`

### 1. What this lecture is about

1. Three possible outcomes of an LP: feasible & bounded, feasible & unbounded, infeasible
2. **Standard form** $\min c^Tx$ s.t. $Ax=b$, $x\ge0$, and how to convert any LP into it (slack, surplus, free variables)
3. **Basic solutions** and **basic feasible solutions (b.f.s.)**
4. **Fundamental Theorem of LP:** if an optimum exists, one exists at a b.f.s. (so we only need to search finitely many points)

---

### 2. Background you need first

- Rank, linear independence, invertible matrices (LA block).
- Block matrix multiplication: $[B\ N]\binom{x_B}{x_N}=Bx_B+Nx_N$.

---

### 3. Main content with full derivations

#### 3.1 Possible outcomes

An LP is exactly one of:
- **Feasible & bounded:** has an optimal value (and an optimal solution).
- **Feasible & unbounded:** objective can go to $-\infty$ (for min).
- **Infeasible:** no point satisfies the constraints.

#### 3.2 Standard form

$$\min c^Tx\quad\text{s.t.}\quad Ax=b,\ x\ge0,\qquad A\in\mathbb{R}^{m\times n},\ x\in\mathbb{R}^n,\ b\in\mathbb{R}^m .$$

**Conversions:**
- **Slack variables** for $\le$: $Ax\le b,\ x\ge0$ ⟺ $Ax+z=b,\ x\ge0,\ z\ge0$, $z\in\mathbb{R}^m$. ($z_i$ = unused amount, the "slack".)
- **Surplus variables** for $\ge$: $Ax\ge b,\ x\ge0$ ⟺ $Ax-z=b,\ x\ge0,\ z\ge0$.
- **Free variables** (no sign constraint): $x=x^+-x^-$ with $x^+\ge0$, $x^-\ge0$.

🔍 *Example:* $\min x_1-x_2$ s.t. $x_1+x_2\le4$, $x_1\ge0$, $x_2$ free. Becomes: $\min x_1-x_2^++x_2^-$ s.t. $x_1+x_2^+-x_2^-+s=4$, all four variables $\ge0$.

#### 3.3 Basic solutions

Assume $Ax=b$ with $A\in\mathbb{R}^{m\times n}$, $m<n$ (more variables than equations), and **$A$ full rank**: $\operatorname{rank}(A)=m$.

Pick $m$ **linearly independent columns** of $A$ and (after reordering) write
$$A=[\,B\mid N\,],\qquad B\in\mathbb{R}^{m\times m}\ (\text{invertible}),\quad N\in\mathbb{R}^{m\times(n-m)},\qquad x=\begin{pmatrix}x_B\\x_N\end{pmatrix}.$$
$x_B$ = **basic variables**, $x_N$ = **non-basic variables**. Then
$$Ax=b\iff Bx_B+Nx_N=b .$$
Set the non-basic variables to zero, $x_N=0$: then $Bx_B=b$, so $x_B=B^{-1}b$. The resulting
$$x=(x_1,x_2,\dots,x_m,0,0,\dots,0)^T$$
is a **basic solution**. The basic columns need not be the first $m$: e.g. $x=(x_1,0,\dots,0,x_{k+1},x_{k+2},\dots,0)$, where the non-zero positions correspond to the chosen columns. The order inside $x_B$ ($x_1,x_2,x_3$ vs $x_2,x_1,x_3$) doesn't matter as long as it matches the column order of $B$.

**Basic feasible solution (b.f.s.):** a basic solution that additionally has $x_B\ge0$ (all $x_i\ge0$).

🔍 *How many?* At most $\binom{n}{m}$ ways to choose the basis, so **finitely many** basic solutions.

#### 3.4 Fundamental Theorem of Linear Programming

> For an LP in standard form (with $\operatorname{rank}A=m$):
> 1. If it has a **feasible** solution, it has a **basic feasible** solution.
> 2. If it has an **optimal** feasible solution, it has an **optimal basic feasible** solution.

So to solve an LP it suffices to search among the finitely many b.f.s. (This is what the simplex method does.) The proof is in the **Bfs-oct12 / Fundamental Theorem** blocks.

---

### 4. Worked example

$x_1+x_2+x_3=4$, $x_1-x_2+x_4=2$, $x\ge0$ ($m=2$, $n=4$). $A=\begin{pmatrix}1&1&1&0\\1&-1&0&1\end{pmatrix}$.
- Basis $\{x_3,x_4\}$: $B=I$, $x_B=(4,2)$ ⇒ $x=(0,0,4,2)$, a b.f.s.
- Basis $\{x_1,x_2\}$: $B=\begin{pmatrix}1&1\\1&-1\end{pmatrix}$; solve $x_1+x_2=4$, $x_1-x_2=2$ ⇒ $x_1=3,x_2=1$ ⇒ $x=(3,1,0,0)$, a b.f.s.
- Basis $\{x_2,x_4\}$: $x_2=4$, $-x_2+x_4=2\Rightarrow x_4=6$ ⇒ $(0,4,0,6)$, a b.f.s.
- Basis $\{x_1,x_3\}$: $x_1=2$ (from row 2), $x_3=4-2=2$ ⇒ $(2,0,2,0)$, a b.f.s.
- Basis $\{x_2,x_3\}$: row 2 gives $-x_2=2\Rightarrow x_2=-2<0$ ⇒ basic but **not feasible**.

---

### 5. Summary

✅ Outcomes: optimal / unbounded / infeasible.
✅ Standard form via slack ($+z$), surplus ($-z$) and free ($x^+-x^-$) variables.
✅ Basic solution: choose $m$ independent columns $B$, set $x_N=0$, $x_B=B^{-1}b$. B.f.s.: also $x_B\ge0$.
✅ Fundamental theorem: feasible ⇒ a b.f.s. exists; optimal ⇒ an optimal b.f.s. exists.

<!-- LECTURE END: lp-bfs -->

---

<!-- LECTURE START: lp-geometric -->
<a id="lp-geometric"></a>

## Lecture 25: Linear Programs: the geometrical (graphical) method

*Source file:* `Lp-geometric.pdf`

### 1. What this lecture is about

For LPs with 2 variables we can **draw** the feasible region and slide the objective line to find the optimum. This builds intuition for the key fact: **an LP optimum (when it exists) occurs at a corner (vertex) of the feasible region.**

---

### 2. Background you need first

- Plotting a line $a_1x_1+a_2x_2=b$: find its intercepts with the axes. The inequality $\le b$ is one side of the line (test the origin to see which side).
- Level sets (contours) of a linear function are parallel lines.

---

### 3. Main content with full derivation

#### 3.1 The example

$$\min\ -x_1-x_2\quad\text{s.t.}\quad x_1+2x_2\le3,\quad 2x_1+x_2\le3,\quad x_1\ge0,\ x_2\ge0 .$$

**Step 1: draw each constraint.**
- $x_1+2x_2=3$: intercepts $(3,0)$ and $(0,1.5)$. Origin: $0\le3$ ✓, so the feasible side contains the origin.
- $2x_1+x_2=3$: intercepts $(1.5,0)$ and $(0,3)$. Origin side again.
- $x_1,x_2\ge0$: first quadrant.

**Step 2: feasible region** = intersection (the hatched quadrilateral). 🔍 Corners:
- $(0,0)$;
- $(1.5,0)$ (on $2x_1+x_2=3$ and $x_2=0$);
- $(0,1.5)$ (on $x_1+2x_2=3$ and $x_1=0$);
- intersection of the two lines: subtract $2\times$(first) from (second)... simpler: from $x_1+2x_2=3$ and $2x_1+x_2=3$, adding gives $3x_1+3x_2=6\Rightarrow x_1+x_2=2$; subtracting gives $x_1-x_2=0$; so $(1,1)$.

**Step 3: objective contours.** $-x_1-x_2=k$ ⟺ $x_2=-x_1-k$: parallel lines of slope $-1$.
- $k=0$: $x_2=-x_1$ (through the origin);
- $k=1$: $x_2=-x_1-1$; $k=-1$: $x_2=-x_1+1$; and so on.
Decreasing $k$ moves the line **up and to the right**, in the direction $(1,1)$ (the arrow drawn), which is $-\nabla(-x_1-x_2)=-(-1,-1)$.

**Step 4: slide** the line in the direction $(1,1)$ as far as possible while still touching the feasible region. The last point touched is the corner $(1,1)$.

🔍 *Check by evaluating at every corner:*

| corner | $-x_1-x_2$ |
|---|---|
| $(0,0)$ | $0$ |
| $(1.5,0)$ | $-1.5$ |
| $(0,1.5)$ | $-1.5$ |
| $(1,1)$ | $\mathbf{-2}$ |

**Optimal solution $x^*=(1,1)$, optimal value $-2$.** This LP is **feasible and bounded**.

#### 3.2 What can happen in general

- **Unique optimum at a vertex** (as here).
- **Infinitely many optima along an edge**, if the objective line is parallel to an edge of the region.
- **Unbounded**, if the region extends forever in the improving direction.
- **Infeasible**, if the constraints have no common point.

---

### 4. Summary

✅ Draw the constraints, shade the feasible polygon, slide the objective line in the direction $-c$ (for min).
✅ The optimum is at a corner; here $x^*=(1,1)$, value $-2$.
✅ This geometric "corner" idea becomes "basic feasible solution" in algebra.

<!-- LECTURE END: lp-geometric -->

---

<!-- LECTURE START: bfs-oct12 -->
<a id="bfs-oct12"></a>

## Lecture 26: Basic Feasible Solutions = Extreme Points (with examples and proof)

*Source file:* `Bfs-oct12.pdf`

### 1. What this lecture is about

1. Recap of standard form, slack/surplus/free variables and basic solutions (see the LP Solutions block)
2. A full worked example computing basic solutions of a small LP (feasible and infeasible ones)
3. **Extreme points** of a set
4. **Theorem: b.f.s. ⟺ extreme point** of the feasible set $\{x:Ax=b,x\ge0\}$, proved in both directions. This links the algebra (b.f.s.) to the geometry (corners).

---

### 2. Background you need first

- Basic solution: choose $m$ independent columns $B$, set $x_N=0$, solve $Bx_B=b$. B.f.s. if also $x_B\ge0$ (LP Solutions block).
- If columns $a_1,\dots,a_p$ are linearly **dependent**, there are $\beta_i$ not all zero with $\sum\beta_ia_i=0$; if independent, $\sum\beta_ia_i=b$ has at most one solution.

---

### 3. Main content with full derivations

#### 3.1 Worked example: computing basic solutions

$$\min3x_1+5x_2\quad\text{s.t.}\quad x_1+2x_2\le4,\quad2x_1+x_2\le6,\quad x_1,x_2\ge0 .$$
(Feasible region: the hatched quadrilateral with corners $(0,0)$, $(3,0)$, $(0,2)$ and the intersection of the two lines.)

**Standard form** (slack $x_3,x_4$):
$$\min3x_1+5x_2+0x_3+0x_4\ \text{ s.t. }\ x_1+2x_2+x_3=4,\ \ 2x_1+x_2+x_4=6,\ \ x_i\ge0 .$$
$$A=\begin{pmatrix}1&2&1&0\\2&1&0&1\end{pmatrix},\quad b=\begin{pmatrix}4\\6\end{pmatrix},\quad m=2,\ n=4 .$$

**(1) Basis = columns 1,2** ($x_B=(x_1,x_2)$, $x_N=(x_3,x_4)=0$):
$x_1+2x_2=4$, $2x_1+x_2=6$. 🔍 From the first, $x_1=4-2x_2$; substitute: $8-4x_2+x_2=6\Rightarrow x_2=\frac23$, $x_1=4-\frac43=\frac83$. So $x=(\frac83,\frac23,0,0)\ge0$: **a b.f.s.** (the intersection corner).

**(2) Basis = columns 3,4** ($B=I$): $x_1=x_2=0$, $x_3=4$, $x_4=6$: **a b.f.s.** (the origin).

**(3) Basis = columns 2,3:** $x_1=x_4=0$: $2x_2+x_3=4$, $x_2=6$ ⇒ $x_3=4-12=-8<0$. **A basic solution, but not a b.f.s.** (Geometrically $(0,6)$ is where the line $2x_1+x_2=6$ meets the $x_2$-axis, outside the feasible region.)

🔍 The remaining bases give $(0,2,0,4)$ (columns 2,4: b.f.s.), $(3,0,1,0)$ (columns 1,4? check: $x_1=4$ from row 1 with $x_2=x_3=0$ gives $x_4=6-8=-2$, not feasible; columns 1,3: $2x_1=6\Rightarrow x_1=3$, $x_3=4-3=1$, b.f.s. $(3,0,1,0)$). So the 4 b.f.s. are exactly the 4 corners. The optimum of this min problem with positive costs is $(0,0)$, value 0.

#### 3.2 Extreme points

$x\in\Omega$ is an **extreme point** of $\Omega$ if it **cannot** be written as
$$x=\alpha y+(1-\alpha)z,\qquad y,z\in\Omega,\ y\neq z,\ \alpha\in(0,1).$$
I.e. $x$ is not strictly inside any segment of $\Omega$. Corners of a polygon are extreme points; edge midpoints are not.

#### 3.3 Theorem: b.f.s. ⟺ extreme point of $\Omega=\{x:Ax=b,\ x\ge0\}$

**(⇒) Every b.f.s. is an extreme point.**
1. Let $x=(x_1,\dots,x_m,0,\dots,0)^T$ be a b.f.s. with $A=[B\ N]$, columns $a_1,\dots,a_m$ of $B$ linearly independent, $x_i\ge0$.
2. Suppose $x$ is **not** extreme: $x=\alpha y+(1-\alpha)z$ with $0<\alpha<1$, $y,z\in\Omega$ ($Ay=b$, $Az=b$, $y,z\ge0$).
3. For each non-basic index $i>m$: $0=x_i=\alpha y_i+(1-\alpha)z_i$, a sum of two non-negative terms with positive weights, so $y_i=z_i=0$. Hence $y=(y_1,\dots,y_m,0,\dots,0)^T$ and $z=(z_1,\dots,z_m,0,\dots,0)^T$.
4. Then $Ay=b$ reads $\sum_{i=1}^ma_iy_i=b$ and $Az=b$ reads $\sum_{i=1}^ma_iz_i=b$, i.e. $By_B=b$ and $Bz_B=b$.
5. But $\operatorname{rank}(B)=m$ (invertible), so $Bq=b$ has a **unique** solution: $y_B=z_B=x_B$. So $y=z=x$, contradicting $y\ne z$. ✗ Hence $x$ is extreme. ∎

**(⇐) Every extreme point is a b.f.s.**
1. Let $x$ be an extreme point of the constraint set; reorder so $x=(x_1,\dots,x_p,0,\dots,0)^T$ with $x_1,\dots,x_p>0$. Then $\sum_{i=1}^px_ia_i=b$.
2. If $a_1,\dots,a_p$ are linearly independent, then $p\le m$ and (adding more independent columns to make $m$ if needed) $x$ **is already a b.f.s.** ✓
3. If not independent: $\exists\beta_1,\dots,\beta_p$, not all 0, with $\sum_{i=1}^p\beta_ia_i=0$. Let $\beta=(\beta_1,\dots,\beta_p,0,\dots,0)^T$, so $A\beta=0$.
4. Consider $x+\varepsilon\beta$ and $x-\varepsilon\beta$. Both satisfy $A(x\pm\varepsilon\beta)=b\pm0=b$ ✓. Are they $\ge0$? Need $x_i\pm\varepsilon\beta_i\ge0$ for $i\le p$. Since every $x_i>0$, this holds for small enough $\varepsilon$. (HW: a valid choice is $\varepsilon<\min_{i:\beta_i\ne0}\frac{x_i}{|\beta_i|}$.)
5. Then $x=\frac12(x+\varepsilon\beta)+\frac12(x-\varepsilon\beta)$ with both points feasible and different ($\beta\neq0$). So $x$ is **not** extreme, a contradiction. ✗ So the columns must be independent, and $x$ is a b.f.s. ∎

---

### 4. Summary

✅ Compute basic solutions by choosing $m$ columns, zeroing the rest, solving $Bx_B=b$; feasible iff $x_B\ge0$.
✅ Extreme point: not a strict convex combination of two other points of the set.
✅ **b.f.s. ⟺ extreme point** (uniqueness of $Bx_B=b$ for ⇒; the perturbation $x\pm\varepsilon\beta$ for ⇐).
✅ So LP corners = b.f.s., and the Fundamental Theorem says: look only at corners.

<!-- LECTURE END: bfs-oct12 -->

---

<!-- LECTURE START: oct14-fundatheorem -->
<a id="oct14-fundatheorem"></a>

## Lecture 27: Proof of the Fundamental Theorem of Linear Programming

*Source file:* `Oct14-fundaTheorem.pdf`

### 1. What this lecture is about

We prove, for the standard-form LP
$$\min c^Tx\quad\text{s.t.}\quad Ax=b,\ x\ge0,\qquad A\in\mathbb{R}^{m\times n},\ \operatorname{rank}(A)=m,\ m<n,$$
1. **Part 1:** if there is a feasible solution, there is a basic feasible solution (b.f.s.);
2. **Part 2:** if there is an optimal solution, there is an optimal b.f.s.

Both proofs use the same trick: **if the columns used are dependent, slide along a null-space direction until one variable hits zero**, and repeat.

---

### 2. Background you need first

- Basic feasible solution (LP Solutions and BFS blocks).
- Columns $a_i$ of $A$: $Ax=\sum_ix_ia_i$.
- If more than $m$ vectors live in $\mathbb{R}^m$, they are linearly dependent.

---

### 3. Main content with full derivations

#### 3.1 Part 1: feasible ⇒ b.f.s. exists

1. Let $x$ be feasible. Reorder variables so that its non-zero entries come first:
$$x=(x_1,x_2,\dots,x_p,0,\dots,0)^T,\qquad x_i>0\ (i\le p),\quad1\le p\le n .$$
Feasibility: $\sum_{i=1}^px_ia_i=b$.
2. **If $a_1,\dots,a_p$ are linearly independent:** then $p\le m$ and $x$ is a b.f.s. (extend $\{a_1..a_p\}$ to $m$ independent columns using $\operatorname{rank}A=m$; the extra basic variables are 0). **Done.**
3. **If not:** $\exists\alpha_1,\dots,\alpha_p$, not all 0, with
$$\sum_{i=1}^p\alpha_ia_i=0 .$$
Let $\alpha=(\alpha_1,\dots,\alpha_p,0,\dots,0)^T$, so $A\alpha=0$. WLOG at least one $\alpha_i>0$ (otherwise replace $\alpha$ by $-\alpha$).
4. Define $\tilde x=x-\varepsilon\alpha$. Then $A\tilde x=b-\varepsilon\cdot0=b$ ✓. We need $\tilde x\ge0$: $x_i-\varepsilon\alpha_i\ge0$ for all $i$.
   - If $\alpha_i\le0$: automatic for $\varepsilon\ge0$.
   - If $\alpha_i>0$: need $\varepsilon\le\frac{x_i}{\alpha_i}$.
5. **Take** $\displaystyle\varepsilon=\min_{i:\alpha_i>0}\Big\{\frac{x_i}{\alpha_i}\Big\}$ (positive).
6. Then $\tilde x\ge0$, $\sum_{i=1}^p(x_i-\varepsilon\alpha_i)a_i=b$, and **at least one** $x_i-\varepsilon\alpha_i=0$ (the index achieving the min). So $\tilde x=(\tilde x_1,\dots,\tilde x_{p-1},0,\dots,0)$ after reordering: feasible with **at most $p-1$** non-zeros.
7. **Repeat inductively.** Each round removes at least one non-zero. The process stops when the remaining columns are independent (at the latest when $p\le m$... more precisely when independence holds). The result is a b.f.s. ∎

#### 3.2 Part 2: optimal ⇒ optimal b.f.s. exists

1. Let $x^*$ be optimal: $Ax^*=b$, $x^*\ge0$, and $c^Tx^*\le c^Tx$ for all feasible $x$. Write $x^*=(x_1^*,\dots,x_p^*,0,\dots,0)^T$ with positive entries first.
2. If $a_1,\dots,a_p$ are independent, $x^*$ is already an optimal b.f.s.
3. Otherwise $\sum_{i=1}^p\alpha_ia_i=0$ for some non-zero $\alpha$. **Same steps as before**: $x^*-\varepsilon\alpha$ is feasible for small $|\varepsilon|$ (positive *or* negative, since all $x_i^*>0$).
4. **Claim: $c^T\alpha=0$.**
   - By optimality, $c^T(x^*-\varepsilon\alpha)\ge c^Tx^*$, i.e. $c^Tx^*-\varepsilon c^T\alpha\ge c^Tx^*$, i.e. $\varepsilon\,c^T\alpha\le0$.
   - This holds for small $\varepsilon>0$ **and** small $\varepsilon<0$. 🔍 With $\varepsilon>0$: $c^T\alpha\le0$. With $\varepsilon<0$: $c^T\alpha\ge0$. So $c^T\alpha=0$.
   - (In the notes' words: otherwise you could construct $x^*+\tilde\varepsilon\alpha$ with a better objective.)
5. Hence $x^*-\varepsilon\alpha$ has the **same** objective value $c^Tx^*-\varepsilon\cdot0=c^Tx^*$: it is also optimal.
6. Choose $\varepsilon=\min_{i:\alpha_i>0}\frac{x_i^*}{\alpha_i}$ as in Part 1. Then $x^*-\varepsilon\alpha$ is optimal with **$p-1$ non-zeros**. Repeat until the columns are independent: an **optimal b.f.s.** ∎

---

### 4. Why this matters

There are at most $\binom nm$ b.f.s., so an LP can be solved (in principle) by checking finitely many points. The **simplex method** moves cleverly from one b.f.s. to a neighbouring better one.

---

### 5. Summary

✅ Non-zeros of $x$ use columns $a_1..a_p$. If dependent, pick $\alpha$ with $A\alpha=0$ and move to $x-\varepsilon\alpha$ with $\varepsilon=\min_{\alpha_i>0}x_i/\alpha_i$, killing one non-zero.
✅ Feasibility is preserved ⇒ Part 1.
✅ For an optimal $x^*$, $c^T\alpha=0$ (else we could improve), so optimality is preserved ⇒ Part 2.

<!-- LECTURE END: oct14-fundatheorem -->

---

<!-- LECTURE START: oct14 -->
<a id="oct14"></a>

## Lecture 28: LP Duality: introduction (Oct 14 class, partial)

*Source file:* `Oct14.pdf`

> 📌 **This file is the first half of `Oct16.pdf`.** The complete explanation is in the block **"LP Duality: deriving the dual, primal–dual pairs, weak duality, transportation & diet duals"**. That block includes the numerical example solved fully (both optimal values equal 1000), all primal–dual rules, and the transportation and diet duals. This block is a checklist.

### Checklist

1. **Example primal:** $\min2x_1+5x_2$ s.t. $4x_1+10x_2\ge100$, $x_1-7x_2\ge500$, $x\ge0$.
2. Multiplying constraint 1 by $0.5$ gives $2x_1+5x_2\ge50$. Multiplying both by $0.1$ and adding gives $0.5x_1+0.3x_2\ge60$, and since $x\ge0$, $2x_1+5x_2\ge60$.
3. **General multipliers** $y\ge0$ with $4y_1+y_2\le2$, $10y_1-7y_2\le5$ give OPT $\ge100y_1+500y_2$. Maximising this bound is the **dual**.
4. $\min c^Tx$, $Ax\ge b$, $x\ge0$ ⟷ $\max b^Ty$, $A^Ty\le c$, $y\ge0$, because $c^Tx\ge y^TAx\ge y^Tb$.
5. $\max c^Tx$, $Ax\le b$, $x\ge0$ ⟷ $\min b^Ty$, $A^Ty\ge c$, $y\ge0$, because $c^Tx\le y^TAx\le b^Ty$.
6. The class ended with the standard form $\min c^Tx$, $Ax=b$, $x\ge0$, whose dual has a free $y$.

<!-- LECTURE END: oct14 -->

---

<!-- LECTURE START: oct16 -->
<a id="oct16"></a>

## Lecture 29: LP Duality: deriving the dual, primal–dual pairs, weak duality, transportation & diet duals

*Source file:* `Oct16.pdf`

### 1. What this lecture is about

Every LP (the **primal**) has a partner LP (the **dual**) whose optimal value bounds it. This lecture:
1. Discovers the dual from scratch on a numerical example: combine the constraints with multipliers $y\ge0$ to get a lower bound on the objective, then find the *best* such bound
2. General primal–dual pairs for different constraint/sign types
3. **Weak duality:** $c^Tx\ge b^Ty$ (with a one-line proof)
4. Economic interpretation: the duals of the transportation and diet problems

---

### 2. Background you need first

- LP forms (LP blocks). For vectors, $u\le v$ means componentwise.
- If $x\ge0$ and $u\le v$ componentwise, then $u^Tx\le v^Tx$ (each term $u_ix_i\le v_ix_i$).

---

### 3. Main content with full derivations

#### 3.1 Discovering the dual: an example

**Primal:**
$$\min\ 2x_1+5x_2\quad\text{s.t.}\quad4x_1+10x_2\ge100,\quad x_1-7x_2\ge500,\quad x_1,x_2\ge0 .$$
Goal: find a **lower bound** on the optimal value without solving.

- **Try multiplier $y_1=0.5$ on constraint 1, $y_2=0$:** $0.5(4x_1+10x_2)\ge0.5\cdot100$, i.e. $2x_1+5x_2\ge50$. That is *exactly* the objective, so OPT $\ge50$.
- **Try $0.1$ on each:** $0.4x_1+x_2\ge10$ and $0.1x_1-0.7x_2\ge50$. Adding: $0.5x_1+0.3x_2\ge60$. Since $x\ge0$ and $2\ge0.5$, $5\ge0.3$:
$$2x_1+5x_2\ \ge\ 0.5x_1+0.3x_2\ \ge\ 60 .$$
A better bound: OPT $\ge60$.

**General multipliers $y_1,y_2$.** Multiply constraint $i$ by $y_i$ and add:
$$(4y_1+y_2)x_1+(10y_1-7y_2)x_2\ \ge\ 100y_1+500y_2 .$$
Conditions for this to bound the objective:
- $y_1,y_2\ge0$ (multiplying a $\ge$ inequality by a negative number flips it);
- the coefficients must not exceed the objective's: $4y_1+y_2\le2$ and $10y_1-7y_2\le5$ (then, since $x\ge0$, $2x_1+5x_2\ge(4y_1+y_2)x_1+(10y_1-7y_2)x_2$).

Then OPT $\ge100y_1+500y_2$. **The best such bound** is the **dual LP**:
$$\max\ 100y_1+500y_2\quad\text{s.t.}\quad4y_1+y_2\le2,\quad10y_1-7y_2\le5,\quad y_1,y_2\ge0 .$$

🔍 *Solving both (not in notes).* Primal: constraint 2 gives $x_1\ge500+7x_2$; cost $2x_1+5x_2\ge1000+19x_2$, minimised at $x_2=0,x_1=500$: **OPT = 1000** (constraint 1: $2000\ge100$ ✓). Dual: $y=(0,2)$: $4\cdot0+2=2\le2$ ✓, $0-14\le5$ ✓, value $500\cdot2=$ **1000**. The best lower bound equals the optimum! (Strong duality, next lecture.)

#### 3.2 General form

$$\textbf{(P)}\ \min c^Tx\ \text{ s.t. }Ax\ge b,\ x\ge0\qquad\Longleftrightarrow\qquad\textbf{(D)}\ \max b^Ty\ \text{ s.t. }A^Ty\le c,\ y\ge0\quad(y\in\mathbb{R}^m).$$
One dual variable $y_i$ per primal constraint; one dual constraint per primal variable.

**Weak duality proof.** For any primal-feasible $x$ and dual-feasible $y$:
$$c^Tx\ \ge\ (A^Ty)^Tx=y^TAx\ \ge\ y^Tb .$$
- First $\ge$: $c\ge A^Ty$ and $x\ge0$.
- Second $\ge$: $Ax\ge b$ and $y\ge0$.

So **every dual-feasible value ≤ every primal-feasible value**.

**Max form:**
$$\max c^Tx\ \text{ s.t. }Ax\le b,\ x\ge0\qquad\Longleftrightarrow\qquad\min b^Ty\ \text{ s.t. }A^Ty\ge c,\ y\ge0,$$
with $c^Tx\le y^T(Ax)\le b^Ty$.

#### 3.3 Other primal forms (table of rules)

| | Primal (min) | Dual (max) |
|---|---|---|
| 1 | $Ax\ge b,\ x\ge0$ | $A^Ty\le c,\ y\ge0$ |
| 2 | $Ax=b,\ x\ge0$ | $A^Ty\le c,\ y$ **unrestricted** |
| 3 | $Ax\ge b,\ x$ unrestricted | $A^Ty=c,\ y\ge0$ |
| 4 | $Ax=b,\ x$ unrestricted | $A^Ty=c,\ y$ unrestricted |

**Pattern:** equality constraint ⟷ free dual variable; free primal variable ⟷ equality dual constraint.

🔍 *Why (case 2):* we want $b^Ty=y^TAx\le c^Tx$ for every feasible $x$. Since $Ax=b$ is an equality, $y$ can have any sign (no inequality is flipped). We only need $y^TAx\le c^Tx$ for $x\ge0$, i.e. $A^Ty\le c$.
🔍 *Why (case 3):* $x$ can be negative, so $(c-A^Ty)^Tx\ge0$ for *all* $x$ forces $c-A^Ty=0$.

#### 3.4 Dual of the transportation problem

$m$ factories (supply $a_i$), $n$ warehouses (demand $b_j$), unit cost $c_{ij}$ (Maruti Suzuki's problem).
**Primal** (write $\le a_i$ as $-\sum_jx_{ij}\ge-a_i$ to get all "$\ge$"):
$$\min\sum_{i,j}c_{ij}x_{ij}\ \text{ s.t. }\ \sum_ix_{ij}\ge b_j\ [v_j],\quad-\sum_jx_{ij}\ge-a_i\ [u_i],\quad x_{ij}\ge0 .$$
**Dual** (multipliers $v_j\ge0$, $u_i\ge0$). Variable $x_{ij}$ appears in demand row $j$ (coefficient $+1$) and supply row $i$ (coefficient $-1$), so its dual constraint is $v_j-u_i\le c_{ij}$:
$$\max\ \sum_{j=1}^nv_jb_j-\sum_{i=1}^ma_iu_i\quad\text{s.t.}\quad v_j-u_i\le c_{ij},\quad u_i\ge0,\ v_j\ge0 .$$
**Interpretation** (an external shipping agency's problem): $v_j$ = what the agency **charges** to deliver one unit to warehouse $j$; $u_i$ = what it **pays** factory $i$ per unit picked up. Its margin on route $(i,j)$ is $v_j-u_i$, which must not exceed Maruti's own cost $c_{ij}$ (else Maruti ships itself). The agency maximises its profit.

#### 3.5 Dual of the diet problem

$m$ foods with cost $c_i$; $n$ nutrients with requirement $b_j$; one unit of food $i$ gives $a_{ij}$ of nutrient $j$.
**Primal:** $\min\sum_ic_ix_i$ s.t. $\sum_ia_{ij}x_i\ge b_j\ \forall j$, $x_i\ge0$ ($x_i$ = quantity of food).
**Dual:**
$$\max\sum_{j=1}^nb_jy_j\quad\text{s.t.}\quad\sum_{j=1}^na_{ij}y_j\le c_i\ \ \forall i,\quad y_j\ge0 .$$
**Interpretation:** a pill-maker sells pure nutrients at price $y_j$ per unit. To compete, the "pill price" of each food's nutrient content, $\sum_ja_{ij}y_j$, must not exceed that food's price $c_i$. The maker maximises revenue $\sum_jb_jy_j$.

---

### 4. Summary

✅ Dual = best lower bound obtained by non-negative combinations of the constraints.
✅ $\min c^Tx$, $Ax\ge b$, $x\ge0$ ⟷ $\max b^Ty$, $A^Ty\le c$, $y\ge0$.
✅ Weak duality: $c^Tx\ge y^TAx\ge b^Ty$.
✅ Equality ⟷ free variable. Dual variables are prices (transportation: charges/payments; diet: nutrient prices).

<!-- LECTURE END: oct16 -->

---

<!-- LECTURE START: dualitytheoremlp -->
<a id="dualitytheoremlp"></a>

## Lecture 30: Strong Duality Theorem for LPs (with proof via Farkas)

*Source file:* `DualityTheoremLP.pdf`

### 1. What this lecture is about

Weak duality says dual value ≤ primal value. **Strong duality** says that when an optimum exists, the two values are **equal** (zero duality gap). We:
1. Define the duality gap
2. Tabulate which primal/dual outcome combinations are possible
3. State a Farkas-type alternative system
4. Prove the **Duality Theorem for LPs** step by step

---

### 2. Background you need first

- Dual pair: $p^*=\min\{c^Tx:Ax\ge b,\ x\ge0\}$ and $d^*=\max\{b^Ty:A^Ty\le c,\ y\ge0\}$ (LP Duality block).
- Farkas' lemma: exactly one of $\{Ax=b,x\ge0\}$ and $\{A^Ty\le0,b^Ty>0\}$ is feasible (FL block).
- Weak duality: $c^Tx\ge y^TAx\ge y^Tb$.

---

### 3. Main content with full derivations

#### 3.1 Duality gap

$$p^*=\min c^Tx\ \text{ s.t. }Ax\ge b,\ x\ge0\ \ [\text{multipliers }y],\qquad d^*=\max b^Ty\ \text{ s.t. }A^Ty\le c,\ y\ge0 .$$
**Duality gap** $=p^*-d^*\ge0$ **always** (weak duality: $c^Tx\ge y^TAx\ge y^Tb$ for all feasible pairs, so $p^*\ge d^*$).

#### 3.2 Possible outcome combinations

| Primal (min) \ Dual (max) | Finite opt | Unbounded | Infeasible |
|---|---|---|---|
| **Finite optimum** | ✓ (equal values) | ✗ | ✗ |
| **Unbounded** | ✗ | ✗ | ✓ |
| **Infeasible** | ✗ | ✓ | ✓ |

🔍 *Why some are impossible:* if the primal is unbounded ($p^*=-\infty$), any dual-feasible $y$ would give $b^Ty\le-\infty$, impossible, so the dual is **infeasible**. Similarly a dual unbounded ($+\infty$) forces the primal infeasible. Both infeasible *can* happen (e.g. contradictory constraints on both sides).

#### 3.3 An alternative system (Farkas variant)

Exactly one of these is feasible:
1. $\{y\in\mathbb{R}^m:\ A^Ty\le c,\ y\ge0\}$
2. $\{q\in\mathbb{R}^n:\ Aq\ge0,\ q\ge0,\ c^Tq<0\}$

🔍 *Derivation from Farkas.* Add slacks: (1) ⟺ $A^Ty+s=c$, $y,s\ge0$, i.e. $[A^T\ \ I]\binom ys=c$, $\binom ys\ge0$. Farkas with matrix $M=[A^T\ I]$ and right side $c$: the alternative is $\exists q$: $M^Tq\le0$, $c^Tq>0$, i.e. $Aq\le0$, $q\le0$, $c^Tq>0$. Replace $q$ by $-q$: $Aq\ge0$, $q\ge0$, $c^Tq<0$. ✓

#### 3.4 Duality Theorem for LPs

> 1. If either the primal or the dual has a finite optimum, then so does the other, and the **duality gap is 0** ($p^*=d^*$).
> 2. If the primal is unbounded, the dual is infeasible (and vice versa).

**Proof of 1** (primal has finite optimum $p^*$, attained at $x^*$).

*Goal:* show that
$$\{y\in\mathbb{R}^m:\ A^Ty\le c,\ y\ge0,\ b^Ty\ge p^*\}$$
is **feasible**. Then such a $y$ has $b^Ty\ge p^*$, while weak duality gives $b^Ty\le p^*$. So $d^*=p^*$ and the gap is 0.

1. **Assume it is infeasible.** Rewrite it as one system "$\tilde A^Ty\le\tilde c$, $y\ge0$" by turning $b^Ty\ge p^*$ into $-b^Ty\le-p^*$:
$$\tilde A^T=\begin{pmatrix}A^T\\-b^T\end{pmatrix}\in\mathbb{R}^{(n+1)\times m},\qquad\tilde c=\begin{pmatrix}c\\-p^*\end{pmatrix}\in\mathbb{R}^{n+1}.$$
2. By the alternative (§3.3) with $\tilde A,\tilde c$: $\exists q\in\mathbb{R}^{n+1}$ with $\tilde Aq\ge0$, $q\ge0$, $\tilde c^Tq<0$.
3. Write $q=\binom{w}{q_0}$, $w\in\mathbb{R}^n$, $q_0\in\mathbb{R}$. Since $\tilde A=[A\ \ -b]$ ($m\times(n+1)$):
$$Aw-bq_0\ge0\ (\text{i.e. }Aw\ge bq_0),\qquad w\ge0,\ q_0\ge0,\qquad c^Tw-q_0p^*<0 .$$
4. **Case $q_0=0$:** then $Aw\ge0$, $w\ge0$, $c^Tw<0$. Take the optimal $x^*$ ($Ax^*\ge b$, $x^*\ge0$, $c^Tx^*=p^*$). Then $x^*+w$ satisfies $A(x^*+w)\ge b+0$ and $x^*+w\ge0$: **feasible**. But $c^T(x^*+w)=p^*+c^Tw<p^*$: **better than optimal**. ✗
5. **Case $q_0>0$:** divide by $q_0$. Let $x=w/q_0$: $Ax\ge b$, $x\ge0$ (feasible), and $c^Tx=\frac{c^Tw}{q_0}<p^*$. Again better than optimal. ✗
6. Both cases are impossible, so the system **is feasible**, and duality gap $=0$. ∎

(The "dual has finite optimum" case is symmetric: the dual of the dual is the primal.)

---

### 4. Worked check

From the LP Duality block: primal $\min2x_1+5x_2$ ($4x_1+10x_2\ge100$, $x_1-7x_2\ge500$, $x\ge0$) has $p^*=1000$ at $(500,0)$; the dual $y=(0,2)$ gives $1000$. Gap $0$ ✓.

---

### 5. Summary

✅ Duality gap $p^*-d^*\ge0$ always.
✅ Strong duality: finite optimum on one side ⇒ finite optimum on the other, with equal values.
✅ Proof: if no dual $y$ reached $p^*$, Farkas would give a direction $w$ or point $w/q_0$ beating the primal optimum.
✅ Unbounded primal ⇒ infeasible dual (and vice versa).

<!-- LECTURE END: dualitytheoremlp -->

---

<!-- LECTURE START: compslack-lp -->
<a id="compslack-lp"></a>

## Lecture 31: Complementary Slackness, certificates of optimality, and fractional knapsack

*Source file:* `CompSlack-LP.pdf`

### 1. What this lecture is about

Strong duality ($c^Tx^*=b^Ty^*$) has a powerful consequence, **complementary slackness (CS)**: at optimality, a positive primal variable forces its dual constraint to be tight, and vice versa. This lecture:
1. Derives CS from strong duality
2. Uses CS as a **certificate of optimality**
3. Solves a numerical LP's dual using only CS
4. Solves the **fractional knapsack** problem and proves the greedy solution optimal via its dual

---

### 2. Background you need first

- Primal $p^*=\min\{c^Tx:Ax\ge b,\ x\ge0\}$, dual $d^*=\max\{b^Ty:A^Ty\le c,\ y\ge0\}$; weak & strong duality (previous blocks).
- A sum of non-negative terms is 0 only if every term is 0.
- $(y^TA)_i$ = $i$-th entry of the row vector $y^TA$ = $a_i^Ty$, where $a_i$ is column $i$ of $A$.

---

### 3. Main content with full derivations

#### 3.1 Deriving complementary slackness

For any feasible $x,y$: $c^Tx\ge y^TAx\ge y^Tb$. At optimal $x^*,y^*$, strong duality ($p^*=d^*$) forces **both** inequalities to be equalities:
$$c^Tx^*=y^{*T}Ax^*=y^{*T}b .$$

**From the first equality:**
$$\sum_{i=1}^nc_ix_i^*=\sum_{i=1}^n(y^{*T}A)_ix_i^*\ \Rightarrow\ \sum_{i=1}^n\underbrace{\big(c_i-(y^{*T}A)_i\big)}_{\ge0\ (\text{dual feas.})}\underbrace{x_i^*}_{\ge0}=0 .$$
A sum of non-negative terms equals 0 ⇒ each term is 0:
$$\boxed{\big(c_i-(y^{*T}A)_i\big)\,x_i^*=0\quad\forall i=1..n}$$
- If $x_i^*>0$ ⇒ $c_i-(y^{*T}A)_i=0$ (dual constraint $i$ tight).
- If $c_i-(y^{*T}A)_i>0$ (slack) ⇒ $x_i^*=0$.

**From the second equality** $y^{*T}Ax^*=y^{*T}b$:
$$\sum_{i=1}^m\underbrace{y_i^*}_{\ge0}\underbrace{(Ax^*-b)_i}_{\ge0}=0\ \Rightarrow\ \boxed{y_i^*\,(Ax^*-b)_i=0\quad\forall i=1..m}$$
- $y_i^*>0$ ⇒ $(Ax^*-b)_i=0$ (primal constraint $i$ tight).
- $(Ax^*-b)_i>0$ ⇒ $y_i^*=0$.

"Complementary": in each pair, at least one of the two is zero (a variable and the slack of its partner constraint).

#### 3.2 Certificate of optimality

> If $x$ is primal feasible ($x\ge0$, $Ax\ge b$), $y$ is dual feasible ($y\ge0$, $A^Ty\le c$), and they satisfy CS:
> $$\sum_{i=1}^n\big(c_i-(y^TA)_i\big)x_i=0\quad\text{and}\quad\sum_{i=1}^m(Ax-b)_iy_i=0,$$
> then $x$ is optimal for the primal **and** $y$ is optimal for the dual (and conversely).

🔍 *Why:* CS gives $c^Tx=y^TAx=y^Tb$. By weak duality every primal value is $\ge y^Tb=c^Tx$, so $x$ is optimal; similarly for $y$. So a pair $(x,y)$ is a short **proof** that $x$ is optimal: anyone can check feasibility and CS.

#### 3.3 Worked example

**Primal:**
$$\min3x_1+8x_2\ \text{ s.t. }\ x_1-2x_2\ge10\ [y_1],\quad5x_1+8x_2\ge2\ [y_2],\quad-x_1\ge-15\ [y_3],\quad x_1\ge0,\ x_2\in\mathbb{R}.$$
Given primal optimum $x_1^*=15$, $x_2^*=-\frac{73}{8}$.
🔍 *Check constraints:* $15+\frac{73}{4}=33.25>10$ (**slack**); $75-73=2$ (**tight**); $-15=-15$ (**tight**).

**Dual** ($x_1\ge0$ gives "$\le$", $x_2$ free gives "$=$"):
$$\max10y_1+2y_2-15y_3\ \text{ s.t. }\ y_1+5y_2-y_3\le3\ [x_1],\quad-2y_1+8y_2=8\ [x_2],\quad y_i\ge0 .$$
**Apply CS:**
- $x_1^*>0$ ⇒ dual constraint 1 tight: $y_1^*+5y_2^*-y_3^*=3$.
- Dual constraint 2 is an equality anyway: $-2y_1^*+8y_2^*=8$ (for any feasible $y$).
- Primal constraint 1 is slack ($33.25>10$) ⇒ $y_1^*=0$.

Solve: $8y_2^*=8\Rightarrow y_2^*=1$; $5-y_3^*=3\Rightarrow y_3^*=2$. So $y^*=(0,1,2)\ge0$ ✓.
**Values:** dual $10\cdot0+2\cdot1-15\cdot2=2-30=-28$; primal $3\cdot15+8\cdot(-\frac{73}{8})=45-73=-28$. Equal ✓, so both are optimal.

#### 3.4 Knapsack and fractional knapsack

**0/1 knapsack:** $n$ items; item $i$ has weight $w_i$, value $v_i$; bag capacity $B$. $x_i=1$ if item $i$ is included, else 0.
$$\max\sum_iv_ix_i\ \text{ s.t. }\sum_iw_ix_i\le B,\ x_i\in\{0,1\}\quad(\text{NP-hard}).$$

**Fractional knapsack** (LP relaxation: items can be cut):
$$\max\sum_{i=1}^nv_ix_i\ \text{ s.t. }\ \sum_iw_ix_i\le B\ [\lambda],\qquad0\le x_i\le1\ [\mu_i\text{ for }x_i\le1].$$
**Dual:**
$$\min\ \lambda B+\sum_{i=1}^n\mu_i\ \text{ s.t. }\ v_i\le\lambda w_i+\mu_i\ \ \forall i,\qquad\lambda\ge0,\ \mu_i\ge0 .$$

**Greedy solution.** Sort by value per weight: $\frac{v_1}{w_1}\ge\frac{v_2}{w_2}\ge\dots\ge\frac{v_n}{w_n}$. Fill greedily: $x_1=1$; if $w_1+w_2\le B$, $x_2=1$; … Let $k$ be the first item that does not fit: $w_1+\dots+w_{k-1}\le B<w_1+\dots+w_k$. Take a fraction of it:
$$\sum_{i=1}^{k-1}w_i+w_kx_k=B\ \Rightarrow\ x_k=\frac{B-\sum_{i=1}^{k-1}w_i}{w_k}<1 .$$
So $x=(1,\dots,1,x_k,0,\dots,0)$.

**Proving optimality via CS** (the notes start this; 🔍 completed here):
- $x_k<1$ ⇒ constraint $x_k\le1$ slack ⇒ $\mu_k=0$. Likewise $x_{k+1}=\dots=x_n=0<1$ ⇒ $\mu_{k+1}=\dots=\mu_n=0$.
- $x_k>0$ ⇒ dual constraint $k$ tight: $v_k=\lambda w_k+\mu_k=\lambda w_k$ ⇒ $\boxed{\lambda=\frac{v_k}{w_k}}$.
- For $i<k$, $x_i=1>0$ ⇒ tight: $\mu_i=v_i-\lambda w_i=w_i\big(\frac{v_i}{w_i}-\frac{v_k}{w_k}\big)\ge0$ ✓ (sorted order).
- Dual feasibility for $i>k$: $v_i\le\lambda w_i\iff\frac{v_i}{w_i}\le\frac{v_k}{w_k}$ ✓ (sorted).
- Values: dual $=\lambda B+\sum_{i<k}(v_i-\lambda w_i)=\sum_{i<k}v_i+\lambda\big(B-\sum_{i<k}w_i\big)=\sum_{i<k}v_i+\frac{v_k}{w_k}w_kx_k$ = primal value ✓.

Feasible primal and dual with equal values ⇒ the greedy solution is optimal.

---

### 4. Summary

✅ At optimality $c^Tx^*=y^{*T}Ax^*=y^{*T}b$ ⇒ $(c_i-(A^Ty^*)_i)x_i^*=0$ and $y_i^*(Ax^*-b)_i=0$.
✅ Feasible $x,y$ + CS ⟺ both optimal (certificate).
✅ Given $x^*$, CS often determines $y^*$ by linear equations.
✅ Fractional knapsack: sort by $v_i/w_i$, fill greedily; $\lambda=v_k/w_k$ certifies optimality.

<!-- LECTURE END: compslack-lp -->

---

<!-- LECTURE START: sensitivity-lp -->
<a id="sensitivity-lp"></a>

## Lecture 32: Sensitivity Analysis: dual variables as shadow prices

*Source file:* `Sensitivity-lp.pdf`

### 1. What this lecture is about

If we change a right-hand side $b_i$ of an LP a little, how much does the optimal value change? Answer: **by $\Delta\cdot y_i^*$**, where $y_i^*$ is the optimal dual variable. So dual variables are **shadow prices**: the marginal value of one more unit of resource $i$. We verify this on a 2-variable example.

---

### 2. Background you need first

- Graphical LP method (LP Geometric block); dual of $\max c^Tx$, $Ax\le b$, $x\ge0$ is $\min b^Ty$, $A^Ty\ge c$, $y\ge0$ (LP Duality block).
- Complementary slackness (previous block).

---

### 3. Main content with full derivations

#### 3.1 The example

$$\max3x_1+5x_2\quad\text{s.t.}\quad x_1\le4\ [y_1],\qquad2x_2\le12\ [y_2],\qquad3x_1+2x_2\le18\ [y_3],\qquad x_1,x_2\ge0 .$$
**Graph:** corners $(0,0)$, $(4,0)$, $(4,3)$, $(2,6)$, $(0,6)$.
🔍 Corner values: $0$; $12$; $12+15=27$; $6+30=\mathbf{36}$; $30$. **Optimum $x^*=(2,6)$, value 36.**
(At $(2,6)$: $x_1=2<4$, slack; $2x_2=12$, tight; $3\cdot2+2\cdot6=18$, tight.)

**Dual:**
$$\min4y_1+12y_2+18y_3\quad\text{s.t.}\quad y_1+3y_3\ge3\ [x_1],\qquad2y_2+2y_3\ge5\ [x_2],\qquad y_i\ge0 .$$
🔍 *Getting $y^*$ by CS:* constraint 1 is slack ⇒ $y_1^*=0$. $x_1^*,x_2^*>0$ ⇒ both dual constraints tight: $3y_3=3\Rightarrow y_3^*=1$; $2y_2+2=5\Rightarrow y_2^*=\frac32=1.5$.
Check: dual value $12(1.5)+18(1)=18+18=36$ ✓.

#### 3.2 Perturbing $b_2$

Change $b=(4,12,18)$ to $\bar b=(4,\ 12+\Delta,\ 18)$.
- Original optimal value: $Z_1=b_1y_1^*+b_2y_2^*+b_3y_3^*$.
- If $\Delta$ is small enough that the **same constraints stay tight** (the same basis stays optimal), the same $y^*$ stays dual-optimal, so
$$Z_2=b_1y_1^*+(b_2+\Delta)y_2^*+b_3y_3^*=Z_1+\Delta y_2^* .$$
- Hence
$$\boxed{y_2^*=\frac{Z_2-Z_1}{\Delta}}$$
the rate of change of the optimal value per unit of resource 2.

**Numerical check, $\Delta=1$** (constraint becomes $2x_2\le13$): the new optimum is where $2x_2=13$ meets $3x_1+2x_2=18$: $x_2=6.5$, $x_1=\frac{18-13}{3}=\frac53$. Value $3\cdot\frac53+5\cdot6.5=5+32.5=37.5$. Change $=37.5-36=1.5=y_2^*$ ✓.

🔍 *Validity range:* this holds while $(x_1,x_2)$ stays feasible: $x_1=\frac{18-(12+\Delta)}{3}\ge0\Rightarrow\Delta\le6$, and $x_1\le4\Rightarrow\Delta\ge-6$. Outside $[-6,6]$ the basis changes and so does the shadow price.

**Interpretation:** $y_1^*=0$: extra capacity in constraint 1 is worthless (it is not binding). $y_2^*=1.5$: one more unit of resource 2 is worth 1.5. $y_3^*=1$: one more unit of resource 3 is worth 1.

---

### 4. Summary

✅ Optimal value $Z=b^Ty^*$; for small changes of $b_i$, $\Delta Z=\Delta\cdot y_i^*$.
✅ Dual variable = shadow price = marginal value of a resource.
✅ Non-binding constraint ⇒ shadow price 0 (complementary slackness).

<!-- LECTURE END: sensitivity-lp -->

---

<!-- LECTURE START: co-optimality -->
<a id="co-optimality"></a>

## Lecture 33: Constrained Optimality (Part 1): feasible & descent directions, active constraints, toward Lagrange

*Source file:* `Co-optimality.pdf`

### 1. What this lecture is about

How do we recognise a local minimum when there are constraints? We can no longer demand $\nabla f=0$ (the minimum may sit on the boundary). This lecture builds the ideas:
1. A generic descent algorithm and local minima in a constrained set
2. **Descent directions** and **feasible directions**. At a local min, no direction is both.
3. **Active constraints**
4. Interior local min ⇒ $\nabla f=0$
5. With an equality constraint $h(x)=0$: move along directions with $\nabla h^Td=0$, the start of the Lagrange condition

---

### 2. Background you need first

- Descent direction: $\nabla f(x)^Td<0$ (Descent Methods block).
- Ball $B(x^*,\delta)$; local minimum.
- Gradients are perpendicular to level sets (Calculus Review).

---

### 3. Main content with full derivations

#### 3.1 Generic algorithm and local minima

For $\min f(x)$ s.t. $x\in\Omega$:
```
Initialise x^0, k = 0
do { x^{k+1} = x^k + α d }  while ( f(x^{k+1}) < f(x^k) )
```
The method stops at points where no further improvement is possible: **local minima**. With a constraint like $a\le x\le b$, these include interior valleys *and* boundary points (e.g. at $x=a$ the function may be increasing into the interval, making $a$ a local min).

**Local min (constrained):** $x^*$ such that $\exists\delta>0$ with
$$f(x^*)\le f(x)\quad\forall x\in B(x^*,\delta)\cap\Omega .$$
(Only feasible neighbours count.)

#### 3.2 Descent and feasible directions

- **Descent direction** $d$ at $x$: $\exists r>0$ with $f(x+\alpha d)\le f(x)$ for all $\alpha\in(0,r)$. (Sufficient: $\nabla f(x)^Td<0$.)
- **Feasible direction** $d$ at $x$: $\exists r>0$ with $x+\alpha d\in\Omega$ for all $\alpha\in(0,r)$ (you can move a little along $d$ and stay feasible).

✅ **Necessary condition:** if $x^*$ is a local min, there is **no direction that is both feasible and (strictly) descent**.
🔍 *Why:* if such a $d$ existed, small steps along it would stay in $\Omega$, stay inside $B(x^*,\delta)$, and strictly decrease $f$, contradicting local minimality.

#### 3.3 Active constraints

For an inequality $g_1(x)\le0$: at a point $\bar x$, if $g_1(\bar x)=0$ we say **$g_1$ is active at $\bar x$** (the point is on that constraint's boundary). If $g_1(\bar x)<0$ it is **inactive**: small moves in any direction keep it satisfied.

**Example:**
$$\min(x_1-4)^2+(x_2-8)^2\ \text{ s.t. }\ g_1:x_1\le4,\ g_2:x_2\le6,\ g_3:3x_1+2x_2\le18,\ g_4:-x_1\le0,\ g_5:-x_2\le0 .$$
(Contours: circles around $(4,8)$, which lies **outside** the feasible pentagon.)
- At $(0,6)$: $g_4$ and $g_2$ active; $g_1,g_3$ inactive.
- At $(2,0)$: only $g_5$ active.
- At an interior point: none active.

🔍 *Solution of this example (not in notes):* we want the feasible point closest to $(4,8)$. Projecting $(4,8)$ onto the line $3x_1+2x_2=18$ gives $(\frac{22}{13},\frac{84}{13})\approx(1.69,6.46)$, which violates $x_2\le6$. Along the edge $x_2=6$, $(x_1-4)^2+4$ is smallest at the largest allowed $x_1$, which is $x_1=2$ (from $g_3$). So $x^*=(2,6)$, $f^*=4+4=8$, with **$g_2,g_3$ active**.

#### 3.4 Interior local minimum

If $x^*$ is an **interior point** (no constraint active), then **all directions are feasible**. So no direction can be descent, which forces
$$\nabla f(x^*)=0$$
(otherwise $d=-\nabla f(x^*)$ would be feasible and descent). This is just the unconstrained FONC.

#### 3.5 Local min with an active (equality) constraint

**Example:** $\min f(x)=x_1+x_2$ s.t. $h(x)=x_1^2+x_2^2-1=0$ (the unit circle).
Points: $A=(-\frac1{\sqrt2},\frac1{\sqrt2})$, $B=(\frac1{\sqrt2},\frac1{\sqrt2})$, $C=(\frac1{\sqrt2},-\frac1{\sqrt2})$, $D=(-\frac1{\sqrt2},-\frac1{\sqrt2})$.

- Descent: $\nabla f(x)^Td<0$.
- To stay on $h(x)=0$ (to first order), move so that $h$ does not change: **$\nabla h(x)^Td=0$**, i.e. $d$ tangent to the circle.

$\nabla f(x)=\binom11$, $\nabla h(x)=\binom{2x_1}{2x_2}$.
**At $A$:** $\nabla h(A)=(-\sqrt2,\sqrt2)$; tangent directions $d\propto(1,1)$ or $(-1,-1)$. $\nabla f^T(1,1)=2>0$ but $\nabla f^T(-1,-1)=-2<0$: **there is a feasible descent direction**, so $A$ is not a local min.
🔍 **At $D$:** $\nabla h(D)=(-\sqrt2,-\sqrt2)$; tangent $d\propto(1,-1)$: $\nabla f^Td=1-1=0$ for both $\pm(1,-1)$. No tangent direction decreases $f$ to first order, so $D$ is a candidate. Here $\nabla f(D)=(1,1)$ is **parallel** to $\nabla h(D)$: $\nabla f(D)=-\frac1{\sqrt2}\nabla h(D)$.

✅ **The pattern:** at a constrained optimum, $\nabla f$ is a multiple of $\nabla h$, i.e. $\nabla f+\lambda\nabla h=0$. This is the **Lagrange condition**, developed in the next lectures.

---

### 4. Summary

✅ Local min in $\Omega$: best among feasible points nearby.
✅ At a local min, no direction is both feasible and descent.
✅ Active constraint: $g_i(\bar x)=0$. Interior min ⇒ $\nabla f=0$.
✅ On $h(x)=0$, allowed directions satisfy $\nabla h^Td=0$; at an optimum $\nabla f\parallel\nabla h$.

<!-- LECTURE END: co-optimality -->

---

<!-- LECTURE START: co-oct26 -->
<a id="co-oct26"></a>

## Lecture 34: General Constrained Optimization: examples, special classes (LP, QP, QCQP), eliminating constraints

*Source file:* `Co-oct26.pdf`

### 1. What this lecture is about

We move from LPs to **general** constrained problems
$$\min f(x)\quad\text{s.t.}\quad g_i(x)\le0\ (i=1..m),\qquad h_j(x)=0\ (j=1..l).$$
This lecture covers:
1. Where such problems come from
2. Special classes: LP, QP, QCQP, and which are convex
3. A first trick, **eliminating constraints by substitution**, and its pitfalls (worked example $\min x_1+x_2$ on the unit circle)

---

### 2. Background you need first

- Convex optimization (Convexity blocks); PSD matrices.
- 1-D calculus: derivative of $\sqrt{1-t^2}$ is $\frac{-t}{\sqrt{1-t^2}}$ (chain rule).

---

### 3. Main content with full derivations

#### 3.1 Applications

1. Estimators in machine learning (e.g. constrained least squares, SVMs).
2. Designing a box with minimum surface area for a given volume $V$.
3. Design of electric circuits.

🔍 *Box example written out:* sides $x_1,x_2,x_3>0$: $\min\ 2(x_1x_2+x_2x_3+x_1x_3)$ s.t. $x_1x_2x_3=V$. (Answer: a cube, $x_i=V^{1/3}$.)

#### 3.2 Special classes

| Class | Form | Convex? |
|---|---|---|
| **LP** | $\min c^Tx$ s.t. $Ax\ge b$, $x\ge0$ | always (convex optimization) |
| **QP** | $\min x^TQx+c^Tx$ s.t. $Ax\ge b$, $x\ge0$, $Q\in\mathbb{R}^{n\times n}$ | only if $Q\succeq0$; **not necessarily** convex |
| **QCQP** (quadratically constrained QP) | $\min x^TQ_0x+c^Tx$ s.t. $x^TQ_ix+c_i^Tx+b_i\ge0$, $i=1..n$ | **not necessarily** convex |

🔍 For QCQP: a constraint $x^TQ_ix+\dots\ge0$ defines a convex set when $-Q_i\succeq0$ (i.e. it is "$\text{convex}\le0$" after negating). With $Q_0\succeq0$ as well, the QCQP is convex.

#### 3.3 Eliminating constraints by substitution

**Sign constraints:**
- $\min f(x_1,x_2)$ s.t. $x_1\ge0$, $x_2\in\mathbb{R}$: substitute $x_1=x_3^2$ (any real $x_3$ gives $x_1\ge0$), giving an unconstrained problem in $(x_3,x_2)$.
- $x_1>0$: substitute $x_1=e^y$, $y\in\mathbb{R}$.

**Equality constraint example:**
$$\min x_1+x_2\quad\text{s.t.}\quad x_1^2+x_2^2=1 .$$
Solve the constraint for $x_1$: $x_1=\pm\sqrt{1-x_2^2}$. **Two branches**, which must be handled separately.

*Branch $x_1=+\sqrt{1-x_2^2}$:* minimise $\phi(x_2)=x_2+\sqrt{1-x_2^2}$ over $x_2\in[-1,1]$.
🔍 $\phi'(x_2)=1-\frac{x_2}{\sqrt{1-x_2^2}}=0\Rightarrow\sqrt{1-x_2^2}=x_2$, which requires $x_2\ge0$; squaring gives $1-x_2^2=x_2^2\Rightarrow x_2=\frac1{\sqrt2}$. (The notes list $\pm\frac1{\sqrt2}$; $-\frac1{\sqrt2}$ is an extraneous root introduced by squaring.) Point $(\frac1{\sqrt2},\frac1{\sqrt2})$, value $\sqrt2$, a **maximum**. Endpoints: $x_2=\pm1$ give $\pm1$.

*Branch $x_1=-\sqrt{1-x_2^2}$:* minimise $\psi(x_2)=x_2-\sqrt{1-x_2^2}$.
🔍 $\psi'=1+\frac{x_2}{\sqrt{1-x_2^2}}=0\Rightarrow\sqrt{1-x_2^2}=-x_2$, needs $x_2\le0$, so $x_2=-\frac1{\sqrt2}$, $x_1=-\frac1{\sqrt2}$. Value $-\sqrt2$.

**Answer:** $x^*=\big(-\frac1{\sqrt2},-\frac1{\sqrt2}\big)$, $f^*=-\sqrt2$.
The four candidate points $(\pm\frac1{\sqrt2},\pm\frac1{\sqrt2})$ in the notes are where the lines $x_1+x_2=k$ and the axes-diagonals meet the circle. Geometrically, slide the line $x_1+x_2=k$ downward (toward $-\sqrt2$) until it last touches the circle, at $x_1+x_2=-\sqrt2$.

⚠️ **Lessons:** substitution can create several branches, square roots and extraneous solutions, and it gets messy in many dimensions. That motivates a systematic method: **Lagrange multipliers / KKT** (next lectures).

---

### 4. Summary

✅ General form: $\min f$ s.t. $g_i\le0$, $h_j=0$.
✅ LP is convex; QP and QCQP are convex only with PSD conditions.
✅ Substitutions ($x=t^2$, $x=e^y$, solving an equality) remove constraints but can be messy; on the circle the answer is $(-\frac1{\sqrt2},-\frac1{\sqrt2})$.

<!-- LECTURE END: co-oct26 -->

---

<!-- LECTURE START: oct28 -->
<a id="oct28"></a>

## Lecture 35: Constrained Optimality (Part 2): Lagrange condition for equality constraints, derived via Farkas

*Source file:* `Oct28.pdf`

### 1. What this lecture is about

1. Formal definitions: feasible direction, descent direction, active (binding) constraints
2. Analysis of $\min x_1+x_2$ on the circle at the four points $A,B,C,D$
3. Why "feasible" means $\nabla h^Td=0$ (to first order)
4. **Lagrange condition:** at an optimum with equality constraints $h_1=\dots=h_l=0$,
$$\nabla f(x^*)=\sum_{i=1}^l\mu_i\nabla h_i(x^*),\qquad\mu\in\mathbb{R}^l,$$
proved with **Farkas' lemma**
5. A first look at inequality constraints

---

### 2. Background you need first

- Constrained Optimality Part 1 (descent/feasible directions, active constraints).
- **Farkas:** exactly one of $\{z:Az=b,\ z\ge0\}$ and $\{y:A^Ty\le0,\ b^Ty>0\}$ is non-empty (FL block).
- First-order Taylor: $h(x+\alpha d)=h(x)+\alpha\nabla h(x)^Td+o(\|\alpha d\|)$.

---

### 3. Main content with full derivations

#### 3.1 Setting and definitions

$$\min f(x)\quad\text{s.t.}\quad g_i(x)\le0\ \ \forall i\in I,\qquad h_j(x)=0\ \ \forall j\in E\qquad(\text{together: }\Omega).$$
- $x^*$ is a **local min** if $\exists\delta>0$: $f(x^*)\le f(x)$ $\forall x\in B(x^*,\delta)\cap\Omega$.
- $d$ is a **feasible direction** at $x\in\Omega$ if $\exists\delta>0$: $x+\alpha d\in\Omega$ $\forall\alpha\in(0,\delta)$.
- $d$ is a **descent direction** if $\exists\delta>0$: $f(x+\alpha d)<f(x)$ $\forall\alpha\in(0,\delta)$.
- ✅ $x^*$ local min ⇒ **no $d$ is both feasible and descent.**
- **Active / binding constraints** at $\bar x$: $\{i\in I:g_i(\bar x)=0\}\cup E$. (Equality constraints are always active.)
- If the optimizer is in the **interior** ($g_i(x^*)<0$ for all $i$, no equalities), then $\nabla f(x^*)=0$, the same as unconstrained.

#### 3.2 The circle example, point by point

$$\min x_1+x_2\ \text{ s.t. }\ h(x)=x_1^2+x_2^2-1=0;\qquad\nabla f=\binom11,\ \nabla h=\binom{2x_1}{2x_2}.$$
- To remain feasible (first order): $\{d:\nabla h(x)^Td=0\}$. To descend: $\{d:\nabla f(x)^Td<0\}$.

**A $=(-\frac1{\sqrt2},\frac1{\sqrt2})$:** $\nabla h=(-\sqrt2,\sqrt2)$. $d=(-1,-1)$ has $\nabla h^Td=\sqrt2-\sqrt2=0$ and $\nabla f^Td=-2<0$: both feasible and descent ⇒ **A is not optimal**.
**B $=(\frac1{\sqrt2},\frac1{\sqrt2})$:** $\nabla h=(\sqrt2,\sqrt2)$, parallel to $\nabla f$. Any $d$ with $\nabla h^Td=0$ also has $\nabla f^Td=0$: no feasible descent direction to first order. **There is no $d$ with $\nabla f^Td<0$ and $\nabla h^Td=0$.** (B is actually the *maximum*: the first-order test cannot tell max from min.)
**C $=(\frac1{\sqrt2},-\frac1{\sqrt2})$:** $\nabla h=(\sqrt2,-\sqrt2)$; tangent $d=(-1,-1)$: $\nabla f^Td=-2<0$ ⇒ feasible descent exists ⇒ **C is not optimal**.
**D $=(-\frac1{\sqrt2},-\frac1{\sqrt2})$:** $\nabla h=(-\sqrt2,-\sqrt2)$, anti-parallel to $\nabla f$. **No $d$ is feasible and descent.** (The minimum.)

**Observation:** at B and D, $\nabla f(x)=\mu\nabla h(x)$ for some $\mu\in\mathbb{R}$. In fact
$$\{d:\nabla f^Td<0,\ \nabla h^Td=0\}=\emptyset\iff\exists\mu\in\mathbb{R}:\ \nabla f(x)=\mu\nabla h(x).$$
🔍 *Why (one constraint):* decompose any $d=d_1+d_2$ with $d_1=\beta\nabla h$ and $d_2\perp\nabla h$. If $\nabla f$ were not parallel to $\nabla h$, it would have a non-zero component $p$ perpendicular to $\nabla h$; then $d=-p$ satisfies $\nabla h^Td=0$ and $\nabla f^Td=-\|p\|^2<0$.

#### 3.3 Why feasibility ⇒ $\nabla h^Td=0$

If $x$ and $x+\alpha d$ both satisfy $h=0$ for all small $\alpha$:
$$h(x+\alpha d)=h(x)+\alpha\nabla h(x)^Td+q(\alpha,d),\qquad\lim_{\|\alpha d\|\to0}\frac{q(\alpha,d)}{\|\alpha d\|}=0 .$$
Subtract $h(x)=0$, use $h(x+\alpha d)=0$, divide by $\alpha$, let $\alpha\to0$:
$$0=\lim_{\alpha\to0}\frac{h(x+\alpha d)-h(x)}{\alpha}=\nabla h(x)^Td+\lim_{\alpha\to0}\frac{q(\alpha,d)}{\alpha}=\nabla h(x)^Td+0 .$$
(E.g. on the circle, $(x_1+\alpha d_1)^2+(x_2+\alpha d_2)^2=1$ to first order gives $2x_1d_1+2x_2d_2=0$.) On a curved constraint, moving exactly along $d$ leaves the curve slightly, so in general "feasible" is taken to first order: the **tangent directions**.

#### 3.4 Lagrange condition for several equality constraints

$$\min f(x)\ \text{ s.t. }\ h_1(x)=0,\ \dots,\ h_l(x)=0 .$$
Let $x^*$ be an optimum. Then no first-order feasible descent direction exists:
$$\{d:\ \nabla f(x^*)^Td<0,\ \nabla h_1(x^*)^Td=0,\dots,\nabla h_l(x^*)^Td=0\}=\emptyset .$$
**Rewrite each equality as two inequalities:**
$$\{d:\ \nabla f(x^*)^Td<0,\ \ \nabla h_i(x^*)^Td\ge0,\ \ -\nabla h_i(x^*)^Td\ge0\ \ \forall i\}=\emptyset .$$
**Match it to Farkas.** Define
$$A=\big[\nabla h_1(x^*)\ \cdots\ \nabla h_l(x^*)\ \ -\nabla h_1(x^*)\ \cdots\ -\nabla h_l(x^*)\big]\in\mathbb{R}^{n\times2l},\qquad b=\nabla f(x^*).$$
The empty set is $\{d:A^Td\ge0,\ b^Td<0\}$. Put $y=-d$: $\{y:A^Ty\le0,\ b^Ty>0\}=\emptyset$, which is Farkas system 2. **By Farkas, system 1 is non-empty:**
$$\exists z\in\mathbb{R}^{2l},\ z\ge0:\quad Az=b .$$
Write it out:
$$\nabla f(x^*)=\sum_{i=1}^lz_i\nabla h_i(x^*)+\sum_{i=1}^lz_{l+i}\big(-\nabla h_i(x^*)\big)=\sum_{i=1}^l(z_i-z_{l+i})\nabla h_i(x^*).$$
Let $\mu_i=z_i-z_{l+i}$ (any sign, since it is a difference of non-negatives):
$$\boxed{x^*\text{ optimum}\ \Rightarrow\ \exists\mu\in\mathbb{R}^l:\ \nabla f(x^*)=\sum_{i=1}^l\mu_i\nabla h_i(x^*)}$$
This is the **Lagrange multiplier condition**. (Rigorous versions need a *regularity* assumption, e.g. the $\nabla h_i(x^*)$ linearly independent, so that the tangent directions are really the feasible ones.)

**Check on the circle:** at $D$, $\nabla f=(1,1)$, $\nabla h=(-\sqrt2,-\sqrt2)$, so $\mu=-\frac1{\sqrt2}$ ✓.

#### 3.5 Inequality constraints (preview)

$$\min f(x)\ \text{ s.t. }\ g(x)\le0 .$$
Example: $\min x_1+x_2$ s.t. $x_1^2+x_2^2\le1$ (the whole disc). Descent needs $\nabla f(x)^Td<0$. **What must $\nabla g(x)^Td$ satisfy?** If $g$ is inactive (interior), anything. If active, we need $g$ not to increase: $\nabla g(x)^Td\le0$. This one-sided condition leads to **sign-restricted multipliers** ($\lambda\ge0$), the **KKT conditions** (next lectures).

---

### 4. Summary

✅ At a local min: no direction is both feasible and descent.
✅ For $h(x)=0$, first-order feasible directions satisfy $\nabla h^Td=0$.
✅ $\{d:\nabla f^Td<0,\ \nabla h_i^Td=0\ \forall i\}=\emptyset$ ⟺ (Farkas) $\nabla f(x^*)=\sum_i\mu_i\nabla h_i(x^*)$, $\mu_i\in\mathbb{R}$.
✅ Inequalities give one-sided conditions, leading to KKT with $\lambda\ge0$.

<!-- LECTURE END: oct28 -->

---

<!-- LECTURE START: nov2 -->
<a id="nov2"></a>

## Lecture 36: KKT Conditions: statement, regular points, worked examples (and why regularity matters)

*Source file:* `Nov2.pdf`

### 1. What this lecture is about

The **Karush–Kuhn–Tucker (KKT)** conditions are the first-order necessary conditions for general constrained problems (inequalities and equalities). They extend "$\nabla f=0$" and the Lagrange condition. We:
1. Define the **Lagrangian** and the KKT conditions
2. Define **regular points** (constraint qualification)
3. Solve an example by **case analysis on active constraints**
4. See an example where the optimum is **not** a KKT point because it is not regular

---

### 2. Background you need first

- Lagrange condition $\nabla f=\sum\mu_i\nabla h_i$ (previous block).
- Active constraints; linear independence; complementary slackness from LP.

---

### 3. Main content with full derivations

#### 3.1 Problem, Lagrangian, KKT

$$\min f(x)\quad\text{s.t.}\quad g_i(x)\le0\ \ i=1..m\ [\lambda_i],\qquad h_j(x)=0\ \ j=1..l\ [\mu_j],\qquad x\in\mathbb{R}^n .$$
**Lagrangian:**
$$L(x,\lambda,\mu)=f(x)+\sum_{i=1}^m\lambda_ig_i(x)+\sum_{j=1}^l\mu_jh_j(x).$$
$\lambda_i,\mu_j$ are **Lagrange multipliers**.

> **First-order necessary conditions (KKT).** If $x^*$ is a local minimum **and $x^*$ is regular**, then $\exists\lambda\in\mathbb{R}^m$, $\mu\in\mathbb{R}^l$ such that $(x^*,\lambda,\mu)$ is a **KKT point**:
> 1. **Stationarity:** $\nabla_xL(x^*,\lambda,\mu)=\nabla f(x^*)+\sum_i\lambda_i\nabla g_i(x^*)+\sum_j\mu_j\nabla h_j(x^*)=0$
> 2. **Complementary slackness (CS):** $\lambda_ig_i(x^*)=0\ \ \forall i$
> 3. **Primal feasibility:** $g_i(x^*)\le0$, $h_j(x^*)=0$
> 4. **Dual feasibility:** $\lambda_i\ge0\ \ \forall i$ ($\mu_j$ any sign)

🔍 *Meaning of each:*
- CS: if constraint $i$ is **inactive** ($g_i<0$) then $\lambda_i=0$: it plays no role. Only active constraints can "push".
- $\lambda_i\ge0$: an active inequality can only push $x$ **back into** the feasible side. At the optimum, $-\nabla f$ (the downhill direction) points *out* of the feasible set, so it must be a non-negative combination of the outward normals $\nabla g_i$.
- Equality constraints push both ways, so $\mu_j$ is free (matches the Lagrange condition with sign flipped: $\nabla f=-\sum\mu_j\nabla h_j$).

#### 3.2 Regular points

$\bar x$ is **regular** if:
- all interior points of $\Omega$ are regular; and
- the gradients of the active constraints $\{\nabla g_i(\bar x),\nabla h_j(\bar x):\ i,j\in\mathcal A(\bar x)\}$ are **linearly independent** (LICQ). For a problem without equalities: $\{\nabla g_i(\bar x):i\in\mathcal A(\bar x)\}$ linearly independent.

Regularity ensures the linearised (tangent) feasible directions describe the true feasible set near $\bar x$, which the Farkas argument needs.

#### 3.3 Example 1 (case analysis)

$$\min x_1+x_2\quad\text{s.t.}\quad x_1^2+x_2^2\le1,\quad x_2\ge0 .$$
Standard form: $g_1=x_1^2+x_2^2-1\le0\ [\lambda_1]$, $g_2=-x_2\le0\ [\lambda_2]$. (Upper half-disc.)
$$L=x_1+x_2+\lambda_1(x_1^2+x_2^2-1)-\lambda_2x_2,\qquad\nabla_xL=\begin{pmatrix}1+2\lambda_1x_1\\1+2\lambda_1x_2-\lambda_2\end{pmatrix}=\begin{pmatrix}0\\0\end{pmatrix}.$$
With 2 inequality constraints there are $2^2=4$ cases.

**Case 1: both inactive** ($x_1^2+x_2^2<1$, $x_2>0$). CS ⇒ $\lambda_1=\lambda_2=0$ ⇒ $\nabla_xL=(1,1)\neq(0,0)$. ✗ No KKT point.

**Case 2: both active:** $x_1^2+x_2^2=1$, $x_2=0$ ⇒ $x_1=\pm1$: points $(1,0)$ and $(-1,0)$.
- $(1,0)$: $1+2\lambda_1=0\Rightarrow\lambda_1=-\frac12<0$ ✗ (dual infeasible).
- $(-1,0)$: $1-2\lambda_1=0\Rightarrow\lambda_1=\frac12$; second row $1+0-\lambda_2=0\Rightarrow\lambda_2=1$. Both $\ge0$ ✓.
**$\big((-1,0),\ \lambda_1=\frac12,\ \lambda_2=1\big)$ is a KKT point.**

**Case 3: (1) inactive, (2) active:** $x_1^2+x_2^2<1$, $x_2=0$, $\lambda_1=0$. Then row 1: $1+0=0$ ✗.

**Case 4: (1) active, (2) inactive:** $x_1^2+x_2^2=1$, $x_2>0$, $\lambda_2=0$. Rows: $1+2\lambda_1x_1=0$, $1+2\lambda_1x_2=0$ ⇒ $x_1=x_2$ and $\lambda_1=-\frac1{2x_2}<0$ (since $x_2>0$) ✗.

**Conclusion:** the only KKT point is $x^*=(-1,0)$, value $-1$. (Picture: the line $x_1+x_2=k$ slides down-left and last touches the half-disc at $(-1,0)$.) Both active gradients there, $\nabla g_1=(-2,0)$ and $\nabla g_2=(0,-1)$, are independent, so it is regular.

🔍 *Check the sign intuition:* $-\nabla f=(-1,-1)=\lambda_1\nabla g_1+\lambda_2\nabla g_2=\frac12(-2,0)+1\cdot(0,-1)$ ✓.

#### 3.4 Example 2: an optimum that is not a KKT point

$$\min x_1+x_2\quad\text{s.t.}\quad h_1=(x_1-1)^2+x_2^2-1=0,\quad h_2=(x_1+1)^2+x_2^2-1=0 .$$
Two unit circles centred at $(\pm1,0)$ touch only at the origin, so **the feasible set is $\{(0,0)\}$** and $x^*=(0,0)$ is (trivially) optimal.
- $\nabla h_1=\binom{2(x_1-1)}{2x_2}\to\binom{-2}{0}$, $\nabla h_2=\binom{2(x_1+1)}{2x_2}\to\binom20$ at the origin: **linearly dependent** ⇒ **not regular**.
- Linearised feasible directions: $\{d:\nabla h_1^Td=0,\ \nabla h_2^Td=0\}=\{d:d_1=0\}$, the whole vertical line. But the true feasible set is a single point.
- KKT would need $(1,1)+\mu_1(-2,0)+\mu_2(2,0)=0$; second component $1=0$ is impossible. **No multipliers exist.**

✅ Lesson: the KKT theorem **requires regularity**; without it an optimum may fail KKT.

---

### 4. Summary

✅ $L=f+\sum\lambda_ig_i+\sum\mu_jh_j$.
✅ KKT: $\nabla_xL=0$; $\lambda_ig_i=0$; feasibility; $\lambda\ge0$.
✅ Regular (LICQ): active constraint gradients linearly independent.
✅ Solve by enumerating active sets; discard cases with $\lambda<0$ or infeasibility.
✅ Non-regular optima may have no KKT multipliers.

<!-- LECTURE END: nov2 -->

---

<!-- LECTURE START: nov4 -->
<a id="nov4"></a>

## Lecture 37: KKT Examples: KKT points need not be minimizers, and minimizers need not be KKT points

*Source file:* `Nov4.pdf`

### 1. What this lecture is about

Three worked KKT problems that show the limits of the KKT conditions:
1. $\min x_1^2-x_2^2$ on the unit disc: one KKT point is a **saddle**, not a minimum
2. $\min x_2$ s.t. $x_2\ge-x_1^2$, $-1\le x_1\le1$: three KKT points, only two are minima
3. $\min x$ s.t. $x^2+y^2\le0$: the minimizer is **not** a KKT point (it is not regular)

---

### 2. Background you need first

- KKT conditions and regularity (previous block).
- Case analysis on which constraints are active.

---

### 3. Main content with full derivations

#### 3.1 Example 1

$$\min x_1^2-x_2^2\quad\text{s.t.}\quad x_1^2+x_2^2\le1\ [\lambda].$$
$$L(x,\lambda)=x_1^2-x_2^2+\lambda(x_1^2+x_2^2-1),\qquad\nabla_xL=\begin{pmatrix}2x_1+2\lambda x_1\\-2x_2+2\lambda x_2\end{pmatrix}=0\iff\begin{cases}2x_1(1+\lambda)=0\\2x_2(-1+\lambda)=0\end{cases}$$

**Constraint inactive** ($x_1^2+x_2^2<1$): $\lambda=0$ ⇒ $2x_1=0$, $-2x_2=0$ ⇒ $x=(0,0)$. KKT point $((0,0),0)$, $f=0$.

**Constraint active** ($x_1^2+x_2^2=1$), need $\lambda\ge0$:
- From row 1: since $\lambda\ge0$, $1+\lambda>0$, so $x_1=0$. Then $x_2=\pm1$.
- Row 2 with $x_2\neq0$: $\lambda=1$ ✓.
KKT points $((0,1),1)$ and $((0,-1),1)$, both with $f=-1$.
(🔍 The other option $x_2=0$, $x_1=\pm1$ would need $\lambda=-1<0$ ✗.)

**Which are minima?** $(0,\pm1)$ with $f=-1$ are the global minima. **$(0,0)$ is a KKT point but a saddle:** moving along $x_2$ decreases $f=-x_2^2$. ✅ **A KKT point need not be a minimizer.**

#### 3.2 Example 2

$$\min x_2\quad\text{s.t.}\quad x_2\ge-x_1^2,\quad-1\le x_1\le1 .$$
Standard form:
$$g_1=-x_1^2-x_2\le0\ [\lambda_1],\quad g_2=-x_1-1\le0\ [\lambda_2],\quad g_3=x_1-1\le0\ [\lambda_3].$$
(Feasible set: the region above the downward parabola $x_2=-x_1^2$, between the vertical lines $x_1=\pm1$.)
$\nabla g_1=\binom{-2x_1}{-1}$, $\nabla g_2=\binom{-1}{0}$, $\nabla g_3=\binom10$.
$$L=x_2+\lambda_1(-x_1^2-x_2)+\lambda_2(-x_1-1)+\lambda_3(x_1-1),\qquad\nabla_xL=\begin{pmatrix}-2\lambda_1x_1-\lambda_2+\lambda_3\\1-\lambda_1\end{pmatrix}=0 .$$
- Row 2: $\lambda_1=1>0$ ⇒ by CS, $g_1$ active: $\boxed{x_2=-x_1^2}$.

**Both (2) and (3) inactive** ($-1<x_1<1$): $\lambda_2=\lambda_3=0$; row 1: $-2x_1=0\Rightarrow x_1=0$, $x_2=0$. KKT point $(0,0)$, $f=0$.

**(2) active, (3) inactive:** $x_1=-1$, $\lambda_3=0$, $x_2=-1$. Row 1: $-2(1)(-1)-\lambda_2=0\Rightarrow\lambda_2=2\ge0$ ✓. KKT point $(-1,-1)$, $f=-1$.

**(3) active, (2) inactive:** $x_1=1$, $\lambda_2=0$, $x_2=-1$. Row 1: $-2(1)(1)+\lambda_3=0\Rightarrow\lambda_3=2\ge0$ ✓. KKT point $(1,-1)$, $f=-1$.

(Both (2),(3) active is impossible: $x_1$ can't be $-1$ and $1$.)

**Minima:** $(\pm1,-1)$ with $f=-1$. The KKT point $(0,0)$ (top of the parabola) is **not** a minimum. ✅ Again: KKT points need not be minimizers.

#### 3.3 Example 3: a minimizer that is not a KKT point

$$\min x\quad\text{s.t.}\quad x^2+y^2\le0\ [\lambda].$$
The only feasible point is $(0,0)$, so it is the minimizer.
$$L=x+\lambda(x^2+y^2),\qquad\nabla L=\binom{1+2\lambda x}{2\lambda y}=0 .$$
- Inactive ($x^2+y^2<0$): impossible ✗.
- Active ($x=y=0$): row 1 becomes $1=0$, i.e. $\lambda=-\frac1{2x}$ is undefined. **No multiplier exists.** ??

**Why:** $\nabla g(0,0)=(2x,2y)=(0,0)$. A zero gradient is not linearly independent, so the point is **not regular**. The linearised feasible directions $\{d:\nabla g^Td\le0\}=\{d:0\le0\}=\mathbb{R}^2$ (every direction looks feasible), but in reality no movement is possible.

✅ **A minimizer need not be a KKT point** (when it is not regular).

---

### 4. Summary

✅ Solve KKT by enumerating active sets; keep cases with feasibility and $\lambda\ge0$.
✅ KKT is **necessary** (for regular points), **not sufficient**: saddles and maxima can be KKT points.
✅ Without regularity, even the true minimizer can fail KKT.
✅ (For **convex** problems KKT is also sufficient: next block.)

<!-- LECTURE END: nov4 -->

---

<!-- LECTURE START: nov4-convex -->
<a id="nov4-convex"></a>

## Lecture 38: KKT for Convex Problems: Slater's condition and KKT sufficiency

*Source file:* `Nov4-convex.pdf`

### 1. What this lecture is about

For **convex** problems the KKT story becomes clean:
1. **Slater's condition** (a strictly feasible point exists) replaces regularity: optimum ⇒ KKT point
2. **KKT is sufficient:** any KKT point of a convex problem is a global minimizer (full proof)

So for convex problems, **KKT ⟺ optimal** (given Slater).

---

### 2. Background you need first

- KKT conditions (KKT blocks).
- First-order convexity: $f$ convex ⇒ $f(x)\ge f(x^*)+\nabla f(x^*)^T(x-x^*)$ for all $x$ (Convexity Part 2).
- For a linear (affine) $h$: $h(x)=h(x^*)+\nabla h(x^*)^T(x-x^*)$ **exactly**.

---

### 3. Main content with full derivations

#### 3.1 Convex optimization problem

$$\min f(x)\quad\text{s.t.}\quad g_i(x)\le0\ (i=1..m),\qquad h_j(x)=0\ (j=1..l)\qquad(\text{set }\Omega),$$
where **$f,g_i$ are convex functions** and **$h_j$ are linear (affine)**. (Then $\Omega$ is convex: sublevel sets of convex functions are convex, hyperplanes are convex, and intersections of convex sets are convex.)

#### 3.2 Slater's condition

$$\exists\bar x:\quad g_i(\bar x)<0\ \ \forall i=1..m,\qquad h_j(\bar x)=0\ \ \forall j .$$
$\bar x$ is a **Slater point**: strictly inside all inequality constraints.

> **Theorem.** For a convex problem, if a Slater point exists and $x^*$ is optimal, then $\exists\lambda,\mu$ such that $(x^*,\lambda,\mu)$ is a KKT point.

(This is a constraint qualification like regularity, but much easier to check. Example 3 of the previous block, $x^2+y^2\le0$, has **no** Slater point, since nothing is strictly $<0$. That is exactly why KKT failed there.)

#### 3.3 KKT is sufficient for convex problems

> **Theorem.** If the problem is convex and $(x^*,\lambda,\mu)$ is a KKT point, then $x^*$ is a (global) minimizer.

**Proof, line by line.** Take any feasible $x\in\Omega$.
1. Convexity of $f$ (first-order condition):
$$f(x)\ge f(x^*)+\nabla f(x^*)^T(x-x^*).$$
2. Convexity of each $g_i$, multiplied by $\lambda_i\ge0$ (non-negative, so the inequality keeps its direction):
$$\lambda_ig_i(x)\ge\lambda_i\big[g_i(x^*)+\nabla g_i(x^*)^T(x-x^*)\big].$$
3. Linearity of each $h_j$, multiplied by $\mu_j$ (any sign is fine for an equality):
$$\mu_jh_j(x)=\mu_j\big[h_j(x^*)+\nabla h_j(x^*)^T(x-x^*)\big].$$
4. Add everything:
$$f(x)+\sum_i\lambda_ig_i(x)+\sum_j\mu_jh_j(x)\ \ge\ f(x^*)+\sum_i\underbrace{\lambda_ig_i(x^*)}_{=0\ (\text{CS})}+\sum_j\mu_j\underbrace{h_j(x^*)}_{=0}+\underbrace{\Big(\nabla f(x^*)+\sum_i\lambda_i\nabla g_i(x^*)+\sum_j\mu_j\nabla h_j(x^*)\Big)^T}_{=0\ (\text{stationarity})}(x-x^*).$$
So the right side is just $f(x^*)$.
5. On the left: $x$ feasible ⇒ $g_i(x)\le0$ and $\lambda_i\ge0$ ⇒ $\sum_i\lambda_ig_i(x)\le0$; and $h_j(x)=0$. Hence
$$f(x)\ \ge\ f(x)+\underbrace{\sum_i\lambda_ig_i(x)}_{\le0}+\underbrace{\sum_j\mu_jh_j(x)}_{=0}\ \ge\ f(x^*).$$
6. This holds for every feasible $x$ ⇒ **$x^*$ is optimal.** ∎

(Picture: for convex $f$ the tangent line at $x^*$ lies below the graph, so nothing feasible can be lower.)

---

### 4. Worked example

$\min x_1^2+x_2^2$ s.t. $x_1+x_2\ge1$, i.e. $g=1-x_1-x_2\le0$. Convex; Slater point $(1,1)$ ($g=-1<0$).
KKT: $\nabla L=(2x_1-\lambda,\ 2x_2-\lambda)=0\Rightarrow x_1=x_2=\lambda/2$. If $\lambda=0$, then $x=0$, infeasible ✗. So active: $x_1+x_2=1\Rightarrow\lambda=1$, $x^*=(\frac12,\frac12)$, $\lambda=1\ge0$ ✓. By sufficiency, $x^*$ is the global min, $f^*=\frac12$.

---

### 5. Summary

✅ Convex problem: $f,g_i$ convex, $h_j$ affine.
✅ Slater: a point with all $g_i<0$ (and $h=0$). Then optimum ⇒ KKT.
✅ KKT point of a convex problem ⇒ global minimum (tangent-plane inequalities + CS + stationarity).

<!-- LECTURE END: nov4-convex -->

---

<!-- LECTURE START: nov5 -->
<a id="nov5"></a>

## Lecture 39: Convex KKT Applications: entropy maximisation, maximum likelihood, LP via KKT, projected gradient descent

*Source file:* `Nov5.pdf`

> 📌 The first half of `Nov5.pdf` repeats Slater's condition and the proof that KKT is sufficient for convex problems. Those are explained in the block **"KKT for Convex Problems: Slater's condition and KKT sufficiency"**. This block covers the **applications** that are new on Nov 5.

### 1. What this lecture is about

1. **Entropy maximisation:** the uniform distribution has maximum entropy
2. **Maximum likelihood estimation (MLE)** of word probabilities: estimate = relative frequency
3. (HW) KKT conditions of an LP give back LP duality and complementary slackness
4. Constrained least squares and **projected gradient descent**

---

### 2. Background you need first

- KKT conditions; for convex problems with a Slater point, KKT ⟺ optimal.
- $\frac{d}{dp}(p\log p)=\log p+1$; $\frac{d}{dq}\log q=\frac1q$.
- Probability: a discrete random variable $X=x_i$ with probability $p_i$, $\sum p_i=1$, $p_i\ge0$.
- Indicator $\mathbb 1\{A\}=1$ if $A$ is true, else 0.

---

### 3. Main content with full derivations

#### 3.1 Entropy maximisation

$X=x_i$ with probability $p_i$, $i=1..n$. **Entropy** $H(X)=-\sum_{i=1}^np_i\log p_i$ measures the randomness of $X$.
$$\max_p\ -\sum_ip_i\log p_i\ \text{ s.t. }\sum_ip_i=1,\ p_i\ge0\qquad\Longleftrightarrow\qquad\min_p\ \sum_ip_i\log p_i\ \text{ s.t. }\sum_ip_i=1\ [\mu],\ p_i\ge0\ [\lambda_i].$$
**Convexity check:** $f(p)=\sum_ip_i\log p_i$,
$$\nabla f(p)=\begin{pmatrix}1+\log p_1\\\vdots\\1+\log p_n\end{pmatrix},\qquad\nabla^2f(p)=\operatorname{diag}\Big(\frac1{p_1},\dots,\frac1{p_n}\Big)\succ0\ (p_i>0).$$
So $f$ is convex, and the feasible set $\Omega$ is convex (a linear equality plus half-spaces). **Slater:** $p_i=\frac1n>0$ works. So the minimiser satisfies KKT.

**Solve** (the solution turns out to have $p_i>0$, so $\lambda_i=0$ by CS; we drop them):
$$\mathcal L(p,\mu)=\sum_ip_i\log p_i+\mu\Big(\sum_ip_i-1\Big),\qquad\frac{\partial\mathcal L}{\partial p_i}=1+\log p_i+\mu=0\quad\forall i .$$
So $\log p_i=-1-\mu$ for **every** $i$: all $p_i$ are equal, $p_i=\alpha$. From $\sum_ip_i=n\alpha=1$: $\alpha=\frac1n$.
$$\boxed{p_i^*=\frac1n\ \ \forall i}\quad(\ge0\ ✓)$$
**The uniform distribution maximises entropy** (maximum entropy $\log n$).

#### 3.2 Maximum likelihood estimation of word probabilities

Words $w_1,\dots,w_n$ with unknown probabilities $q_1,\dots,q_n$. A dataset $\mathcal D=\{x_1,\dots,x_m\}$ of $m$ words is generated **i.i.d.** (independent, identically distributed). Estimate $q$.

**Likelihood** (probability of seeing the data given $q$):
$$P(\mathcal D\mid q)=P(x_1\mid q)\cdots P(x_m\mid q)=\prod_{i=1}^m\prod_{j=1}^nq_j^{\mathbb 1\{x_i=w_j\}} .$$
🔍 For each data point $x_i$, exactly one indicator is 1 (the word it equals), so $\prod_jq_j^{\mathbb1\{x_i=w_j\}}=q_{\text{that word}}$.

**Log-likelihood** (log turns products into sums; $\log a^b=b\log a$):
$$\log P(\mathcal D\mid q)=\sum_{i=1}^m\sum_{j=1}^nz_{ij}\log q_j,\qquad z_{ij}=\mathbb1\{x_i=w_j\}.$$
**MLE problem:**
$$\max\ \sum_{i,j}z_{ij}\log q_j\ \text{ s.t. }\sum_jq_j=1,\ q_j\ge0\ \ \Longleftrightarrow\ \ -\min\Big[-\sum_{i,j}z_{ij}\log q_j\Big]\ \text{ s.t. }\sum_jq_j=1\ [\mu],\ q_j\ge0 .$$
($-\log$ is convex, so this is a convex problem.)
$$\mathcal L(q,\mu)=-\sum_{i,j}z_{ij}\log q_j+\mu\Big(\sum_jq_j-1\Big),\qquad\frac{\partial\mathcal L}{\partial q_j}=-\sum_{i=1}^m\frac{z_{ij}}{q_j}+\mu=0\ \ \forall j .$$
So $q_j=\frac{\sum_iz_{ij}}{\mu}$. Use $\sum_jq_j=1$: $\mu=\sum_j\sum_iz_{ij}=\sum_{i=1}^m1=m$ (each data point counts once). Therefore
$$\boxed{q_j^*=\frac{\sum_iz_{ij}}{m}=\frac{\#\text{ occurrences of word }j}{m}}$$
The MLE is the **relative frequency**.

#### 3.3 HW: KKT for an LP

$$\min c^Tx\ \text{ s.t. }\ Ax\ge b\ (a_i^Tx\ge b_i,\ i=1..m),\quad x\ge0\ (x_j\ge0,\ j=1..n).$$
🔍 *Solution.* Write as $b-Ax\le0\ [y\ge0]$, $-x\le0\ [s\ge0]$. Lagrangian $L=c^Tx+y^T(b-Ax)-s^Tx$.
- Stationarity: $c-A^Ty-s=0\Rightarrow s=c-A^Ty\ge0$, i.e. **$A^Ty\le c$** (dual feasibility!), $y\ge0$.
- CS: $y_i(b-Ax)_i=0$ and $s_jx_j=(c-A^Ty)_jx_j=0$: exactly **LP complementary slackness**.
- LPs are convex (and linear constraints need no Slater), so KKT ⟺ optimal: **LP duality is KKT**.

For $\min c^Tx$, $Ax=b$, $x\ge0$: the same steps with a free $y$ give $A^Ty\le c$, $(c-A^Ty)_jx_j=0$.

#### 3.4 Constrained least squares and projected gradient descent

$$\min_x\ \|b-Ax\|^2\ \text{ s.t. }\ x_i\ge0\ \ \forall i\qquad(A\in\mathbb{R}^{m\times n}).$$
A general method for $\min f(x)$ s.t. $x\in\Omega$ (with $\Omega$ convex) is **projected gradient descent**:
```
init x^(0), k = 0
until convergence:  x^(k+1) = Proj_Ω( x^(k) − α ∇f(x^(k)) );  k = k+1
```
Take an ordinary gradient step, then **project back** onto $\Omega$:
$$\operatorname{Proj}_\Omega(z)=\arg\min_{x\in\Omega}\|z-x\|^2\quad(\text{closest feasible point}).$$

**Example: projecting onto $\{x\ge0\}$.** $\min\|z-x\|^2$ s.t. $x_i\ge0$. It separates by coordinate: $\min(z_i-x_i)^2$ s.t. $x_i\ge0$.
🔍 *KKT:* $L=\sum_i(z_i-x_i)^2-\sum_i\lambda_ix_i$; $\partial_{x_i}:\ -2(z_i-x_i)-\lambda_i=0$.
- If $x_i>0$: $\lambda_i=0\Rightarrow x_i=z_i$ (valid when $z_i>0$).
- If $x_i=0$: $\lambda_i=-2z_i\ge0\Rightarrow z_i\le0$.

So $\boxed{\operatorname{Proj}(z)_i=\max(z_i,0)}$: just clip negatives to zero. For constrained least squares, $\nabla f=2A^T(Ax-b)$, so each step is $x\leftarrow\max\big(x-2\alpha A^T(Ax-b),\ 0\big)$.

---

### 4. Summary

✅ Max entropy ⇒ uniform $p_i=1/n$ (from $1+\log p_i+\mu=0$).
✅ MLE of a categorical distribution ⇒ relative frequencies ($\mu=m$).
✅ KKT of an LP = dual feasibility + complementary slackness.
✅ Projected GD: gradient step, then nearest feasible point; onto $x\ge0$ that is $\max(\cdot,0)$.

<!-- LECTURE END: nov5 -->

---

<!-- LECTURE START: nov6 -->
<a id="nov6"></a>

## Lecture 40: Practical Algorithms for Constrained Optimization: projected gradient, feasible directions (SLP), Frank–Wolfe

*Source file:* `Nov6.pdf`

### 1. What this lecture is about

Three practical algorithms for $\min f(x)$ s.t. $x\in\Omega$ ($g_i(x)\le0$, $h_j(x)=0$):
1. **Projected gradient descent**, fully worked for non-negative least squares (projection derived with KKT)
2. **Feasible-direction method**: pick the best descent direction among feasible ones by solving a small LP (**SLP**, sequential linear programming)
3. **Frank–Wolfe (conditional gradient):** minimise the linearisation over $\Omega$, then move toward that point

---

### 2. Background you need first

- Projection $\operatorname{Proj}_\Omega(z)=\arg\min_{x\in\Omega}\|z-x\|^2$ (previous block).
- KKT and Slater for convex problems.
- First-order Taylor: $f(x_k+\alpha d)\approx f(x_k)+\alpha\nabla f(x_k)^Td$.
- Least squares: $\nabla\|Ax-b\|^2=2A^T(Ax-b)$; unconstrained minimiser $x_u^*=(A^TA)^{-1}A^Tb$ (Linear Regression block).

---

### 3. Main content with full derivations

#### 3.1 Projected gradient descent

```
Initialise x^(0), k = 0
Till convergence:  x^(k+1) = P_Ω( x^(k) − α_k ∇f(x^(k)) );  k = k+1
```
$P_\Omega(z)=\arg\min_{x\in\Omega}\|z-x\|^2$, the **projection operator** (nearest point of $\Omega$ to $z$).

**Example: non-negative least squares.**
$$\min\|Ax-b\|^2\ \text{ s.t. }\ x_i\ge0,\qquad A\in\mathbb{R}^{m\times n},\ b\in\mathbb{R}^m .$$
- $f(x)=(Ax-b)^T(Ax-b)=x^TA^TAx-2b^TAx+b^Tb$, $g_i(x)=-x_i$.
- $\nabla f(x)=2A^TAx-2A^Tb=2A^T(Ax-b)$.
- Step: $x^{k+1}=P_\Omega\big(x^k-\alpha\cdot2A^T(Ax^k-b)\big)$.

**Computing $P_\Omega(z)$ for $\Omega=\{x\ge0\}$:**
$$P_\Omega(z)=\arg\min\ \tfrac12\|z-x\|^2\ \text{ s.t. }x_i\ge0\ [\lambda_i].$$
A Slater point exists (e.g. $x=\mathbf 1$), so the optimum satisfies KKT.
$$\mathcal L(x,\lambda)=\tfrac12(z-x)^T(z-x)-\sum_i\lambda_ix_i=\tfrac12(z^Tz-2z^Tx+x^Tx)-\sum_i\lambda_ix_i .$$
$\nabla_x\mathcal L=-z+x-\lambda=0\ \Rightarrow\ \lambda_i=x_i-z_i\ \forall i$.
KKT system: $\lambda_i=x_i-z_i$, $\lambda_ix_i=0$, $x_i\ge0$, $\lambda_i\ge0$. Substituting: $(x_i-z_i)x_i=0$, so $x_i=z_i$ or $x_i=0$.
- $x_i=z_i$ ⇒ $\lambda_i=0$; feasible only if $z_i\ge0$.
- $x_i=0$ ⇒ $\lambda_i=-z_i\ge0$; possible only if $z_i\le0$.
$$\boxed{x_i=\begin{cases}z_i&z_i\ge0\ (\lambda_i=0)\\0&z_i<0\ (\lambda_i=-z_i)\end{cases}=\max(z_i,0)}$$
**Algorithm:** $x_i^{k+1}=\max\big((x^k-\alpha\nabla f_k)_i,\ 0\big)$ for all $i$.

**HW:** if $x^*$ minimises $\|Ax-b\|^2$ over $x\ge0$, is $x^*=P_\Omega(x_u^*)$, i.e. just clip the unconstrained solution?
🔍 *Answer: not in general.* Take $f(x)=(x-u)^TQ(x-u)$ (a least-squares objective with $A^TA=Q$) with $Q=\begin{pmatrix}1&0.9\\0.9&1\end{pmatrix}$, $u=(1,-1)$. Clipping gives $(1,0)$. But with $x_2=0$, minimising over $x_1$: $\frac{d}{dx_1}\big[(x_1-1)^2+1.8(x_1-1)+1\big]=2(x_1-1)+1.8=0\Rightarrow x_1=0.1$. KKT check for $x_2$: $\partial f/\partial x_2=2[0.9(0.1-1)+(0+1)]=0.38\ge0$ ✓. So $x^*=(0.1,0)\neq(1,0)$. Clipping is exact only when the variables don't interact (e.g. $A^TA$ diagonal).

#### 3.2 Feasible-direction method (SLP)

At $x$, with active inequality set $\mathcal A(x)$:
- **Descent directions:** $\{d:\nabla f(x)^Td<0\}$.
- **Feasible directions (first order):** $\{d:\nabla g_i(x)^Td\le0\ \forall i\in\mathcal A(x),\ \ \nabla h_j(x)^Td=0\ \forall j\in E\}$.

Since $f(x_k+\alpha d)\approx f(x_k)+\alpha\nabla f(x_k)^Td$, pick the feasible direction with the most negative slope:
$$d_k=\arg\min_d\ \nabla f(x_k)^Td\ \text{ s.t. }\ \nabla g_i(x_k)^Td\le0\ (i\in\mathcal A(x_k)),\ \ \nabla h_j(x_k)^Td=0\ (j\in E),\ \ \|d\|_p\le1 .$$
(The norm bound stops $d$ from being scaled to $-\infty$.)
```
x^0;  loop:
   solve for d_k (above)
   if ∇f(x_k)ᵀd_k < 0 :  x_{k+1} = x_k + α d_k ;  k = k+1
   else : return x_k        # no feasible descent direction ⇒ (first-order) optimal
```
If $p=1$ (or $\infty$), $\|d\|_p\le1$ can be written with linear inequalities, so **$d_k$ is obtained from an LP**. Solving a sequence of LPs = **Sequential Linear Programming (SLP)**.
🔍 *E.g. $p=\infty$:* $-1\le d_i\le1$ for all $i$. *$p=1$:* write $d=d^+-d^-$, $d^\pm\ge0$, $\sum(d_i^++d_i^-)\le1$.

#### 3.3 Frank–Wolfe method

$$s_k=\arg\min_{s\in\Omega}\ \nabla f_k^Ts,\qquad x_{k+1}=(1-\alpha_k)x_k+\alpha_ks_k=x_k+\alpha_k\underbrace{(s_k-x_k)}_{d_k}.$$
- Minimise the **linear approximation** of $f$ over the whole feasible set to get a target point $s_k$, then move part of the way toward it.
- Since $\nabla f_k^T(s-x_k)=\nabla f_k^Ts-\nabla f_k^Tx_k$ and the last term is a constant,
$$\arg\min_{s\in\Omega}\nabla f_k^Ts\ \equiv\ \arg\min_{s\in\Omega}\nabla f_k^T(s-x_k):$$
$s_k-x_k$ is the direction of steepest first-order decrease **within $\Omega$**.
- If $\Omega$ is convex and $\alpha_k\in[0,1]$, $x_{k+1}$ is a convex combination of feasible points, hence **feasible**: no projection needed.
- If $\Omega$ is a polyhedron, each step is an **LP**. (More in the Frank–Wolfe slides block.)

---

### 4. Summary

✅ Projected GD: $x\leftarrow P_\Omega(x-\alpha\nabla f)$; onto $x\ge0$, $P$ = clip at 0 (derived from KKT).
✅ Clipping the unconstrained least-squares solution is **not** the constrained solution in general.
✅ Feasible-direction / SLP: best first-order feasible descent direction via an LP; stop when none exists.
✅ Frank–Wolfe: $s_k=\arg\min_{s\in\Omega}\nabla f_k^Ts$, $x_{k+1}=(1-\alpha_k)x_k+\alpha_ks_k$; stays feasible.

<!-- LECTURE END: nov6 -->

---

<!-- LECTURE START: yinyuye-lecture12 -->
<a id="yinyuye-lecture12"></a>

## Lecture 41: Constrained Optimization Algorithms I (Yinyu Ye, Stanford MS&E211 #12): gradient projection, feasible directions, reduced gradient, Newton on KKT

*Source file:* `YinyuYe_lecture12.pdf`

### 1. What this lecture is about

These are borrowed slides (Yinyu Ye, Stanford) that survey algorithms for constrained problems. This block explains each one, with the steps the slides skip:
1. Overview of method families
2. **Approximate KKT** as a stopping rule
3. **Gradient projection** for $x\ge0$ and for $Ax=b$ (the projection matrix $P_A$)
4. Feasible-direction methods and **SLP**
5. **Reduced gradient** (eliminate basic variables)
6. **Newton's method on the KKT equations** (portfolio QP and spherical-constraint examples)
7. An application: the Wasserstein barycenter LP

---

### 2. Background you need first

- KKT conditions; projected gradient (Nov 5–6 blocks).
- Newton's method for root finding: $z_{k+1}=z_k-[\nabla g(z_k)]^{-1}g(z_k)$, the multivariate version of $x-f/f'$ (OneDim block).
- Matrix inverse; $(AA^T)^{-1}$ exists when $A$ has independent rows.

---

### 3. Main content with full derivations

#### 3.1 Families of methods

- **Gradient projection:** take a gradient step, project onto the feasible set (first order).
- **Feasible-direction:** search along a direction that is both feasible and descent (first or second order).
- **Reduced gradient:** like simplex: change non-basic variables, let the basic ones adjust (first order).
- **Newton on KKT:** solve the KKT equations with Newton's method (second order).
- **Lagrangian methods:** relax equality constraints into the objective and update multipliers (see Dual Ascent).
- **Barrier / interior-point:** keep iterates strictly inside while improving the objective.

#### 3.2 Approximate KKT conditions

For $\min f(x)$ s.t. $Ax=b$, $x\ge0$, the KKT conditions are (with $y$ for $Ax=b$, $r$ for $x\ge0$):
$$\nabla f(x)-A^Ty-r=0,\qquad Ax=b,\qquad(x,r)\ge0,\qquad x_jr_j=0\ \forall j .$$
In practice we look for an **$\varepsilon$-solution**: $(x,r)\ge0$ and
$$\|\nabla f(x)-A^Ty-r\|\le\varepsilon,\qquad\|Ax-b\|\le\varepsilon,\qquad x_jr_j\le\varepsilon\ \forall j .$$

#### 3.3 Gradient projection for $x\ge0$

$\min f(x)$ s.t. $x\ge0$. KKT: $\nabla f(x)-y=0$, $y\ge0$; $x\,.\!*\,y=0$ (componentwise product = CS).
Projected step (as derived in the Nov 6 block):
$$x^{k+1}=\max\{0,\ x^k-\beta^{-1}\nabla f(x^k)\},\qquad y^{k+1}=\max\{0,\nabla f(x^k)\}\ (\text{dual estimate}).$$
($\beta^{-1}$ = step size; $\beta$ is a Lipschitz constant of $\nabla f$.) MATLAB loop for the QP $\min\frac12x^TQx+c^Tx$, $x\ge0$: `g=Q*x+c; x=max(0,x-(1/beta)*g);`, then `y=max(0,g)` and check `norm(x.*y)`.
*Upper bounds* $x\le u$: project by clipping to $[0,u]$, i.e. $\min(u,\max(0,\cdot))$.

#### 3.4 Gradient projection for $Ax=b$

$\min f(x)$ s.t. $Ax=b$; KKT: $\nabla f(x)-A^Ty=0$.
$$x^{k+1}=x^k-\beta^{-1}\big(I-A^T(AA^T)^{-1}A\big)\nabla f(x^k),\qquad y^{k+1}=(AA^T)^{-1}A\nabla f(x^k).$$
🔍 *Where $P_A=I-A^T(AA^T)^{-1}A$ comes from:* to stay on $Ax=b$, the step $d$ must satisfy $Ad=0$. Project $g=\nabla f$ onto $\mathcal N(A)$: write $g=A^Ty+p$ with $Ap=0$. Multiply by $A$: $Ag=AA^Ty\Rightarrow y=(AA^T)^{-1}Ag$, so $p=g-A^Ty=(I-A^T(AA^T)^{-1}A)g=P_Ag$. Check $AP_A=A-AA^T(AA^T)^{-1}A=0$ ✓. Also $P_A^2=P_A$ (projecting twice changes nothing).
If the start is feasible ($Ax^0=b$), all iterates stay feasible.

**Convergence:** if $\nabla f$ is $\beta$-Lipschitz and the step is $\beta^{-1}$, then within $\frac{2\beta(f(x^0)-f(x^*))}{\varepsilon^2}$ steps, $\|\nabla f(x^k)-A^Ty^k\|\le\varepsilon$ with $Ax^k=b$.

**Example.** $\min0.5(x_1-1)^2+0.5(x_2-1)^2$ s.t. $x_1+x_2=1$. $A=(1\ 1)$, start $x=(1,0)$, $\nabla f=(x_1-1,x_2-1)=(0,-1)$.
$$P=I-\tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix}=\begin{pmatrix}0.5&-0.5\\-0.5&0.5\end{pmatrix},\qquad-P\nabla f=-\begin{pmatrix}0.5\\-0.5\end{pmatrix}=\begin{pmatrix}-0.5\\0.5\end{pmatrix}.$$
A full step gives $(0.5,0.5)$, the minimiser, in one step (line search finds it).

#### 3.5 Application: Wasserstein barycenter

The **Wasserstein distance** between a supply distribution $s$ and a demand distribution $d$ is the minimum transportation cost:
$$WD(s,d)=\min\sum_{i,j}c_{ij}x_{ij}\ \text{ s.t. }\sum_jx_{ij}=s_i,\ \sum_ix_{ij}=d_j,\ x_{ij}\ge0 .$$
**Barycenter:** find $s$ minimising $\sum_kWD(s,d^k)$ over a set of distributions $d^k$, with a total-mass constraint (e.g. $s_1+\dots+s_4=9$, $s\ge0$ for 4 cities with cost matrix $C$). All of it is one big LP.

#### 3.6 Feasible-descent direction methods

$x_{k+1}=x_k+\alpha_kd_k$ with $d_k$ **feasible and descent**; choose $\alpha_k$ by line search while keeping feasibility (the simplex method is an example: a negative reduced cost gives the direction, and the step is the largest feasible one).

**Feasible direction space** $F_x=\{d\neq0:\ x+\lambda d\in F\ \text{for all small }\lambda>0\}$:
- Unconstrained: $F_x=\mathbb{R}^n$.
- $\{Ax=b\}$: $F_x=\{d:Ad=0\}$.
- $\{Ax\ge b\}$: $F_x=\{d:a_id\ge0\ \forall i\in B(x)\}$, $B(x)=\{i:a_ix=b_i\}$ (binding set).
- $\{Ax=b,\ x\ge0\}$: $F_x=\{d:Ad=0,\ d_i\ge0\ \forall i\in B(x)\}$, $B(x)=\{i:x_i=0\}$.

**Finding $d_k$ (SLP):** $\min\nabla f(x_k)^Td$ s.t. $d\in F_{x_k}$, $\|d\|\le1$. With the 1- or $\infty$-norm this is an **LP** (SLP / Frank–Wolfe style); with the 2-norm a QP (closed form if only equalities). If the optimal value is $<0$, we found a feasible descent direction.

#### 3.7 Reduced gradient method

For $Ax=b$, split $x=(x_B,x_N)$ with $A_B$ invertible: $x_B=A_B^{-1}(b-A_Nx_N)$. Substituting gives an **unconstrained** problem in $x_N$:
$$\min_{x_N}f\big(A_B^{-1}(b-A_Nx_N),\ x_N\big).$$
Its gradient w.r.t. $x_N$ is the **reduced gradient**, the same as the reduced cost in the simplex method.

#### 3.8 Newton's method on the KKT equations

**Linear equality constraints.** KKT for $\min f(x)$ s.t. $Ax=b$:
$$g(x,y)=\begin{pmatrix}\nabla f(x)-A^Ty\\Ax-b\end{pmatrix}=0,\qquad\nabla g(x,y)=\begin{pmatrix}\nabla^2f(x)&-A^T\\A&0\end{pmatrix}.$$
Newton: $\binom{x_{k+1}}{y_{k+1}}=\binom{x_k}{y_k}-\nabla g(x_k,y_k)^{-1}g(x_k,y_k)$.
For a **QP** $\min\frac12x^TQx+c^Tx$ s.t. $Ax=b$, the KKT system is linear and one solve gives the answer:
$$\begin{pmatrix}Q&-A^T\\A&0\end{pmatrix}\begin{pmatrix}x\\y\end{pmatrix}=\begin{pmatrix}-c\\b\end{pmatrix}.$$

**Portfolio QP example:** $\min x_1^2+2x_2^2-2x_1x_2-0.5x_1-0.5x_2$ s.t. $x_1+x_2=1$.
Here $Q=\begin{pmatrix}2&-2\\-2&4\end{pmatrix}$, $c=(-0.5,-0.5)$, giving
$$\begin{pmatrix}2&-2&-1\\-2&4&-1\\1&1&0\end{pmatrix}\begin{pmatrix}x_1\\x_2\\y\end{pmatrix}=\begin{pmatrix}0.5\\0.5\\1\end{pmatrix}.$$
🔍 *Solving:* subtract row 2 from row 1: $4x_1-6x_2=0\Rightarrow x_1=1.5x_2$. Row 3: $2.5x_2=1\Rightarrow x_2=0.4$, $x_1=0.6$. Row 1: $y=1.2-0.8-0.5=-0.1$. So $x^*=(0.6,0.4)$.

**Nonlinear equality constraints** $h(x)=0$:
$$g(x,y)=\begin{pmatrix}\nabla f(x)-\sum_iy_i\nabla h_i(x)\\h(x)\end{pmatrix},\qquad\nabla g=\begin{pmatrix}\nabla^2f(x)-\sum_iy_i\nabla^2h_i(x)&-\nabla h(x)^T\\\nabla h(x)&0\end{pmatrix}.$$
**Spherical example:** $\min(x_1-1)^2+(x_2-1)^2$ s.t. $x_1^2+x_2^2-1=0$:
$$g=\begin{pmatrix}2(x_1-1)-2yx_1\\2(x_2-1)-2yx_2\\x_1^2+x_2^2-1\end{pmatrix},\qquad\nabla g=\begin{pmatrix}2-2y&0&-2x_1\\0&2-2y&-2x_2\\2x_1&2x_2&0\end{pmatrix}.$$
Newton converges from a good initial $(x,y)$ to $x^*=(\frac1{\sqrt2},\frac1{\sqrt2})$ (the point of the circle closest to $(1,1)$). 🔍 Check: from rows 1–2, $x_1=x_2=\frac{1}{1-y}$; on the circle $x_i=\frac1{\sqrt2}$, so $y=1-\sqrt2$.

---

### 4. Summary

✅ Approximate KKT residuals are the practical stopping rule.
✅ Projection onto $x\ge0$: clip. Onto $Ax=b$: multiply by $P_A=I-A^T(AA^T)^{-1}A$.
✅ Feasible directions: $Ad=0$ plus $d_i\ge0$ on binding bounds; SLP finds the best one by an LP.
✅ Reduced gradient = eliminate basic variables.
✅ Newton on KKT: solve $\nabla g\,\Delta=-g$; for an equality-constrained QP, one linear solve.

<!-- LECTURE END: yinyuye-lecture12 -->

---

<!-- LECTURE START: frankwolfe -->
<a id="frankwolfe"></a>

## Lecture 42: The Frank–Wolfe Algorithm (reading: S. Pokutta, "The Frank–Wolfe algorithm: a short introduction")

*Source file:* `FrankWolfe.pdf`

### 1. What this lecture is about

This is a research-article reading (arXiv 2311.05313). Its core ideas, explained from scratch:
1. Three ways to handle constraints: projection, barrier, **linear minimization oracle (LMO)**
2. The Frank–Wolfe (conditional gradient) algorithm, which is **projection-free**
3. Smoothness ($L$) and strong convexity ($\mu$): quadratic upper/lower bounds
4. **Frank–Wolfe gap**: a computable optimality certificate
5. The **$O(1/t)$ convergence proof**, every step
6. A matching lower bound, sparsity, and applications (approximate Carathéodory)

---

### 2. Background you need first

- Convexity: $f(y)-f(x)\ge\nabla f(x)^T(y-x)$ (Convexity Part 2).
- Projected gradient descent (Nov 5–6 blocks); Frank–Wolfe step from the Nov 6 block.
- Extreme points / vertices of a polytope (BFS block). $\langle a,b\rangle=a^Tb$.
- Diameter of a set $P$: $D=\max_{x,y\in P}\|x-y\|$.

---

### 3. Main content with full derivations

#### 3.1 The problem and three access models

$$\min_{x\in P}f(x),\qquad P\subseteq\mathbb{R}^n\ \text{compact convex (e.g. a polytope)},\quad f\ \text{smooth convex}.$$
How can an algorithm "see" $P$?
1. **Projection** $\Pi_P(x)=\arg\min_{y\in P}\|x-y\|$, used by projected GD: $x_{t+\frac12}=x_t-\gamma_t\nabla f(x_t)$, then $x_{t+1}=\Pi_P(x_{t+\frac12})$ (typically $\gamma_t=1/L$). For complicated $P$, each projection is itself a hard optimization problem.
2. **Barrier function**, e.g. $-\sum_i\log(b_i-A_ix)$ for $\{Ax\le b\}$, used by interior-point methods: $x_\mu=\arg\min_x\langle c,x\rangle+\mu\sum\log(\dots)$, $\mu\to0$. This needs an explicit description of $P$ and usually second-order information.
3. **Linear Minimization Oracle (LMO):** given $c$, return $v\in\arg\min_{x\in P}\langle c,x\rangle$. Optimising a *linear* function over $P$ is often far cheaper than projecting (e.g. over the nuclear-norm ball, matching polytopes). WLOG $v$ is an **extreme point** ("atom").

We also assume a **first-order oracle**: given $x$, it returns $f(x)$ and $\nabla f(x)$.

#### 3.2 The algorithm

```
Input: x_0 ∈ P (an atom), smooth convex f
for t = 0, 1, 2, ...:
    v_t  = argmin_{v ∈ P} ⟨∇f(x_t), v⟩      # LMO call
    γ_t  = 2/(t+2)                           # "open-loop" step size
    x_{t+1} = (1 − γ_t) x_t + γ_t v_t        # convex combination ⇒ stays in P
```
- **Why it works:** $v_t$ minimises the linear model $f(x_t)+\langle\nabla f(x_t),x-x_t\rangle$ over $P$, so $v_t-x_t$ is a descent direction aligned with $-\nabla f$ that never leaves $P$.
- **Feasibility for free:** $x_{t+1}$ is a convex combination of points of $P$.
- **Sparsity:** $x_t$ is a convex combination of at most $t+1$ atoms.
- **Affine invariance:** rescaling the problem doesn't change the iterates.

#### 3.3 Smoothness and strong convexity

- **Convex:** $f(y)-f(x)\ge\langle\nabla f(x),y-x\rangle$ (linear **lower** bound).
- **$\mu$-strongly convex:** $f(y)-f(x)\ge\langle\nabla f(x),y-x\rangle+\frac\mu2\|y-x\|^2$ (quadratic **lower** bound).
- **$L$-smooth:** $f(y)-f(x)\le\langle\nabla f(x),y-x\rangle+\frac L2\|y-x\|^2$ (quadratic **upper** bound).
- Gradient versions: $\langle\nabla f(y)-\nabla f(x),y-x\rangle\le L\|y-x\|^2$ (smooth) and $\ge\mu\|y-x\|^2$ (strongly convex).
- For twice differentiable $f$: $L$ bounds the largest eigenvalue of $\nabla^2f$, and $\mu$ the smallest.

🔍 *Example:* $f=\|x\|^2$ has $\nabla^2f=2I$, so $L=\mu=2$.

#### 3.4 Primal gap, dual gap and Frank–Wolfe gap

**Lemma 4.1.** For convex $f$ and all $x\in P$:
$$\underbrace{f(x)-f(x^*)}_{\text{primal gap}}\ \le\ \underbrace{\langle\nabla f(x),x-x^*\rangle}_{\text{dual gap}}\ \le\ \underbrace{\max_{v\in P}\langle\nabla f(x),x-v\rangle}_{\text{Frank–Wolfe gap}} .$$
🔍 *Proof:* the first inequality is convexity with $y=x^*$, rearranged: $f(x^*)\ge f(x)+\langle\nabla f(x),x^*-x\rangle$. The second holds because $x^*$ is one particular $v\in P$, and the max is at least as big.

The FW gap is computed for free in every iteration (it uses the $v_t$ from the LMO). It is an **upper bound on the error** and a natural **stopping criterion**.

**Lemma 4.2 (first-order optimality).** $x^*$ is optimal ⟺ $\langle\nabla f(x^*),x^*-v\rangle\le0$ for all $v\in P$ ⟺ FW gap at $x^*$ is 0. (⚠️ With constraints, $\nabla f(x^*)$ need **not** be 0.)

#### 3.5 Progress per step (from smoothness)

**Lemma 4.3.** With $x_{t+1}=(1-\gamma_t)x_t+\gamma_tv_t$:
$$f(x_t)-f(x_{t+1})\ \ge\ \gamma_t\langle\nabla f(x_t),x_t-v_t\rangle-\frac{L}{2}\gamma_t^2\|x_t-v_t\|^2 .$$
🔍 *Proof:* in the smoothness inequality take $x=x_t$, $y=x_{t+1}$, so $y-x=\gamma_t(v_t-x_t)$:
$$f(x_{t+1})-f(x_t)\le\gamma_t\langle\nabla f(x_t),v_t-x_t\rangle+\tfrac L2\gamma_t^2\|v_t-x_t\|^2 .$$
Multiply by $-1$ and use $\langle\nabla f,v_t-x_t\rangle=-\langle\nabla f,x_t-v_t\rangle$.

#### 3.6 Theorem 4.4: $f(x_t)-f(x^*)\le\dfrac{2LD^2}{t+2}$

So for accuracy $\varepsilon$, $t\ge\frac{2LD^2}{\varepsilon}$ iterations suffice (**sublinear**, $O(1/t)$).

**Proof.**
1. From Lemma 4.3, using that $v_t$ maximises the FW gap, so $\langle\nabla f(x_t),x_t-v_t\rangle\ge f(x_t)-f(x^*)$ (Lemma 4.1):
$$f(x_t)-f(x_{t+1})\ge\gamma_t\big(f(x_t)-f(x^*)\big)-\gamma_t^2\tfrac L2\|x_t-v_t\|^2 .$$
2. Let $h_t=f(x_t)-f(x^*)$ (primal gap). Subtract $f(x^*)$ on both sides, rearrange, and bound $\|x_t-v_t\|\le D$:
$$h_{t+1}\le(1-\gamma_t)h_t+\gamma_t^2\frac{LD^2}{2}. \tag{4.2}$$
3. **Induction, base case $t=0$:** $\gamma_0=1$, so $h_1\le0+\frac{LD^2}{2}\le\frac{2LD^2}{3}$ ✓ (the claim for $t=1$ is $\frac{2LD^2}{3}$).
4. **Step:** assume $h_t\le\frac{2LD^2}{t+2}$. With $\gamma_t=\frac2{t+2}$, so $1-\gamma_t=\frac{t}{t+2}$:
$$h_{t+1}\le\frac{t}{t+2}\cdot\frac{2LD^2}{t+2}+\frac{4}{(t+2)^2}\cdot\frac{LD^2}{2}=\frac{2LD^2(t+1)}{(t+2)^2}.$$
🔍 (Combine: $\frac{2LD^2t}{(t+2)^2}+\frac{2LD^2}{(t+2)^2}=\frac{2LD^2(t+1)}{(t+2)^2}$.)
5. We need $\frac{t+1}{(t+2)^2}\le\frac{1}{t+3}$ ⟺ $(t+1)(t+3)\le(t+2)^2$ ⟺ $t^2+4t+3\le t^2+4t+4$ ✓. So $h_{t+1}\le\frac{2LD^2}{t+3}$. ∎

**Theorem 4.5:** the running minimum of the FW gaps satisfies $\min_{0\le\tau\le t}\text{gap}_\tau\le\frac{6.75\,LD^2}{t+2}$. So the computable gap also shrinks like $O(1/t)$.

#### 3.7 Lower bound: $O(1/t)$ cannot be improved in general

$\min\|x\|^2$ over the probability simplex $P=\operatorname{conv}\{e_1,\dots,e_n\}$ ($L=2$, $D=\sqrt2$). Optimum $x^*=(\frac1n,\dots,\frac1n)$, $f^*=\frac1n$.
🔍 After $t<n$ LMO calls, an iterate uses at most $t$ (or $t+1$) vertices, so it has at most that many non-zeros. The best point on $k$ coordinates is uniform on them, with value $\frac1k$. Hence $f(x_t)-f^*\ge\frac1t-\frac1n$: **$\Omega(LD^2/\varepsilon)$ LMO calls are necessary** for any LMO-based method.

#### 3.8 Other topics in the article (brief)

- **Nonconvex $f$:** FW still drives the FW gap to 0 (at rate $O(1/\sqrt t)$), finding stationary points.
- **Dual prices:** the LMO answers give price information (like LP duals).
- **Adaptive step sizes:** estimate $L$ on the fly and use the "short step" $\gamma_t=\min\big\{1,\frac{\langle\nabla f(x_t),x_t-v_t\rangle}{L\|x_t-v_t\|^2}\big\}$ (minimiser of the smoothness upper bound), often much faster than $\frac{2}{t+2}$.
- **Application 1: approximate Carathéodory.** Any $x^*\in P$ can be approximated within $\varepsilon$ by a convex combination of only $O(D^2/\varepsilon^2)$ vertices. 🔍 *Proof idea:* run FW on $f(x)=\|x-x^*\|^2$ (so $f^*=0$, $L=2$) until $f\le\varepsilon^2$, which takes $O(D^2/\varepsilon^2)$ steps; each step adds one vertex. (Compare exact Carathéodory, which needs $n+1$ points: Sep 11 block.)
- **Application 2: separating hyperplanes.** If $x^*\notin P$, the FW gap on $\|x-x^*\|^2$ yields a separating hyperplane.
- Software: FrankWolfe.jl (Julia).

---

### 4. Worked example

$\min(x_1-0.2)^2+(x_2-0.3)^2$ over the simplex $\{x\ge0,\ x_1+x_2\le1\}$ (vertices $(0,0),(1,0),(0,1)$), $x_0=(1,0)$.
- $t=0$: $\nabla f=(1.6,-0.6)$. LMO: $\langle\nabla f,v\rangle$ at the vertices is $0,\ 1.6,\ -0.6$, so $v_0=(0,1)$. $\gamma_0=1$ ⇒ $x_1=(0,1)$.
- $t=1$: $\nabla f=(-0.4,1.4)$; values $0,-0.4,1.4$ ⇒ $v_1=(1,0)$; $\gamma_1=\frac23$ ⇒ $x_2=\frac13(0,1)+\frac23(1,0)=(0.667,0.333)$.
- $t=2$: $\nabla f=(0.933,0.067)$; values $0,0.933,0.067$ ⇒ $v_2=(0,0)$; $\gamma_2=\frac12$ ⇒ $x_3=(0.333,0.167)$. …

The iterates approach $x^*=(0.2,0.3)$, which is inside $P$. Convergence is slow there (the zig-zag that the lower bound describes).

---

### 5. Summary

✅ FW: $v_t=\text{LMO}(\nabla f(x_t))$, $x_{t+1}=(1-\gamma_t)x_t+\gamma_tv_t$, $\gamma_t=\frac2{t+2}$. Projection-free and always feasible.
✅ FW gap $\max_v\langle\nabla f(x),x-v\rangle\ge f(x)-f^*$: a free certificate.
✅ Rate $f(x_t)-f^*\le\frac{2LD^2}{t+2}$ (proof: smoothness + gap + induction); optimal for LMO methods.
✅ Iterates are sparse combinations of vertices ⇒ approximate Carathéodory.

<!-- LECTURE END: frankwolfe -->

---

<!-- LECTURE START: duality1 -->
<a id="duality1"></a>

## Lecture 43: Lagrangian Duality (Part 1): dual function, dual problem, examples, LP & QP duals, a duality gap

*Source file:* `Duality1.pdf`

### 1. What this lecture is about

LP duality generalises to **any** constrained problem through the Lagrangian. This lecture:
1. Defines the **dual function** $\theta(\lambda,\mu)=\inf_{x\in X}L(x,\lambda,\mu)$ and proves it is a lower bound on $p^*$
2. Works three examples in full, including one where we keep simple constraints in $X$
3. Rederives the **LP dual** and derives the **QP dual**
4. Shows $\theta$ is **concave** (so the dual is always a convex problem), weak duality, and a **discrete example with a positive duality gap**

---

### 2. Background you need first

- Lagrangian $L=f+\sum\lambda_ig_i+\sum\mu_jh_j$ (KKT blocks); LP duality.
- "inf" = greatest lower bound (Topology block). Minimising a 1-D parabola $x^2-\lambda x$: unconstrained minimiser $x=\lambda/2$; on an interval, clip $\lambda/2$ into the interval (the parabola is convex).
- Pointwise minimum of affine (linear + constant) functions is concave.

---

### 3. Main content with full derivations

#### 3.1 Dual function and dual problem

**Primal:**
$$p^*=\min f(x)\ \text{ s.t. }\ g_i(x)\le0\ (i\in I)\ [\lambda_i\ge0],\quad h_j(x)=0\ (j\in E)\ [\mu_j],\quad x\in X .$$
Constraints are split: the "hard" ones get multipliers; the simple set $X$ (e.g. bounds, integrality) is kept as is. $\Omega$ = all of them together; $f(x^*)=p^*$.

**Lagrangian:** $L(x,\lambda,\mu)=f(x)+\sum_{i\in I}\lambda_ig_i(x)+\sum_{j\in E}\mu_jh_j(x)$.

**Key inequality.** For every feasible $x\in\Omega$ and every $\lambda\ge0$, $\mu$:
$$L(x,\lambda,\mu)=f(x)+\underbrace{\sum\lambda_ig_i(x)}_{\le0\ (\lambda\ge0,\ g\le0)}+\underbrace{\sum\mu_jh_j(x)}_{=0}\ \le\ f(x).$$
Taking the infimum over the **larger** set $X\supseteq\Omega$ only makes the left side smaller:
$$\theta(\lambda,\mu):=\inf_{x\in X}L(x,\lambda,\mu)\ \le\ L(x,\lambda,\mu)\ \le\ f(x)\quad\forall x\in\Omega\ \Rightarrow\ \theta(\lambda,\mu)\le f(x^*)=p^* .$$
$\theta$ is the **dual function**. Every $\lambda\ge0$ gives a lower bound. The best bound is the **dual problem**:
$$d^*=\max\ \theta(\lambda,\mu)\ \text{ s.t. }\ \lambda\ge0\ (\lambda\in\mathbb{R}^{|I|}),\ \mu\in\mathbb{R}^{|E|}.$$

#### 3.2 Example 1

$$\min x_1^2+x_2^2\ \text{ s.t. }\ x_1+x_2-4\ge0,\ \ x_1\ge0,\ x_2\ge0 .$$
(Picture: circles around the origin; the line $x_1+x_2=4$; the closest point is $(2,2)$, so $p^*=8$.)
Dualise only $-x_1-x_2+4\le0\ [\lambda]$; keep $X=\{x\in\mathbb{R}^2:x_1\ge0,x_2\ge0\}$.
$$L(x,\lambda)=x_1^2+x_2^2+\lambda(-x_1-x_2+4)=(x_1^2-\lambda x_1)+(x_2^2-\lambda x_2)+4\lambda .$$
$$\theta(\lambda)=\inf_{x_1\ge0}\{x_1^2-\lambda x_1\}+\inf_{x_2\ge0}\{x_2^2-\lambda x_2\}+4\lambda .$$
(It separates into two 1-D problems.) Ignoring $x_1\ge0$: $2x_1-\lambda=0\Rightarrow x_1^*=\frac\lambda2$, which is $\ge0$ anyway since we only care about $\lambda\ge0$ (and it is the minimiser, since $x_1^2-\lambda x_1$ is convex). Value: $\frac{\lambda^2}{4}-\frac{\lambda^2}{2}=-\frac{\lambda^2}{4}$. Same for $x_2$.
$$\theta(\lambda)=-\frac{\lambda^2}{4}-\frac{\lambda^2}{4}+4\lambda=4\lambda-\frac{\lambda^2}{2}.$$
$d^*=\max_{\lambda\ge0}4\lambda-\frac{\lambda^2}2$: derivative $4-\lambda=0\Rightarrow\lambda^*=4$, $d^*=16-8=8$. So **$p^*\ge8$** from the dual alone, and indeed $p^*=8$ (no gap).

#### 3.3 Example 2 (bounds kept in $X$)

$$\min x_1^2+x_2^2\ \text{ s.t. }\ -x_1-x_2+4\le0\ [\lambda],\qquad X=\{0\le x_1\le1,\ 1\le x_2\le5\}.$$
(Feasible region: the small triangle near $(1,3)$ to $(1,5)$; optimum $(1,3)$, $p^*=10$.)
$$\theta(\lambda)=\underbrace{\inf_{0\le x_1\le1}(x_1^2-\lambda x_1)}_{\theta_1(\lambda)}+\underbrace{\inf_{1\le x_2\le5}(x_2^2-\lambda x_2)}_{\theta_2(\lambda)}+4\lambda .$$
**$\theta_1$:** unconstrained minimiser $\frac\lambda2$. If $0\le\frac\lambda2\le1$ ($\lambda\le2$): $x_1^*=\frac\lambda2$, value $-\frac{\lambda^2}4$. Else ($\lambda\ge2$) the parabola is still decreasing at 1, so $x_1^*=1$: value $1-\lambda$.
$$\theta_1(\lambda)=\begin{cases}-\frac{\lambda^2}{4}&0\le\lambda\le2\\1-\lambda&\lambda\ge2\end{cases}$$
**$\theta_2$:** $x_2^*=\frac\lambda2$ if $1\le\frac\lambda2\le5$; $x_2^*=1$ if $\frac\lambda2\le1$; $x_2^*=5$ if $\frac\lambda2\ge5$:
$$\theta_2(\lambda)=\begin{cases}1-\lambda&\lambda\le2\\-\frac{\lambda^2}{4}&2\le\lambda\le10\\25-5\lambda&\lambda\ge10\end{cases}$$
**Add** (for $\lambda\le2$ and for $2\le\lambda\le10$ the sum $\theta_1+\theta_2$ is the same expression, $1-\lambda-\frac{\lambda^2}4$):
$$\theta(\lambda)=\begin{cases}4\lambda+1-\lambda-\frac{\lambda^2}{4}=3\lambda-\frac{\lambda^2}{4}+1&\lambda\le10\\4\lambda+(1-\lambda)+(25-5\lambda)=26-2\lambda&\lambda\ge10\end{cases}$$
Maximise: on $[0,10]$, $3-\frac\lambda2=0\Rightarrow\lambda=6$, $\theta(6)=18-9+1=10$ ✓. For $\lambda\ge10$, $26-2\lambda\le6$ (e.g. $\theta(10)=6$). So **$d^*=10$** and $p^*\ge10$; with $p^*=10$ at $(1,3)$, no gap.

#### 3.4 Example 3: dual of an LP (recovering LP duality)

$$\min c^Tx\ \text{ s.t. }\ -a_i^Tx+b_i\le0\ [\lambda_i]\ (i=1..m),\quad-x_i\le0\ [\tau_i]\ (i=1..n),\quad X=\mathbb{R}^n .$$
$$L(x,\lambda,\tau)=c^Tx+\sum_i\lambda_i(b_i-a_i^Tx)-\sum_i\tau_ix_i=\Big(c-\sum_i\lambda_ia_i-\tau\Big)^Tx+\sum_i\lambda_ib_i .$$
This is **linear in $x$**. A linear function $v^Tx$ over all of $\mathbb{R}^n$ is unbounded below unless $v=0$:
$$\theta(\lambda,\tau)=\inf_{x\in\mathbb{R}^n}L=\begin{cases}\sum_i\lambda_ib_i&\text{if }c-\sum_i\lambda_ia_i-\tau=0\\-\infty&\text{otherwise}\end{cases}$$
Dual: $\max\sum_i\lambda_ib_i$ s.t. $c-\sum\lambda_ia_i-\tau=0$, $\lambda\ge0$, $\tau\ge0$. Eliminate $\tau=c-\sum\lambda_ia_i\ge0$:
$$\max\ \sum_{i=1}^m\lambda_ib_i\ \text{ s.t. }\ c-\sum_i\lambda_ia_i\ge0\ (\text{i.e. }A^T\lambda\le c),\ \ \lambda\ge0 .$$
Exactly the LP dual ✓.

#### 3.5 Example 4: dual of a QP

$$\min\tfrac12x^TQx+c^Tx\ \text{ s.t. }\ Ax\le b\ [\lambda],\quad x\in\mathbb{R}^n,\qquad Q\succ0\ \text{symmetric},\ A\in\mathbb{R}^{m\times n}.$$
$$L(x,\lambda)=\tfrac12x^TQx+c^Tx+\lambda^T(Ax-b).$$
$L$ is convex in $x$, so its infimum is where $\nabla_xL=Qx+c+A^T\lambda=0$, i.e. $x=Q^{-1}(-c-A^T\lambda)$. Substitute, writing $u=c+A^T\lambda$ so $x=-Q^{-1}u$:
- $\frac12x^TQx=\frac12u^TQ^{-1}QQ^{-1}u=\frac12u^TQ^{-1}u$ (the middle $QQ^{-1}=I$);
- $c^Tx+\lambda^TAx=(c+A^T\lambda)^Tx=u^Tx=-u^TQ^{-1}u$;
- $-\lambda^Tb$ stays.
$$\theta(\lambda)=\tfrac12u^TQ^{-1}u-u^TQ^{-1}u-\lambda^Tb=-\tfrac12(c+A^T\lambda)^TQ^{-1}(c+A^T\lambda)-\lambda^Tb .$$
Dual: $\max_{\lambda\ge0}\theta(\lambda)$, again a **QP** (concave quadratic). With $Q=I$: $\theta(\lambda)=-\frac12\|c+A^T\lambda\|^2-\lambda^Tb$. The dual involves the data only through inner products $a_i^Ta_j$, the **kernel idea in SVMs**.

**HW:** $\min e^x$ s.t. $x^2-1\le0\ [\lambda]$.
🔍 *Solution.* $p^*=e^{-1}$ at $x=-1$. $\theta(\lambda)=\inf_x e^x+\lambda(x^2-1)$. For $\lambda=0$: $\inf e^x=0$. For $\lambda>0$ the minimiser solves $e^x+2\lambda x=0$ (so $x<0$). At $\lambda=\frac{e^{-1}}{2}$ the minimiser is $x=-1$ (check: $e^{-1}-e^{-1}=0$ ✓), giving $\theta=e^{-1}+\lambda(1-1)=e^{-1}$. So $d^*=e^{-1}=p^*$ (convex problem with Slater point $x=0$).

#### 3.6 General facts

$$\theta(\lambda,\mu)=\inf_{x\in X}\Big\{f(x)+\sum_i\lambda_ig_i(x)+\sum_j\mu_jh_j(x)\Big\},\qquad d^*=\max_{\lambda\ge0}\theta(\lambda,\mu).$$
- **The dual is always a convex optimization problem:** $\theta$ is **concave** in $(\lambda,\mu)$ (HW: prove).
  🔍 *Proof:* for each fixed $x$, $L(x,\lambda,\mu)$ is **affine** in $(\lambda,\mu)$. $\theta$ is the pointwise infimum of a family of affine functions, so for $\alpha\in[0,1]$:
  $\theta(\alpha z_1+(1-\alpha)z_2)=\inf_x[\alpha L(x,z_1)+(1-\alpha)L(x,z_2)]\ge\alpha\inf_xL(x,z_1)+(1-\alpha)\inf_xL(x,z_2)$. That is concavity, even if the primal is non-convex.
- **Weak duality:** $d^*\le p^*$ (shown in §3.1).
- **Strong duality** ($d^*=p^*$)? Not always…

#### 3.7 A duality gap: discrete example

$$\min-2x_1+x_2\ \text{ s.t. }\ x_1+x_2=3\ [\mu],\qquad x\in X=\{(0,0),(0,4),(4,4),(4,0),(1,2),(2,1)\}.$$
Feasible points: $(1,2)$ with value $0$, and $(2,1)$ with value $-3$. So $p^*=-3$.
$L(x,\mu)=-2x_1+x_2+\mu(x_1+x_2-3)$. Evaluate at each point of $X$:

| $x$ | $L$ |
|---|---|
| $(0,0)$ | $-3\mu$ |
| $(0,4)$ | $4+\mu$ |
| $(4,4)$ | $-4+5\mu$ |
| $(4,0)$ | $-8+\mu$ |
| $(1,2)$ | $0$ |
| $(2,1)$ | $-3$ |

$\theta(\mu)=\min$ of these lines (the lower envelope in the plot). The maximum of the envelope is where $-3\mu=-8+\mu\Rightarrow\mu^*=2$, $\theta=-6$. (Check at $\mu=2$: other values are $6,6,0,-3$, all larger than $-6$ ✓. For $\mu<2$, $-8+\mu<-6$; for $\mu>2$, $-3\mu<-6$.)
$$d^*=-6\ <\ p^*=-3 .$$
**Duality gap $=3>0$.** Because $X$ is non-convex (discrete), strong duality fails.

---

### 4. Summary

✅ $\theta(\lambda,\mu)=\inf_{x\in X}L$ is a lower bound on $p^*$ for every $\lambda\ge0$; the dual maximises it.
✅ Keep easy constraints in $X$; the inner problem often separates by variable.
✅ LP dual and QP dual ($-\frac12(c+A^T\lambda)^TQ^{-1}(c+A^T\lambda)-\lambda^Tb$) follow mechanically.
✅ $\theta$ is always concave ⇒ the dual is a convex problem. Weak duality always; strong duality can fail (discrete example: gap 3).

<!-- LECTURE END: duality1 -->

---

<!-- LECTURE START: duality2 -->
<a id="duality2"></a>

## Lecture 44: Lagrangian Duality (Part 2): strong duality for convex problems, KKT from strong duality, projection onto a half-space

*Source file:* `Duality2.pdf`

### 1. What this lecture is about

1. **Strong duality** for convex problems under **Slater's condition**: $p^*=d^*$, and the optimal $(x^*,\lambda^*,\mu^*)$ is a KKT point
2. **Proof that strong duality ⇒ complementary slackness and KKT** (the "sandwich" argument)
3. Worked example: minimum-norm point in a half-space, solved via the dual

---

### 2. Background you need first

- Dual function $\theta=\inf_{x\in X}L$, weak duality $p^*\ge d^*$ (Duality Part 1).
- KKT conditions; Slater's condition (KKT for Convex Problems block).

---

### 3. Main content with full derivations

#### 3.1 Convex problem and its dual

$$p^*=\min f(x)\ \text{ s.t. }\ g_i(x)\le0\ (i\in I),\ \ h_j(x)=0\ (j\in E),\ \ x\in X,$$
with $f,g_i$ **convex** and $h_j$ **linear/affine**.
$$\theta(\lambda,\mu)=\inf_{x\in X}\Big\{f(x)+\sum_{i\in I}\lambda_ig_i(x)+\sum_{j\in E}\mu_jh_j(x)\Big\},\qquad d^*=\max_{\lambda\ge0,\ \mu}\theta(\lambda,\mu),\qquad p^*\ge d^*\ (\text{weak}).$$

#### 3.2 When is $p^*=d^*$ (strong duality)?

> **Theorem (Slater).** If the problem is convex and $\exists x^0$ (Slater point) with $g_i(x^0)<0\ \forall i\in I$ and $h_j(x^0)=0\ \forall j\in E$, then $p^*=d^*$. Moreover, if $\lambda^*,\mu^*$ are dual optimal and $x^*$ is primal optimal, then $(x^*,\lambda^*,\mu^*)$ is a KKT point.

#### 3.3 Strong duality ⇒ complementary slackness ⇒ KKT

Assume strong duality, $x^*$ primal optimal and $(\lambda^*,\mu^*)$ dual optimal. Then
$$f(x^*)=\theta(\lambda^*,\mu^*)\qquad(p^*=d^*)$$
$$=\inf_{x\in X}\Big\{f(x)+\sum_i\lambda_i^*g_i(x)+\sum_j\mu_j^*h_j(x)\Big\}\qquad(\text{definition of }\theta)$$
$$\le f(x^*)+\sum_i\lambda_i^*\underbrace{g_i(x^*)}_{\le0}+\sum_j\mu_j^*\underbrace{h_j(x^*)}_{=0}\qquad(\text{inf}\le\text{value at the particular point }x^*)$$
$$\le f(x^*)\qquad(\lambda_i^*\ge0,\ g_i(x^*)\le0\ \Rightarrow\ \text{the sum is}\le0).$$
We started and ended with $f(x^*)$, so **every $\le$ is an equality**:
1. $\sum_i\lambda_i^*g_i(x^*)=0$. Each term is $\le0$, so each term is 0: **$\lambda_i^*g_i(x^*)=0\ \forall i$** (complementary slackness).
2. The infimum of $L(\cdot,\lambda^*,\mu^*)$ over $X$ is **attained at $x^*$**. If $X=\mathbb{R}^n$ (or $x^*$ is interior to $X$) and things are differentiable, the gradient vanishes there:
$$\nabla f(x^*)+\sum_i\lambda_i^*\nabla g_i(x^*)+\sum_j\mu_j^*\nabla h_j(x^*)=0 .$$

Together with primal feasibility ($g_i(x^*)\le0$, $h_j(x^*)=0$) and $\lambda_i^*\ge0$: **these are exactly the first-order KKT conditions.**

#### 3.4 Example: closest point to the origin in a half-space

$$\min\tfrac12\|x\|^2\ \text{ s.t. }\ a^Tx\le b\ [\lambda],\qquad x\in\mathbb{R}^n=X,\quad a\in\mathbb{R}^n,\ b<0 .$$
($b<0$ means the origin is **not** feasible, so the answer is on the boundary line $a^Tx=b$; the dotted perpendicular from the origin in the picture.)
$$L(x,\lambda)=\tfrac12x^Tx+\lambda(a^Tx-b).$$
**Dual function:** $L$ is convex in $x$; set $\nabla_xL=x+\lambda a=0\Rightarrow x=-\lambda a$. Substitute:
$$\theta(\lambda)=\tfrac12\lambda^2a^Ta+\lambda\big(a^T(-\lambda a)-b\big)=\tfrac12\lambda^2a^Ta-\lambda^2a^Ta-\lambda b=-\tfrac12\lambda^2a^Ta-\lambda b .$$
**Dual problem:** $\max_{\lambda\ge0}\ -\frac{\lambda^2}{2}a^Ta-\lambda b$.
Ignore $\lambda\ge0$ first: derivative $-\lambda a^Ta-b=0$ ⇒
$$\lambda^*=\frac{-b}{a^Ta}\ \ge0\quad(\text{since }b<0)\ ✓.$$
**Primal solution** (from the KKT stationarity $x=-\lambda a$):
$$\boxed{x^*=-\lambda^*a=\Big(\frac{b}{a^Ta}\Big)a}$$
🔍 *Checks.* Feasibility: $a^Tx^*=\frac{b}{a^Ta}a^Ta=b$ ✓ (on the boundary; CS holds since $\lambda^*>0$). Values: $p^*=\frac12\|x^*\|^2=\frac12\frac{b^2}{(a^Ta)^2}a^Ta=\frac{b^2}{2a^Ta}$; $d^*=-\frac12\frac{b^2}{(a^Ta)^2}a^Ta+\frac{b^2}{a^Ta}=\frac{b^2}{2a^Ta}$. **$p^*=d^*$** ✓. Distance from the origin to the hyperplane: $\|x^*\|=\frac{|b|}{\|a\|}$, the familiar formula.

---

### 4. Summary

✅ Convex problem + Slater point ⇒ strong duality $p^*=d^*$.
✅ Strong duality: $f(x^*)=\theta(\lambda^*)\le L(x^*,\lambda^*)\le f(x^*)$ ⇒ CS, and $x^*$ minimises $L(\cdot,\lambda^*)$ ⇒ KKT.
✅ Recipe: solve $\nabla_xL=0$ for $x(\lambda)$, plug in to get $\theta(\lambda)$, maximise over $\lambda\ge0$, recover $x^*=x(\lambda^*)$.
✅ Projection of 0 onto $\{a^Tx\le b\}$, $b<0$: $x^*=\frac{b}{a^Ta}a$.

<!-- LECTURE END: duality2 -->

---

<!-- LECTURE START: dualascent -->
<a id="dualascent"></a>

## Lecture 45: Algorithms Using the Dual: dual ascent and dual decomposition

*Source file:* `DualAscent.pdf`

> 📌 `DualAscent.pdf` is Part I of `DualMethods.pdf`. Parts II and III (method of multipliers and ADMM) are in the block **"Dual Methods (full notes)"**.

### 1. What this lecture is about

Instead of solving the primal directly, **solve the dual by gradient ascent** and recover $x$ along the way. This lecture:
1. Sets up equality-constrained problems $\min f(x)$ s.t. $Ax=b$
2. **Dual ascent:** alternate "minimise the Lagrangian in $x$" and "gradient step in $\lambda$". Why the dual gradient is the constraint residual $Ax-b$.
3. **Dual decomposition:** if $f$ is separable into blocks, the $x$-step splits into independent small problems (parallelisable). This is the stepping stone to **ADMM**.

---

### 2. Background you need first

- Dual function $\theta(\lambda)=\inf_xL(x,\lambda)$, concave; strong duality (Duality blocks).
- Gradient **ascent** for maximising: $\lambda\leftarrow\lambda+\alpha\nabla\theta(\lambda)$.
- Block vectors: $x=(x^{(1)},\dots,x^{(B)})$ and $A=[A_1\ A_2\ \cdots\ A_B]$, so $Ax=\sum_iA_ix^{(i)}$.

---

### 3. Main content with full derivations

#### 3.1 Setting

$$p^*=\min f(x)\ \text{ s.t. }\ Ax=b\ [\lambda\in\mathbb{R}^m],\quad x\in\mathbb{R}^n=X,\qquad A\in\mathbb{R}^{m\times n},\ b\in\mathbb{R}^m .$$
$$\theta(\lambda)=\inf_{x\in\mathbb{R}^n}\big\{f(x)+\lambda^T(Ax-b)\big\},\qquad\text{dual: }\max_{\lambda\in\mathbb{R}^m}\theta(\lambda)$$
($\lambda$ is free because the constraint is an equality.)

#### 3.2 Dual ascent

```
Init λ^0, k = 0
Till convergence {
    x^{k+1} = argmin_{x∈R^n}  f(x) + (λ^k)ᵀ(Ax − b)      # minimise the Lagrangian
    λ^{k+1} = λ^k + α ∇θ(λ^k),   where  ∇θ(λ^k) = A x^{k+1} − b
    k = k+1
}
stop when ‖∇θ(λ)‖ = ‖Ax − b‖ < ε  (e.g. 1e−9)
```
🔍 **Why $\nabla\theta(\lambda^k)=Ax^{k+1}-b$.** Let $x(\lambda)$ be the minimiser of $L(\cdot,\lambda)$, so $\theta(\lambda)=f(x(\lambda))+\lambda^T(Ax(\lambda)-b)$. Differentiate by the chain rule:
$$\nabla\theta=\Big[\underbrace{\nabla f(x)+A^T\lambda}_{=\nabla_xL=0\text{ at the minimiser}}\Big]^T\frac{\partial x}{\partial\lambda}+(Ax(\lambda)-b)=Ax(\lambda)-b .$$
(The $x$-dependence drops out because $x(\lambda)$ already minimises $L$: Danskin's/envelope theorem.)

**Interpretation:** $\lambda$ is a *price* on violating $Ax=b$. If $(Ax-b)_i>0$, raise $\lambda_i$ so that the next $x$-minimisation is pushed to reduce $(Ax)_i$. At convergence $Ax=b$ (primal feasible) and $x$ minimises $L$ (stationarity): a KKT point.

🔍 *Small example:* $\min\frac12\|x\|^2$ s.t. $a^Tx=b$. $x$-step: $x+\lambda a=0\Rightarrow x=-\lambda a$. $\lambda$-step: $\lambda\leftarrow\lambda+\alpha(-\lambda a^Ta-b)$. Fixed point: $\lambda^*=-\frac{b}{a^Ta}$, $x^*=\frac{b}{a^Ta}a$ (same as Duality Part 2). It converges if $0<\alpha<\frac{2}{a^Ta}$, since the error is multiplied by $1-\alpha a^Ta$ each step.

⚠️ Dual ascent needs $f$ to behave well (e.g. strictly convex) so that $\arg\min_xL$ is unique and finite. For a linear $f$, $\inf_xL$ is $-\infty$ for most $\lambda$. This weakness motivates augmented Lagrangian methods and **ADMM** (noted in the margin).

#### 3.3 Dual decomposition

Suppose $f$ is **separable across $B$ blocks**:
$$f(x)=\sum_{i=1}^Bf_i(x^{(i)}),\qquad x=\begin{pmatrix}x^{(1)}\\\vdots\\x^{(B)}\end{pmatrix}.$$
*Example* ($B=3$): $\underbrace{x_1^2+x_2^2+x_1x_2}_{f_1(x_1,x_2)}+\underbrace{x_3^2}_{f_2(x_3)}+\underbrace{x_4x_5\cdots x_n}_{f_3}$. Each piece involves only its own variables.

Split $A=[A_1\ A_2\ \cdots\ A_B]$ conformally, so $Ax=\sum_{i=1}^BA_ix^{(i)}$. Then
$$L(x,\lambda)=\sum_{i=1}^Bf_i(x^{(i)})+\lambda^T\Big(\sum_{i=1}^BA_ix^{(i)}-b\Big)=\sum_{i=1}^B\Big[f_i(x^{(i)})+\lambda^TA_ix^{(i)}\Big]-\lambda^Tb .$$
(🔍 The notes write $\lambda^T(A_ix^{(i)}-b)$ inside the sum. Strictly, $\lambda^Tb$ appears once, not $B$ times, but it does not depend on $x$, so it doesn't affect any $\arg\min$.)
**The Lagrangian is a sum of terms, each depending on one block**, so minimising over $x$ splits into $B$ independent problems.

**Dual decomposition via dual ascent:**
```
Init λ^0, k = 0
Loop till convergence {
    for i = 1..B (in parallel):
        x^{(i),k+1} = argmin_{x^{(i)}} { f_i(x^{(i)}) + (λ^k)ᵀ A_i x^{(i)} }
    λ^{k+1} = λ^k + α ( Σ_i A_i x^{(i),k+1} − b )      # gather, update the shared price
    k = k+1
}
```
**Why useful:** a huge problem becomes $B$ small ones solved in parallel (distributed optimisation). They are coordinated only through the price vector $\lambda$, e.g. sharing a resource among $B$ agents.

---

### 4. Summary

✅ Dual ascent: $x^{k+1}=\arg\min_xL(x,\lambda^k)$; $\lambda^{k+1}=\lambda^k+\alpha(Ax^{k+1}-b)$.
✅ $\nabla\theta(\lambda)=Ax(\lambda)-b$ (envelope theorem); stop when the residual $\|Ax-b\|$ is tiny.
✅ Separable $f$ ⇒ the $x$-step decomposes into independent block problems (dual decomposition).
✅ Requires a well-behaved (e.g. strictly convex) $f$; ADMM fixes this.

<!-- LECTURE END: dualascent -->

---

<!-- LECTURE START: lec18-duality-thy -->
<a id="lec18-duality-thy"></a>

## Lecture 46: Duality Theory of Constrained Optimization (reading: R. M. Freund, MIT, 2004)

*Source file:* `lec18_duality_thy.pdf`

### 1. What this lecture is about

A reading that collects all of duality theory in one place. Several parts repeat earlier class notes; this block explains the **new** ideas and finishes the exercises the notes leave open:
1. Why duality matters in practice
2. The 3-step dual construction, with examples: LP, **binary IP**, **log-barrier**, QP; sign rules for $\le,\ge,=$ constraints
3. **Column geometry:** the dual = the highest supporting hyperplane of a resource–cost set
4. Concavity of the dual function; weak duality and corollaries
5. **Saddle-point optimality**
6. **Perturbation function**, the Strong Duality Theorem, stability, Slater ⇒ stability
7. Duality strategies: dualising "bad" constraints, decomposing large problems

---

### 2. Background you need first

- Lagrangian duality (Duality Parts 1–2), KKT, Slater; LP duality; sensitivity (dual = shadow price).
- **Subgradient:** $\gamma$ is a subgradient of a convex $h$ at $\bar y$ if $h(y)\ge h(\bar y)+\gamma^T(y-\bar y)$ for all $y$ (a "tangent plane below the graph" that works even without differentiability). $\partial h(\bar y)$ = the set of all subgradients.

---

### 3. Main content with full derivations

#### 3.1 Why duality matters

Duality appears in electrical networks (currents primal, voltages dual), markets (quantities primal, **prices** dual), and structures (tensions primal, displacements dual). Uses:
- bounding and recognising near-optimal solutions;
- **proving optimality** (a dual solution with equal value);
- **sensitivity**: a dual variable is the change in optimal value per unit increase of the RHS;
- **KKT multipliers** = dual optimal solutions;
- convergence analysis of algorithms; exploiting nicer structure in the dual.

**Convention:** a function may take the value $\pm\infty$. E.g. $\max_{x\in S}f(x)$ equals $\max_{x\in\mathbb{R}^n}h(x)$ with $h=f$ on $S$ and $-\infty$ outside.

#### 3.2 Constructing the dual (3 steps)

$$OP:\ z^*=\min f(x)\ \text{ s.t. }\ g_i(x)\le0\ (i=1..m),\ \ x\in X$$
($X$ can be $\mathbb{R}^n$, $\{x\ge0\}$, $\{Ax\le b\}$, integer vectors, …)
1. **Lagrangian:** $L(x,u)=f(x)+u^Tg(x)$.
2. **Dual function:** $L^*(u)=\min_{x\in X}L(x,u)$ (assumed easy to compute).
3. **Dual problem:** $D:\ v^*=\max_{u\ge0}L^*(u)$.

**Sign rules for mixed constraints:** $u_i\ge0$ for $g_i\le0$; $u_i\le0$ for $g_i\ge0$; $u_i$ free for $g_i=0$.

#### 3.3 Examples of dual construction

**(a) LP** $\min c^Tx$ s.t. $Ax\ge b$ ($X=\mathbb{R}^n$). Write $b-Ax\le0$: $L=c^Tx+u^T(b-Ax)=u^Tb+(c-A^Tu)^Tx$. Linear in free $x$ ⇒
$$L^*(u)=\begin{cases}u^Tb&A^Tu=c\\-\infty&A^Tu\neq c\end{cases}\qquad D:\ \max u^Tb\ \text{ s.t. }A^Tu=c,\ u\ge0 .$$

**(b) Binary integer program** $\min c^Tx$ s.t. $Ax\ge b$, $x_j\in\{0,1\}$ (the notes pose this as a question).
🔍 *Answer.* Keep $X=\{0,1\}^n$. $L=u^Tb+(c-A^Tu)^Tx=u^Tb+\sum_j(c-A^Tu)_jx_j$. Each $x_j$ is chosen independently: $x_j=1$ if its coefficient is negative, else 0. So
$$L^*(u)=u^Tb+\sum_{j=1}^n\min\{0,\ (c-A^Tu)_j\},\qquad D:\ \max_{u\ge0}L^*(u).$$
A concave, piecewise-linear dual, which gives **lower bounds** for the hard IP (the basis of Lagrangian relaxation in branch-and-bound). There may be a duality gap.

**(c) Log-barrier problem**
$$BP:\ \min5x_1+7x_2-4x_3-\sum_{j=1}^3\ln x_j\ \text{ s.t. }x_1+3x_2+12x_3=37,\ x>0 .$$
🔍 *Answer.* Free $u$ for the equality: $L=\sum_j(c_j-ua_j)x_j-\sum_j\ln x_j+37u$ with $c=(5,7,-4)$, $a=(1,3,12)$. For one variable, $\min_{x>0}kx-\ln x$: if $k>0$, set the derivative $k-\frac1x=0\Rightarrow x=\frac1k$, value $1+\ln k$; if $k\le0$ it is $-\infty$. Hence
$$L^*(u)=37u+3+\ln(5-u)+\ln(7-3u)+\ln(-4-12u)\quad\text{for }u<-\tfrac13\ (\text{all three arguments}>0),\ \text{else }-\infty .$$
The dual is a 1-D concave maximisation.

**(d) QP** $\min\frac12x^TQx+c^Tx$ s.t. $Ax\ge b$, $Q\succ0$: $L=u^Tb+(c-A^Tu)^Tx+\frac12x^TQx$, minimised at $Q\tilde x=-(c-A^Tu)$:
$$L^*(u)=-\tfrac12(c-A^Tu)^TQ^{-1}(c-A^Tu)+u^Tb,\qquad D:\ \max_{u\ge0}L^*(u).$$
(The same as the class QP dual with $A\to-A$, $b\to-b$.)

#### 3.4 Column geometry: the dual as a supporting hyperplane

For each $x\in X$, record the resource–cost vector $(s,z)=(g(x),f(x))\in\mathbb{R}^{m+1}$. Define
$$I=\{(s,z):\ \exists x\in X,\ s\ge g(x),\ z\ge f(x)\}.$$
(If $X$ is convex and $f,g_i$ convex, then $I$ is convex.) The primal value $z^*$ is the lowest $z$ on the vertical axis $s=0$ inside $I$.

Hyperplane $H_{u,\alpha}=\{(s,z):u^Ts+z=\alpha\}$ meets the axis $L=\{s=0\}$ at height $\alpha$. It is a **lower support** of $I$ if $u^Ts+z\ge\alpha$ on $I$. Finding the highest such intercept:
$$\max_{u,\alpha}\alpha\ \text{ s.t. }u^Ts+z\ge\alpha\ \forall(s,z)\in I\ =\ \max_{u\ge0,\alpha}\alpha\ \text{ s.t. }u^Tg(x)+f(x)\ge\alpha\ \forall x\in X$$
$$=\max_{u\ge0,\alpha}\alpha\ \text{ s.t. }L^*(u)\ge\alpha\ =\ \max_{u\ge0}L^*(u)=v^* .$$
(🔍 $u\ge0$ is forced: $I$ extends upward in every $s_i$, so a negative $u_i$ would make $u^Ts+z\to-\infty$.)
**The dual = the highest lower-supporting hyperplane of $I$; the primal = the lowest point of $I$ on the axis.** A gap appears when $I$ is non-convex near the axis.

#### 3.5 Concavity and weak duality

**Theorem 7.1: $L^*(u)$ is concave.** For $u=\lambda u_1+(1-\lambda)u_2$:
$$L^*(u)=\min_x\big[\lambda(f+u_1^Tg)+(1-\lambda)(f+u_2^Tg)\big]\ge\lambda\min_x(f+u_1^Tg)+(1-\lambda)\min_x(f+u_2^Tg)$$
(the min of a sum is at least the sum of the mins).

**Theorem 8.1 (weak duality):** if $\bar x$ is feasible for OP and $\bar u\ge0$, then $f(\bar x)\ge f(\bar x)+\bar u^Tg(\bar x)\ge\min_{x\in X}(f+\bar u^Tg)=L^*(\bar u)$. So $z^*\ge v^*$.
- **Cor 8.1:** $f(\bar x)=L^*(\bar u)$ ⇒ both are optimal.
- **Cor 8.2:** $z^*=-\infty$ ⇒ the dual is infeasible. **Cor 8.3:** $v^*=+\infty$ ⇒ the primal is infeasible.
- Unlike LP, strong duality does **not** hold for general nonlinear problems.

#### 3.6 Saddle-point optimality

$(\bar x,\bar u)$ with $\bar x\in X$, $\bar u\ge0$ is a **saddle point** of $L$ if
$$L(\bar x,u)\le L(\bar x,\bar u)\le L(x,\bar u)\qquad\forall x\in X,\ u\ge0$$
(min over $x$, max over $u$: a mountain pass).

**Theorem 9.1.** $(\bar x,\bar u)$ is a saddle point ⟺ (1) $L(\bar x,\bar u)=L^*(\bar u)$, (2) $g(\bar x)\le0$, (3) $\bar u^Tg(\bar x)=0$. ⟺ $\bar x,\bar u$ are optimal for OP and D with **no duality gap**.

*Proof sketch.* (⇒) (1) is the right inequality. The left inequality says $f(\bar x)+u^Tg(\bar x)\le f(\bar x)+\bar u^Tg(\bar x)$ for all $u\ge0$. If some $g_i(\bar x)>0$, a large $u_i$ breaks it, so (2). Putting $u=0$ gives $\bar u^Tg(\bar x)\ge0$, and with (2) and $\bar u\ge0$ this gives (3). Then $f(\bar x)=L(\bar x,\bar u)=L^*(\bar u)$: optimal with no gap (Cor 8.1). (⇐) Reverse the steps. (No gap ⇒) $L^*(\bar u)\le f(\bar x)+\bar u^Tg(\bar x)\le f(\bar x)$ with equal ends forces (1)–(3).

**Theorem 9.2:** for convex differentiable $f,g$ and $X=\mathbb{R}^n$, **KKT ⇒ saddle point ⇒ optimal, no gap.** (Stationarity of a convex $L(\cdot,\bar u)$ means $\bar x$ minimises it, i.e. condition (1).)

#### 3.7 Perturbation function and strong duality

Perturb the right-hand sides:
$$OP_y:\ z(y)=\min f(x)\ \text{ s.t. }g_i(x)\le y_i,\ x\in X;\qquad z^*=z(0),\quad Y=\{y:\exists x\in X,\ g(x)\le y\}.$$
**Lemma 10.1:** (convex case) $Y$ is convex and $z(y)$ is a **convex** function. 🔍 Mix near-optimal $x^1,x^2$ for $y^1,y^2$: $x^3=\lambda x^1+(1-\lambda)x^2$ is feasible for $y^3$ by convexity of $g$, and $f(x^3)\le\lambda f(x^1)+(1-\lambda)f(x^2)$.
Also $z(y)$ is **non-increasing** (a bigger $y$ is a looser constraint).

**Theorem 10.1 (Strong Duality).** $X$ open convex, $f,g_i$ convex, and $\partial z(0)\neq\emptyset$. Then: (1) D attains its optimum at some $\bar u\ge0$; (2) $\bar u$ is dual optimal ⟺ $-\bar u\in\partial z(0)$; (3) **$z^*=v^*$**; (4) primal optima $=S^1\cap S^2$ with $S^1=\{x\in X:g(x)\le0,\ \bar u^Tg(x)=0\}$ and $S^2=\{x:L(x,\bar u)=L^*(\bar u)\}$.

*Key steps of the proof.*
- Take $-\bar u\in\partial z(0)$: $z(y)\ge z(0)-\bar u^Ty$.
- $\bar u\ge0$: since $z$ is non-increasing, $z(0)\ge z(e_i)\ge z(0)-\bar u_i$, so $\bar u_i\ge0$.
- For any $x\in X$: $f(x)\ge z(g(x))\ge z(0)-\bar u^Tg(x)$, so $z(0)\le f(x)+\bar u^Tg(x)$. Taking the inf over $x$: $z^*\le L^*(\bar u)\le v^*$. With weak duality, $z^*=v^*$.
- **Interpretation:** $-\bar u$ is the slope of the perturbation function at 0, i.e. $\bar u_i$ is the **shadow price** of constraint $i$ (as in LP sensitivity).

#### 3.8 Stability and Slater

OP is **stable** if $\partial z(0)\neq\emptyset$. Equivalent (Lemma 10.2): $z(0)$ is finite and $\exists M>0$ with $\frac{z(0)-z(y)}{\|y\|}\le M$ for all $y\neq0$. In words, $z$ doesn't drop infinitely steeply away from 0. (Proved with the separating hyperplane theorem.)

**Lemma 10.3 (Slater ⇒ stable).** For $\min f$ s.t. $g(x)\le0$, $Ax=b$ (convex, $A$ full row rank) with a **Slater point** $x^0$ ($g(x^0)<0$, $Ax^0=b$) and finite $z^*$, the problem is stable. Hence **strong duality holds**. *Idea:* along $y=(\epsilon,\dots,\epsilon)$, $x^0$ stays feasible for small $|\epsilon|$, so $h(\epsilon)=z(\epsilon,\dots,\epsilon)$ is finite and convex around 0, has a subgradient $K$, and $M=|K|$ bounds the slope.

#### 3.9 Duality strategies

1. **Dualise "bad" constraints.** $\min c^Tx$ s.t. $Ax\le b$ (hard, e.g. side constraints), $Nx\le g$ (easy, e.g. network structure). Keep $X=\{Nx\le g\}$:
$$L^*(u)=-u^Tb+\min_{Nx\le g}(c+A^Tu)^Tx,$$
easy to evaluate for each $u$; maximise over $u\ge0$ (e.g. by subgradient ascent).
2. **Split a large problem into small ones.** With coupling constraints $B^1x^1+B^2x^2\le d$ and separate blocks $A^1x^1\le b^1$, $A^2x^2\le b^2$: dualise only the coupling constraint; $L^*(u)$ separates into one problem per block (this is the dual decomposition of the Dual Ascent block).
3. The reading also illustrates Lagrangian duality on discrete problems (with gaps) and introduces conic duality.

---

### 4. Summary

✅ Dual: $L=f+u^Tg$, $L^*(u)=\min_{x\in X}L$, $\max_{u\ge0}L^*$; sign rules $\le/\ge/=$ ⇒ $u\ge0/\le0$/free.
✅ Worked duals: LP, binary IP ($u^Tb+\sum\min\{0,(c-A^Tu)_j\}$), log-barrier, QP.
✅ Geometry: dual = highest lower support of $I=\{(g(x),f(x))+\text{up}\}$.
✅ $L^*$ concave; weak duality; saddle point ⟺ optimal pair with no gap.
✅ Strong duality ⟸ $\partial z(0)\neq\emptyset$ (stability) ⟸ Slater (convex case); $\bar u$ = shadow prices.

<!-- LECTURE END: lec18-duality-thy -->

---

<!-- LECTURE START: newton-equality -->
<a id="newton-equality"></a>

## Lecture 47: Newton's Method for Equality-Constrained Optimization

*Source file:* `Newton-equality.pdf`

### 1. What this lecture is about

Extend Newton's method to
$$\min f(x)\quad\text{s.t.}\quad Ax=b,\qquad x\in\mathbb{R}^n,$$
by minimising the quadratic Taylor model **subject to staying feasible**. Each step solves one linear **KKT system**.

---

### 2. Background you need first

- Unconstrained Newton: minimise $f(x^k)+\nabla f^Td+\frac12d^T\nabla^2fd$, giving $d=-(\nabla^2f)^{-1}\nabla f$ (Descent Methods block).
- Lagrange/KKT for equality constraints (Constrained Optimality blocks).
- Block matrices.

---

### 3. Main content with full derivation

1. **Quadratic model.** At the current $x^k$, a step $d$ gives
$$f(x^k+d)\approx f(x^k)+\nabla f(x^k)^Td+\tfrac12d^T\nabla^2f(x^k)d .$$
2. **Stay feasible.** Suppose $x^k$ is feasible: $Ax^k=b$. We want $A(x^k+d)=b$, i.e. $Ax^k+Ad=b$, which gives
$$\boxed{Ad=0}$$
(The step must lie in the null space of $A$.)
3. **Subproblem** (an equality-constrained QP in $d$):
$$\min_d\ f(x^k)+\nabla f(x^k)^Td+\tfrac12d^T\nabla^2f(x^k)d\quad\text{s.t.}\quad Ad=0\ \ [\lambda].$$
4. **Lagrangian** (write $\nabla f^k=\nabla f(x^k)$, $\nabla^2f^k=\nabla^2f(x^k)$):
$$\mathcal L(d,\lambda)=f(x^k)+\nabla f^{kT}d+\tfrac12d^T\nabla^2f^kd+\lambda^T(Ad).$$
5. **KKT conditions:** $\nabla_d\mathcal L=0$ and feasibility:
$$\nabla f^k+\nabla^2f^k\,d+A^T\lambda=0,\qquad Ad=0 .$$
🔍 ($\nabla_d(\lambda^TAd)=A^T\lambda$, and $\nabla_d(\frac12d^THd)=Hd$ for symmetric $H$.)
6. **As one linear system:**
$$\begin{pmatrix}\nabla^2f^k&A^T\\A&0\end{pmatrix}\begin{pmatrix}d\\\lambda\end{pmatrix}=\begin{pmatrix}-\nabla f^k\\0\end{pmatrix}.$$
Solve it to get $d^k$ and $\lambda$.
7. **Update:** $x^{k+1}=x^k+d^k$. It is still feasible ($Ad^k=0$). Repeat until $d^k\approx0$ (or use a line search $x^{k+1}=x^k+\alpha d^k$ for robustness far from the optimum).

**Remarks.**
- The KKT matrix is invertible if $A$ has full row rank and $\nabla^2f^k$ is positive definite on $\{d:Ad=0\}$.
- If $f$ is quadratic, the model is exact, so **one step reaches the optimum**.
- Compare the Yinyu Ye block: this is "Newton on the KKT equations" with a feasible start (there the sign convention is $-A^Ty$).

---

### 4. Worked example

$\min x_1^2+2x_2^2$ s.t. $x_1+x_2=3$, starting from the feasible $x^0=(3,0)$.
- $\nabla f=(2x_1,4x_2)=(6,0)$; $\nabla^2f=\operatorname{diag}(2,4)$; $A=(1\ 1)$.
- System:
$$\begin{pmatrix}2&0&1\\0&4&1\\1&1&0\end{pmatrix}\begin{pmatrix}d_1\\d_2\\\lambda\end{pmatrix}=\begin{pmatrix}-6\\0\\0\end{pmatrix}.$$
- Row 3: $d_2=-d_1$. Row 2: $-4d_1+\lambda=0\Rightarrow\lambda=4d_1$. Row 1: $2d_1+4d_1=-6\Rightarrow d_1=-1$, so $d_2=1$, $\lambda=-4$.
- $x^1=(3,0)+(-1,1)=(2,1)$.

Check with Lagrange directly: $\nabla f=\mu\nabla h$ ⇒ $2x_1=4x_2$ ⇒ $x_1=2x_2$; with $x_1+x_2=3$: $(2,1)$ ✓. One step, as expected for a quadratic.

---

### 5. Summary

✅ Feasible start + $Ad=0$ keeps every iterate feasible.
✅ Newton step = minimiser of the quadratic model on $\{Ad=0\}$ ⇒ solve $\begin{pmatrix}\nabla^2f&A^T\\A&0\end{pmatrix}\binom d\lambda=\binom{-\nabla f}{0}$.
✅ $x^{k+1}=x^k+d^k$; exact in one step for quadratics.

<!-- LECTURE END: newton-equality -->

---

<!-- LECTURE START: lecture-22 -->
<a id="lecture-22"></a>

## Lecture 48: ADMM: Alternating Direction Method of Multipliers (reading: R. Tibshirani, CMU 10-725, Lecture 22)

*Source file:* `Lecture_22.pdf`

### 1. What this lecture is about

1. Recap: **dual ascent**, and why it needs strong assumptions on $f$
2. **Augmented Lagrangian / method of multipliers (MM):** add $\frac\rho2\|Ax-b\|^2$ to make things well behaved
3. **ADMM:** combine the decomposability of dual ascent with the robustness of MM by splitting variables into $x$ and $z$ and updating them **alternately**
4. **Scaled form** of ADMM (with a full derivation by completing the square)
5. Examples: lasso-type problems via **soft-thresholding**, group lasso, sparse + low-rank decomposition, **consensus ADMM** for distributed learning, and why the choice of splitting matters (2-D fused lasso)

---

### 2. Background you need first

- Dual ascent (Dual Ascent block): $x^{k}=\arg\min_xL(x,u^{k-1})$, $u^{k}=u^{k-1}+t_k(Ax^k-b)$.
- $\|v\|_2^2=v^Tv$; completing the square: $a^Tv+\frac\rho2\|v\|^2=\frac\rho2\|v+a/\rho\|^2-\frac{1}{2\rho}\|a\|^2$.
- **Lasso** penalty $\lambda\|x\|_1$ encourages sparse solutions.
- **Proximal operator** $\operatorname{prox}_{t g}(v)=\arg\min_z\ g(z)+\frac1{2t}\|z-v\|^2$ ("a gradient step for non-smooth $g$").

---

### 3. Main content with full derivations

#### 3.1 Dual ascent recap and its weakness

$\min f(x)$ s.t. $Ax=b$; $L(x,u)=f(x)+u^T(Ax-b)$.
$$x^{(k)}=\arg\min_xL(x,u^{(k-1)}),\qquad u^{(k)}=u^{(k-1)}+t_k(Ax^{(k)}-b).$$
- ✅ If $f$ decomposes over blocks, the $x$-step decomposes (distributed).
- ❌ Primal convergence needs strong assumptions (e.g. strictly convex, closed $f$). For many real problems (e.g. linear $f$), $\arg\min_xL$ is unbounded or not unique.

#### 3.2 Augmented Lagrangian (method of multipliers)

For fixed $\rho>0$, solve the equivalent problem $\min f(x)+\frac\rho2\|Ax-b\|_2^2$ s.t. $Ax=b$. (The extra term is zero on the feasible set, so the optimum is unchanged.)
$$L_\rho(x,u)=f(x)+u^T(Ax-b)+\tfrac\rho2\|Ax-b\|_2^2 .$$
$$x^{(k)}=\arg\min_xL_\rho(x,u^{(k-1)}),\qquad u^{(k)}=u^{(k-1)}+\rho(Ax^{(k)}-b).$$
- ✅ The quadratic term makes the $x$-problem strongly convex (if $A$ has full column rank), so convergence is much more robust. The step size is simply $\rho$.
- ❌ $\|Ax-b\|^2$ couples all variables, so separability of $f$ **no longer** gives separable updates.

🔍 *Why step size $\rho$:* the optimality of $x^{(k)}$ gives $0=\nabla f(x^{(k)})+A^T\big(u^{(k-1)}+\rho(Ax^{(k)}-b)\big)=\nabla f(x^{(k)})+A^Tu^{(k)}$. So with this update, $(x^{(k)},u^{(k)})$ always satisfies the stationarity condition of the original Lagrangian; only feasibility remains to be achieved.

#### 3.3 ADMM

Assume the problem has (or is rewritten into) the form
$$\min_{x,z}\ f(x)+g(z)\quad\text{s.t.}\quad Ax+Bz=c .\tag{22.1}$$
$$L_\rho(x,z,u)=f(x)+g(z)+u^T(Ax+Bz-c)+\tfrac\rho2\|Ax+Bz-c\|_2^2 .$$
**Updates:**
$$x^{(k)}=\arg\min_xL_\rho(x,z^{(k-1)},u^{(k-1)})$$
$$z^{(k)}=\arg\min_zL_\rho(x^{(k)},z,u^{(k-1)})\quad(\text{uses the newest }x)$$
$$u^{(k)}=u^{(k-1)}+\rho(Ax^{(k)}+Bz^{(k)}-c)$$
- Updating $x$ and $z$ **jointly** would be MM. Alternating them (one Gauss–Seidel pass) is what makes each subproblem simple: often only $f$ or only $g$ appears.
- Convergence (under mild conditions, Boyd et al. 2010): residual $r^{(k)}=Ax^{(k)}+Bz^{(k)}-c\to0$; objective $f(x^{(k)})+g(z^{(k)})\to f^*+g^*$; $u^{(k)}\to u^*$.
- Like other first-order methods: a moderately accurate answer in a few dozen iterations, high accuracy slowly. **$\rho$ matters:** too large under-weights the objective, too small under-weights feasibility.

#### 3.4 Scaled ADMM (derivation)

Let $w=u/\rho$ (scaled dual variable) and $r=Ax+Bz-c$. Complete the square:
$$u^Tr+\tfrac\rho2\|r\|^2=\tfrac\rho2\big(\|r\|^2+2w^Tr\big)=\tfrac\rho2\big(\|r+w\|^2-\|w\|^2\big).$$
So $L_\rho=f(x)+g(z)+\frac\rho2\|Ax+Bz-c+w\|_2^2-\frac\rho2\|w\|_2^2$, and since the last term does not involve $x$ or $z$:
$$x^{(k)}=\arg\min_x\ f(x)+\tfrac\rho2\|Ax+Bz^{(k-1)}-c+w^{(k-1)}\|_2^2$$
$$z^{(k)}=\arg\min_z\ g(z)+\tfrac\rho2\|Ax^{(k)}+Bz-c+w^{(k-1)}\|_2^2$$
$$w^{(k)}=w^{(k-1)}+Ax^{(k)}+Bz^{(k)}-c$$
($w$ = running sum of residuals.)

#### 3.5 Example: lasso (the template for the others)

🔍 *Worked out here as the simplest case.* $\min_\beta\frac12\|y-X\beta\|^2+\lambda\|\beta\|_1$. Split $\alpha=\beta$: $f(\beta)=\frac12\|y-X\beta\|^2$, $g(\alpha)=\lambda\|\alpha\|_1$, constraint $\beta-\alpha=0$.
- $\beta$-step: $\min_\beta\frac12\|y-X\beta\|^2+\frac\rho2\|\beta-\alpha+w\|^2$. Gradient: $-X^T(y-X\beta)+\rho(\beta-\alpha+w)=0$ ⇒ $\beta^+=(X^TX+\rho I)^{-1}\big(X^Ty+\rho(\alpha-w)\big)$. ($X^TX+\rho I$ is always invertible.)
- $\alpha$-step: $\min_\alpha\lambda\|\alpha\|_1+\frac\rho2\|\beta^+-\alpha+w\|^2$ separates by coordinate, and its solution is **soft-thresholding**:
$$\alpha_j^+=S_{\lambda/\rho}(\beta_j^++w_j),\qquad S_t(v)=\operatorname{sign}(v)\max(|v|-t,0).$$
(1-D: minimise $t|a|+\frac12(a-v)^2$. For $a>0$ the derivative $t+a-v=0$ gives $a=v-t$, valid if $v>t$; similarly $a=v+t$ if $v<-t$; otherwise $a=0$.)
- $w^+=w+\beta^+-\alpha^+$.

#### 3.6 Group lasso

$\min_\beta\frac12\|y-X\beta\|_2^2+\lambda\sum_{g=1}^Gc_g\|\beta_g\|_2$ (whole groups of coefficients switched on/off). Split $\beta-\alpha=0$:
- $\beta^+=(X^TX+\rho I)^{-1}(X^Ty+\rho(\alpha-w))$ (factor once by Cholesky in $O(p^3)$; afterwards each update costs $O(p^2)$);
- $\alpha_g^+=R_{c_g\lambda/\rho}(\beta_g^++w_g)$, with **group soft-thresholding** $R_t(x)=\big(1-\frac{t}{\|x\|_2}\big)_+x$ (shrink the whole block toward 0, or set it to 0 if $\|x\|_2\le t$);
- $w^+=w+\beta^+-\alpha^+$.

This works for any sum of norms whose proximal operator is known.

#### 3.7 Sparse subspace estimation and sparse + low-rank decomposition

- **Sparse PCA (Fantope):** $\max\operatorname{tr}(SY)-\lambda\|Y\|_1$ s.t. $Y\in\mathcal F_k=\{0\preceq Y\preceq I,\ \operatorname{tr}Y=k\}$. Split $Y=Z$: $Y^+=P_{\mathcal F_k}(Z-W+S/\rho)$ (clip eigenvalues: $\sigma_i\mapsto\min\{\max\{\sigma_i-\theta,0\},1\}$ with $\theta$ chosen so they sum to $k$); $Z^+=S_{\lambda/\rho}(Y^++W)$; $W^+=W+Y^+-Z^+$.
- **Sparse + low rank:** $\min\|L\|_{tr}+\lambda\|S\|_1$ s.t. $L+S=M$. $L^+=$ matrix soft-thresholding of the singular values of $M-S+W$; $S^+=S_{\lambda/\rho}(M-L^++W)$ elementwise; $W^+=W+M-L^+-S^+$. On video, $L$ captures the static background and $S$ the moving objects.

#### 3.8 Consensus ADMM (distributed learning)

$\min_x\sum_{i=1}^Bf_i(x)$ (e.g. the loss on data shard $i$). Give each machine its own copy and force agreement:
$$\min_{x_1,\dots,x_B,x}\sum_if_i(x_i)\ \text{ s.t. }\ x_i=x,\ i=1..B .$$
$$x_i^+=\arg\min_{x_i}f_i(x_i)+\tfrac\rho2\|x_i-x+w_i\|^2\quad(\text{in parallel}),\qquad x^+=\tfrac1B\sum_i(x_i^++w_i),\qquad w_i^+=w_i+x_i^+-x^+ .$$
🔍 *Simplification:* averaging the $w$-update gives $\bar w^+=\bar w+\bar x^+-x^+$, and $x^+=\bar x^++\bar w$, so $\bar w^+=0$ after the first iteration. Then $x^+=\bar x^+$ (just the average):
$$x_i^+=\arg\min f_i(x_i)+\tfrac\rho2\|x_i-\bar x+w_i\|^2,\qquad w_i^+=w_i+x_i^+-\bar x^+ .$$
Intuition: each machine fits its own data with an $\ell_2$ pull toward the average; if $x_i$ is above the average, $w_i$ grows and pulls it back harder next time.

**General consensus** $\min\sum_if_i(a_i^Tx+b_i)+g(x)$: the same, but the $x$-update becomes a prox of $g$: $x^+=\arg\min\frac{B\rho}{2}\|x-\bar x^+-\bar w\|^2+g(x)$, and $\bar w$ no longer vanishes.

#### 3.9 Choosing the decomposition matters (2-D fused lasso)

ADMM resembles block coordinate descent: it is slow if the blocks are highly correlated. For image denoising $\min\frac12\|y-\theta\|^2+\lambda\|D\theta\|_1$ ($D$ = differences between neighbouring pixels):
- Split $z=D\theta$: each step is a sparse linear solve plus soft-thresholding, $O(n)$ per iteration, but **slow convergence**.
- Split into horizontal ($H$) and vertical ($V$) copies with $H=V$: each update becomes many independent **1-D fused-lasso** problems (fast exact solvers), and convergence is **much faster**.

---

### 4. Summary

✅ Dual ascent decomposes but is fragile; the augmented Lagrangian is robust but couples variables.
✅ ADMM on $\min f(x)+g(z)$, $Ax+Bz=c$: $x$-step, $z$-step (with the new $x$), $u\mathrel{+}=\rho\cdot$residual.
✅ Scaled form with $w=u/\rho$ (complete the square).
✅ Lasso/group lasso ⇒ ridge-type solve + (group) soft-thresholding; consensus ADMM ⇒ parallel local solves + averaging.
✅ $\rho$ and the choice of splitting strongly affect speed.

<!-- LECTURE END: lecture-22 -->

---

<!-- LECTURE START: dualmethods -->
<a id="dualmethods"></a>

## Lecture 49: Dual Methods (full notes): dual ascent → method of multipliers → ADMM

*Source file:* `DualMethods.pdf`

> 📌 Part I of `DualMethods.pdf` (dual ascent and dual decomposition) is identical to `DualAscent.pdf` and is explained step by step in the block **"Algorithms Using the Dual: dual ascent and dual decomposition"**. This block explains **Parts II and III**, which are new: the method of multipliers and ADMM.

### 1. What this lecture is about

1. Quick recap of dual ascent on $\min f(x)$ s.t. $Ax=b$
2. **Method of multipliers:** add the penalty $\frac\rho2\|Ax-b\|^2$, run dual ascent on the new problem, and see why the step size is exactly $\rho$ (full derivation)
3. **ADMM:** split into $\min f(x)+g(z)$ s.t. $Ax+Bz=b$ and update $x$, then $z$, then $\lambda$

---

### 2. Background you need first

- Dual ascent: $x^{k+1}=\arg\min_xf(x)+\lambda^{kT}(Ax-b)$; $\lambda^{k+1}=\lambda^k+\alpha(Ax^{k+1}-b)$, since $\nabla\theta(\lambda^k)=Ax^{k+1}-b$.
- KKT for $\min f$ s.t. $Ax=b$: $\nabla f(x)+A^T\lambda=0$ and $Ax=b$.
- $\nabla_x\|Ax-b\|^2=2A^T(Ax-b)$. (From $(Ax-b)^T(Ax-b)=x^TA^TAx-2b^TAx+b^Tb$, with gradient $2A^TAx-2A^Tb$.)

---

### 3. Main content with full derivations

#### 3.1 Part I recap: dual ascent

$p^*=\min f(x)$ s.t. $Ax=b$, $A\in\mathbb{R}^{m\times n}$. $\theta(\lambda)=\min_x f(x)+\lambda^T(Ax-b)$; dual $\max_{\lambda\in\mathbb{R}^m}\theta(\lambda)$.
Loop: $x^{k+1}=\arg\min_xL(x,\lambda^k)$, $\lambda^{k+1}=\lambda^k+\alpha(Ax^{k+1}-b)$; stop when $\|\nabla\theta\|=\|Ax-b\|<\varepsilon$ (e.g. $10^{-9}$). If $f=\sum_if_i(x^{(i)})$ is separable, the $x$-step splits into $B$ independent problems (dual decomposition).

#### 3.2 Part II: method of multipliers

**Two equivalent problems:**
$$P1:\ \min_xf(x)\ \text{ s.t. }Ax=b\ [\lambda],\qquad P2:\ \min_xf(x)+\tfrac\rho2\|Ax-b\|^2\ \text{ s.t. }Ax=b .$$
$x$ is feasible for P1 ⟺ feasible for P2, and on feasible points $f(x)=f(x)+\frac\rho2\underbrace{\|Ax-b\|^2}_{=0}$. Same feasible set, same objective values ⇒ same solutions.

**Dual ascent on P2.** The Lagrangian of P2 is the **augmented Lagrangian** of P1:
$$L_\rho(x,\lambda)=f(x)+\lambda^T(Ax-b)+\tfrac\rho2\|Ax-b\|^2=f(x)+\lambda^T(Ax-b)+\tfrac\rho2\big(x^TA^TAx-2b^TAx+b^Tb\big).$$
```
Init λ^0, k = 0
Loop till convergence {
    x^{k+1} = argmin_{x∈R^n} L_ρ(x, λ^k)
    λ^{k+1} = λ^k + ρ (A x^{k+1} − b)        # step size = ρ
}
```
**Why the step size is $\rho$ (derivation).** $x^{k+1}$ minimises $L_\rho(\cdot,\lambda^k)$, so
$$\nabla_xL_\rho(x^{k+1},\lambda^k)=0:\quad\nabla f(x^{k+1})+A^T\lambda^k+\rho\big(A^TAx^{k+1}-A^Tb\big)=0 .$$
Factor $A^T$ out of the last two terms:
$$\nabla f(x^{k+1})+A^T\underbrace{\big(\lambda^k+\rho(Ax^{k+1}-b)\big)}_{=\lambda^{k+1}}=0\quad\Longrightarrow\quad\boxed{\nabla f(x^{k+1})+A^T\lambda^{k+1}=0}$$
Compare with the KKT conditions of P1: $\nabla f(x)+A^T\lambda=0$ and $Ax=b$. With step size $\rho$, **every iterate $(x^{k+1},\lambda^{k+1})$ satisfies the first KKT condition (stationarity) automatically**. Only primal feasibility $Ax=b$ is left, and that is what the iterations drive to zero.

**Pros/cons.** The penalty makes the $x$-subproblem well behaved (strongly convex when $A$ has full column rank), so convergence needs far weaker assumptions than dual ascent. But if $f(x)=\sum_if_i(x^{(i)})$, the term $\|Ax-b\|^2=\|\sum_iA_ix^{(i)}-b\|^2$ couples the blocks (cross terms $x^{(i)T}A_i^TA_jx^{(j)}$), so the $x$-step **no longer decomposes**. ADMM fixes this.

#### 3.3 Part III: ADMM (Alternating Direction Method of Multipliers)

$$\min f(x)+g(z)\quad\text{s.t.}\quad Ax+Bz=b .$$
$$L_\rho(x,z,\lambda)=f(x)+g(z)+\lambda^T(Ax+Bz-b)+\tfrac\rho2\|Ax+Bz-b\|^2 .$$
```
Initialize λ^0, z^0, k = 0
Loop {
    x^{k+1} = argmin_x  L_ρ(x, z^k, λ^k)            # z frozen at its old value
    z^{k+1} = argmin_z  L_ρ(x^{k+1}, z, λ^k)        # uses the NEW x
    λ^{k+1} = λ^k + ρ (A x^{k+1} + B z^{k+1} − b);  k = k+1
}
```
🔍 *Why it helps:* in the $x$-step only $f(x)$ and a quadratic in $x$ appear; in the $z$-step only $g(z)$ and a quadratic in $z$. Each is often easy (closed form, or a proximal step such as soft-thresholding). Minimising over $(x,z)$ jointly would be the method of multipliers. Alternating (one pass of block coordinate descent) keeps the subproblems simple while keeping the robustness of the augmented Lagrangian.

🔍 *Same stationarity trick for $z$:* optimality of $z^{k+1}$ gives $\nabla g(z^{k+1})+B^T\big(\lambda^k+\rho(Ax^{k+1}+Bz^{k+1}-b)\big)=\nabla g(z^{k+1})+B^T\lambda^{k+1}=0$. So the $z$-part of KKT holds at every iteration; the $x$-part and feasibility hold in the limit.

Examples, the scaled form and consensus/distributed ADMM are worked out in the **ADMM (Lecture 22)** block.

---

### 4. Worked example (method of multipliers)

$\min\frac12x^2$ s.t. $x=1$ ($A=1$, $b=1$), $\rho=1$, $\lambda^0=0$.
- $x$-step: $\min\frac12x^2+\lambda(x-1)+\frac12(x-1)^2$ ⇒ $x+\lambda+(x-1)=0$ ⇒ $x=\frac{1-\lambda}{2}$.
- $k=0$: $x^1=0.5$, $\lambda^1=0+(0.5-1)=-0.5$.
- $k=1$: $x^2=0.75$, $\lambda^2=-0.5-0.25=-0.75$.
- $k=2$: $x^3=0.875$, $\lambda^3=-0.875$. …
Converges to $x^*=1$, $\lambda^*=-1$ (check KKT: $x+\lambda=0$ ✓). At every step $x^{k+1}+\lambda^{k+1}=0$ holds exactly, as derived.

---

### 5. Summary

✅ P2 adds $\frac\rho2\|Ax-b\|^2$; it has the same solutions as P1.
✅ Method of multipliers: $x^{k+1}=\arg\min L_\rho$, $\lambda^{k+1}=\lambda^k+\rho(Ax^{k+1}-b)$; step $\rho$ makes $\nabla f+A^T\lambda=0$ hold at every iterate.
✅ ADMM: minimise over $x$, then $z$ (with the new $x$), then dual update. Robust like MM, decomposable like dual ascent.

<!-- LECTURE END: dualmethods -->

---

<!-- APPENDIX START -->
## Appendix: Assignments & problem sets

- `Assignment1.pdf`, `Assignment2.pdf`, `ProblemSet3.pdf` — not covered in these notes (ask if you want worked solutions).
<!-- APPENDIX END -->
