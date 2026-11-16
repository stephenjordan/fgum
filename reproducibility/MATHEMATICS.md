# Mathematical checks implemented by the certificate

The supported family is the fifteen integer pairs $3\le k<D\le8$.
The notation in this file follows Appendix A. Logarithms are natural.

## Balance equation and positive branch

Write

$$
\varphi(z)=1-(1-z)^{D-1},\qquad
\pi_D(z)=1-(1-z)^D-Dz(1-z)^{D-1}.
$$

For $0<y<1$, set

$$
p=1-(1-y)^{1/(D-1)},\quad u=y^{k-1},\quad
a=\varphi(u),\quad e=\frac{p}{ay^{k-2}},
$$

$$
c= Dpy,\qquad
\mathcal B=e\pi_D(u)+(k-1)\pi_D(p)-c.
$$

The program proposes a positive root of $\mathcal B=0$, constructs an exact
rational bracket of radius $10^{-45}$, and verifies its endpoint signs with
interval arithmetic. It also verifies $e'(y)>0$ and

$$
\mathcal B'(y)=e'(y)[\pi_D(u)-Dua]<0
$$

throughout the bracket.

### Why this is the unique positive balance on the stable nonzero branch

This branch selection has an analytic justification independent of a root
finding heuristic. Let $r=D-1$ and

$$
H(z)=\frac{z\varphi'(z)}{\varphi(z)}.
$$

Writing $q=1-z$ gives

$$
H(z)=\frac{r}{1+q^{-1}+\cdots+q^{-(r-1)}}.
$$

Thus $H$ strictly decreases from one to zero. Parameterize the fixed points
by $p$, with $y=\varphi(p)$ and $u=y^{k-1}$. For the Jacobian
$J=\begin{pmatrix}0&A\\B&C\end{pmatrix}$ of the two-message recursion, define

$$
S(p)=C+AB=H(p)[k-2+(k-1)H(u)].
$$

It strictly decreases from $2k-3$ to zero. Since

$$
\rho(J)<1\ \Longleftrightarrow\ S(p)<1,\qquad
\frac{p e'(p)}{e(p)}=1-S(p),
$$

$e(p)$ has exactly one minimum and then increases along its unique stable
nonzero branch. The balance derivative satisfies

$$
\mathcal B'(p)=e'(p)[\pi_D(u)-Du\varphi(u)].
$$

The bracket equals $1-(1-u)^D-Du<0$. Also

$$
\mathcal B(p)\sim\frac{D(D-1)(k-2)}2p^2
\quad(p\downarrow0),\qquad \mathcal B(1)=k-D<0.
$$

The balance first increases from zero, then strictly decreases, crossing zero
exactly once on the stable nonzero branch. This is the largest positive fixed
point at its erasure probability. The zero fixed point is separate.

The program additionally certifies the whole upper branch numerically.
With $x=(1-y)^{1/(D-1)}$ and $Q(x)=a(y)y^{k-2}$, one has

$$
e=\frac{1-x}{Q(x)},\qquad e_x=-\frac{T(x)}{Q(x)^2},
\qquad T(x)=Q(x)+(1-x)Q_x(x).
$$

$T$ is a polynomial, regular at $x=0$. An interval cover certifies $T>0$
from zero to an outward upper enclosure of the balance value of $x$.

## Profiles, variances, and stability

For $d=2,\ldots,D$, the limiting row profiles are

$$
\rho_{A,d}=e\binom Dd u^d(1-u)^{D-d},\qquad
\rho_{B,d}=\binom Dd p^d(1-p)^{D-d}.
$$

Both incidence sums $\sum_d d\rho_{A,d}$ and $\sum_d d\rho_{B,d}$ equal
$c$. For

$$
g_d(z)=\sum_{\substack{0\le j\le d\\j\ {\rm even}}}\binom dj z^j,
\qquad L_\rho(t)=\sum_d\rho_d\log g_d(e^t),
$$

one has $L_A'(0)=L_B'(0)=c/2$. The variances are

$$
V_A=L_A''(0)=c/4+\rho_{A,2}/2,\qquad
V_B=L_B''(0)=c/4+\rho_{B,2}/2.
$$

The degree-two term is treated separately: its variance is one.
At nonzero tilt the program uses the positive pair-sum variance formula,
with weights $w_j=\binom dj e^{jt}$:

$$
(\log g_d(e^t))''
=\frac{\sum_{i<j}w_iw_j(i-j)^2}{(\sum_jw_j)^2}.
$$

This avoids subtracting nearly equal moments. Stability is checked with
the positive eigenvalue

$$
\rho(J)=\frac{C+\sqrt{C^2+4AB}}2<1.
$$

## Coverage of support fractions

The exponent is

$$
\Phi(\alpha)=\psi_A(c\alpha)+(k-1)\psi_B(c\alpha)
 -(k-1)c h(\alpha),
\qquad
\psi_\rho(z)=\inf_t[L_\rho(t)-zt].
$$

Every part of the feasible interval is covered by one of the following
arguments.

### Small support

For $\delta=10^{-4}$, the program checks

$$
\beta_{k,D}=
\exp(1)(D-1)^{k/2}\delta^{k/2-1}(1-\delta)^{-k/2}<1.
$$

This is the finite-support estimate used in the paper. The following
additional scalar check establishes the literal inequality $\Phi<0$ near
zero, without a grid. For $z\ge0$,

$$
\binom d{2j}\le\frac{[d(d-1)/2]^j}{j!},
\qquad
g_d(z)\le\exp[d(D-1)z^2/2].
$$

Using $\sum_d d\rho_d=c$ and choosing $z^2=\alpha/(D-1)$ gives

$$
\psi_\rho(c\alpha)
\le \frac{c\alpha}{2}[1+\log(D-1)-\log\alpha].
$$

Since $h(\alpha)\ge-\alpha\log\alpha$,

$$
\Phi(\alpha)\le c\alpha
\left[\frac{k}{2}(1+\log(D-1))
 +\left(\frac{k}{2}-1\right)\log\alpha\right].
$$

The bracket increases with $\alpha$. Its value at $\delta$ is
interval-certified negative for every supported pair.

### Exterior intervals and large-support cutoff

Let $\omega_A$ and $\omega_B$ be the sums of the odd-degree profile entries.
The odd-degree argument gives

$$
\alpha_{\max}=1-\frac{\max(\omega_A,\omega_B)}{c}<1.
$$

The program uses an outward upper enclosure of this cutoff, rounded upward
and padded by $10^{-3}$. It checks that the padded endpoint is below one.

For each cell covering
$[\delta,0.4]\cup[0.6,\overline\alpha_{\max}+10^{-3}]$,
it proposes a single pair of fixed tilts and evaluates

$$
U(\alpha)=L_A(t_A)+(k-1)L_B(t_B)
 -c\alpha[t_A+(k-1)t_B]-(k-1)c h(\alpha)
$$

at both endpoints with directed intervals. The same tilts are used at
both endpoints. Since $U''=(k-1)c/[\alpha(1-\alpha)]>0$, a negative upper
bound at both endpoints certifies the entire cell. The recorded cells are
contiguous and cover both ranges exactly.

### Central interval

For $0.4\le\alpha\le0.6$, write $v=c(\alpha-1/2)$ and choose
$t_A=v/V_A$, $t_B=v/V_B$. Here $v$ is the variable called $u$ in the
central-interval proof in Appendix A; the different letter avoids confusing
it with the earlier probability $u=y^{k-1}$ in this note.

Let $\Gamma(v)$ be the resulting Chernoff bound. At exact core balance,
$\Gamma(0)=\Gamma'(0)=0$ analytically. The program certifies

$$
\Gamma''(v)=
\frac{L_A''(v/V_A)}{V_A^2}-\frac2{V_A}
 +(k-1)\left[\frac{L_B''(v/V_B)}{V_B^2}-\frac2{V_B}\right]
 +\frac{k-1}{c\alpha(1-\alpha)}<0
$$

over a complete interval cover of $[0.4,0.6]$. Thus $\Phi\le\Gamma<0$
away from $\alpha=1/2$. The exact zero at one half follows from the
balance identity, rather than from rounding a numerical value to zero.

## Scores and exported data

The score calculation uses the certified root together with

$$
\alpha_{\min}=(1-e_{\rm core})(1-2^{1-D}),\qquad
s=1-\alpha_{\min}\hat I_D^*/D.
$$

The Hamming-weight average $\hat I_D^*$ is calculated as an exact rational
binomial sum, including the half-weight middle layer when $D$ is even.
`generate_tables.py` checks the displayed threshold, beta, and score precision,
and that the computed interval upper bounds imply the manuscript's displayed
conservative bounds.

Every decimal enclosure is rounded outward. The JSON files retain exact
rational cover endpoints, proposed tilts, root brackets, all interval results,
and the parameters required for replay.
