# Exact \(k=8\) target geometry in moment space

**Status:** exact moment-cone and Christoffel algebra. This note does not supply the missing arithmetic evaluations of \(m_7,m_8\); it tells us exactly what those evaluations need to achieve.

The exact moments already available are

\[
m_0=1,\quad m_1=1,\quad m_2=\frac43,\quad m_3=2,\quad
m_4=\frac{13}{4},\quad
m_5=\frac{101}{18},\quad
m_6=\frac{12809}{1260}.
\]

We now treat \(m_7,m_8\) as variables.

## 1. Exact Stieltjes lower bound for \(m_7\)

For a positive measure on \([0,\infty)\), the shifted Hankel matrix

\[
H_3^{(1)}=(m_{i+j+1})_{0\le i,j\le3}
\]

must be positive semidefinite. With the exact \(m_1,\dots,m_6\) above,

\[
\det H_3^{(1)}
=
\frac{1352400\,m_7-25866469}{76204800}.
\]

Therefore

\[
\boxed{
m_7\ge m_7^*
:=
\frac{25866469}{1352400}
=
19.126345016267\ldots
}
\]

This is a stronger exact floor than the approximate \(19.123\) value quoted in the reference outlook.

## 2. Moment-cone floor for \(m_8\)

The ordinary \(5\times5\) Hankel matrix

\[
H_4=(m_{i+j})_{0\le i,j\le4}
\]

must also be positive semidefinite.

Its determinant vanishes at

\[
\boxed{
m_8^{\min}(m_7)
=
\frac{
13335840000\,m_7^2
-504286322400\,m_7
+4794396454037
}{
748818000
}.
}
\]

Thus every admissible moment sequence satisfies

\[
m_8\ge m_8^{\min}(m_7).
\]

Write

\[
g:=m_7-m_7^*\ge0,
\qquad
\delta:=m_8-m_8^{\min}(m_7)\ge0.
\]

The variables \(g,\delta\) measure how far the moment sequence lies from the two relevant Stieltjes moment-cone boundaries.

## 3. Exact degree-four Christoffel function

The degree-four Christoffel origin mass is

\[
\lambda_4(0)
=
\frac{1}{e_0^T H_4^{-1}e_0}.
\]

Direct exact determinant reduction gives

\[
\lambda_4(0)
=
\frac{
13335840000m_7^2
-504286322400m_7
-748818000m_8
+4794396454037
}{
4\left(
24004512000m_7^2
-903891063600m_7
-1837779300m_8
+8574904106497
\right)
}.
\]

The expression becomes far more informative in \((g,\delta)\) coordinates:

\[
\boxed{
\lambda_4(0)
=
\frac{39243610000\,\delta}{
385252994000\,\delta
+
(1352400m_7-25866469)^2
}.
}
\]

Since

\[
\frac{39243610000}{385252994000}
=
\frac{1415}{13891}
=
\lambda_3(0),
\]

we may write

\[
\boxed{
\lambda_4(0)
=
\lambda_3(0)\,
\frac{\delta}{
\delta+c\,g^2
},
}
\]

where

\[
\lambda_3(0)=\frac{1415}{13891},
\qquad
c=\frac{18663120}{3931153}
=4.7474926567\ldots.
\]

This formula is the clean target geometry for the \(k=8\) programme.

### Interpretation

- If \(g=0\), i.e. \(m_7\) sits exactly at its Stieltjes floor, then the fourth degree adds no improvement over the six-moment certificate.
- If \(g>0\), improvement is controlled by how small the \(m_8\) excess \(\delta\) is compared with \(g^2\).
- The consumption side therefore does not require exact \(m_8\) in principle: a sufficiently strong **upper bound** on \(\delta\) is enough.

## 4. Exact target inequalities

Let the desired simple-critical proportion be \(S\). Under the same spectral counting conversion,

\[
1-2\lambda_4(0)\ge S
\]

is equivalent to

\[
\lambda_4(0)\le L:=\frac{1-S}{2}.
\]

For \(L<\lambda_3(0)\), this is exactly

\[
\boxed{
\delta
\le
\frac{L\,c}{\lambda_3(0)-L}\,g^2.
}
\]

Some useful targets are:

### 80%

For \(S=0.8\), \(L=0.1\),

\[
\boxed{
\delta
\le
\frac{2666160}{10471}\,g^2
=
254.6232451\ldots\,g^2.
}
\]

### 80.2%

For \(S=0.802\), \(L=0.099\),

\[
\delta
\le
164.0771689\ldots\,g^2.
\]

### 81%

For \(S=0.81\), \(L=0.095\),

\[
\delta
\le
65.70190286\ldots\,g^2.
\]

### 90%

For \(S=0.9\), \(L=0.05\),

\[
\delta
\le
4.576821466\ldots\,g^2.
\]

The 90% line illustrates why higher moments will eventually be necessary: the \(m_7,m_8\) pair would have to sit implausibly close to the \(H_4\) boundary.

## 5. Direct 80% upper bound for \(m_8\)

Equivalently, 80% is achieved whenever

\[
\boxed{
m_8
\le
B_{80}(m_7)
:=
\frac{
18670176000m_7^2
-713649484800m_7
+6822174057191
}{
68531400
}.
}
\]

Examples:

- at \(m_7=19.130\),
  \[
  m_8^{\min}=37.04741582\ldots,
  \qquad
  B_{80}=37.05081731\ldots;
  \]
- at \(m_7=19.136\),
  \[
  m_8^{\min}=37.09567533\ldots,
  \qquad
  B_{80}=37.11941098\ldots.
  \]

So an \(m_7\) value in the upper part of the reconnaissance band gives noticeably more \(m_8\) slack.

## 6. Research consequence

The \(k=8\) programme should be split into two arithmetic tasks:

1. prove \(m_7\) is strictly above
   \[
   m_7^*=25866469/1352400,
   \]
   by transporting the exact model ledger \(m_7=862/45\) to the arithmetic setting;
2. prove an upper bound on \(m_8\) strong enough to place \((m_7,m_8)\) below the desired target curve.

There is no need to design a new consumption LP after that. The exact Christoffel engine converts the moments automatically.

The finite evaluations \(\{5,2\}=1/8\) and \(C_7=-17/360\) are in the
[joint lattice certificate](spectator_52_exact_lattice.md) and
[pure-cycle flow certificate](pure_cycle_flow_polytopes.md). Their model
ledger exceeds the displayed floor by \(118513/4057200\). The arithmetic
transport remains a separate proof obligation.

The natural first headline target is now

\[
\boxed{80\%}
\]

rather than merely improving the sixth decimal place of \(79.62\%\).

