# Bell(7) derivation of the seventh-moment ledger

**Status:** exact combinatorial reduction, conditional only on extending the same frozen-singleton identities and 3-block vanishing rule already consumed at \(m_5,m_6\).

This note derives the \(m_7\) ledger from the 877 set partitions of a seven-cycle. It does not assume a numerical value for the new \(\{5,2\}\) joint class or for \(C_7\).

## 1. Bell(7) classes

The block-size signatures and multiplicities are:

\[
\begin{array}{c|r}
\text{signature}&\text{count}\\ \hline
1^7&1\\
2\,1^5&21\\
2^2\,1^3&105\\
2^3\,1&105\\
3\,1^4&35\\
3\,2\,1^2&210\\
3\,2^2&105\\
3^2\,1&70\\
4\,1^3&35\\
4\,2\,1&105\\
4\,3&35\\
5\,1^2&21\\
5\,2&21\\
6\,1&7\\
7&1
\end{array}
\]

The total is \(877\).

## 2. Classes containing a 3-block

The lower-moment arithmetic local law kills the classes whose nontrivial connected structure contains a 3-block. Thus the signatures

\[
3\,1^4,\quad
3\,2\,1^2,\quad
3\,2^2,\quad
3^2\,1,\quad
4\,3
\]

do not contribute.

This is the same rule that removes the 3-block classes in the \(m_5,m_6\) ledgers, including the separately audited \(\{3,3\}\) class.

## 3. Pair/four-cycle layer

The surviving classes built only from pairs, a four-block and singletons collapse to the certified four-cycle anchors.

### One pair

There are

\[
\binom72=21
\]

placements, each contributing the two-cycle constant \(1/3\). Hence

\[
21\cdot\frac13=7.
\]

### Two pairs

Choose the four paired vertices in

\[
\binom74=35
\]

ways. For any cyclically ordered four-set, its three pairings consist of two non-crossing pairings and one crossing pairing. Therefore

\[
70\,t_{\rm adj}+35\,t_{\rm opp}.
\]

### One four-block

There are

\[
\binom74=35
\]

placements, each reducing by frozen singletons to \(\Phi_4\).

Together with the all-singleton class, the pair/four layer is

\[
1+7+70t_{\rm adj}+35t_{\rm opp}+35\Phi_4.
\]

Using

\[
t_{\rm adj}=\frac7{60},\qquad
t_{\rm opp}=\frac1{30},\qquad
\Phi_4=-\frac1{60},
\]

this collapses exactly to

\[
\boxed{\frac{67}{4}}.
\]

## 4. Frozen-singleton lifts of lower connected classes

A singleton block repeats one walk value and does not change the overlap range. Thus deleting singleton blocks maps the following seven-cycle classes to already-defined lower-cycle aggregates.

### Three pairs plus one singleton

There are \(105\) partitions of type \(2^3\,1\). The six-cycle aggregate \(\{2,2,2\}\) contains \(15\) placements. Hence the seven-cycle class is

\[
\frac{105}{15}\{2,2,2\}
=
7\{2,2,2\}.
\]

### Four-block, pair, singleton

There are \(105\) partitions of type \(4\,2\,1\). The six-cycle \(\{4,2\}\) aggregate contains \(15\) placements. Hence

\[
7\{4,2\}.
\]

### Five-block plus two singletons

There are

\[
\binom75=21
\]

placements, each reducing to the pure five-cycle constant \(C_5\). Hence

\[
21C_5.
\]

### Six-block plus singleton

There are \(7\) placements, yielding

\[
7C_6.
\]

## 5. New seventh-order classes

Two genuinely seventh-order objects remain:

\[
\{5,2\},
\qquad
C_7.
\]

Here \(\{5,2\}\) already denotes the sum over all \(21\) placements of the five-block/pair joint class; no additional multiplicity is applied.

## 6. Exact ledger identity

Consequently,

\[
\boxed{
m_7
=
\frac{67}{4}
+21C_5
+7\{2,2,2\}
+7\{4,2\}
+\{5,2\}
+7C_6
+C_7.
}
\]

Substituting the lower-order exact constants

\[
C_5=\frac1{36},\qquad
\{2,2,2\}=\frac{131}{420},\qquad
\{4,2\}=-\frac{23}{420},\qquad
C_6=-\frac1{126},
\]

gives the affine form

\[
\boxed{
m_7
=
\frac{17139}{900}
+\{5,2\}+C_7
}
\]

or, more simply after combining the rational part,

\[
\boxed{
m_7
=
\frac{3431}{180}
+\{5,2\}+C_7.
}
\]

With the currently pre-registered/identified candidates

\[
\{5,2\}=\frac18,\qquad
C_7=-\frac{17}{360},
\]

this becomes

\[
\boxed{
m_7=\frac{862}{45}=19.155555\ldots.
}
\]

The last numerical value remains a **candidate** until both new seventh-order constants and their arithmetic transports are certified.
