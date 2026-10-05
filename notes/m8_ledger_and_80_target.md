# Bell(8) reduction: the five new classes controlling the 80% target

**Status:** exact combinatorial reduction conditional on the same frozen-singleton and 3-block-vanishing identities used below order eight, plus the currently pre-registered seventh-order inputs.

Let

\[
A_8:=\{2^4\}+\{4,2,2\}+\{4,4\}+\{6,2\}+C_8.
\]

The main point of this note is that, once lower classes are inserted, **all genuinely new eighth-order information enters only through \(A_8\)**.

## 1. Surviving Bell(8) signatures

Among the \(B_8=4140\) set partitions, the signatures not containing a 3-block are

\[
1^8,\;
2\,1^6,\;
2^2\,1^4,\;
2^3\,1^2,\;
2^4,\;
4\,1^4,\;
4\,2\,1^2,\;
4\,2^2,\;
4^2,
\]
\[
5\,1^3,\;
5\,2\,1,\;
6\,1^2,\;
6\,2,\;
7\,1,\;
8.
\]

Their multiplicities are respectively

\[
1,28,210,420,105,70,420,210,35,56,168,28,28,8,1.
\]

All signatures containing a 3-block are discarded by the same local-law mechanism already used in the lower ledgers.

## 2. Pair/four-cycle base layer

The all-singleton class contributes \(1\).

There are \(28\) one-pair placements, each contributing \(1/3\).

For two pairs, choose four endpoints in \(\binom84=70\) ways. Each four-set has two non-crossing pairings and one crossing pairing, giving

\[
140t_{\rm adj}+70t_{\rm opp}.
\]

There are \(70\) four-block placements, each contributing \(\Phi_4\).

Therefore

\[
B_8^{\rm base}
=
1+\frac{28}{3}
+140t_{\rm adj}
+70t_{\rm opp}
+70\Phi_4.
\]

With

\[
t_{\rm adj}=\frac7{60},\qquad
t_{\rm opp}=\frac1{30},\qquad
\Phi_4=-\frac1{60},
\]

this is exactly

\[
\boxed{B_8^{\rm base}=\frac{167}{6}}.
\]

## 3. Frozen-singleton lifts

Deleting singleton blocks preserves overlap ranges.

Thus

\[
2^3\,1^2:\qquad
\frac{420}{15}\{2,2,2\}=28\{2,2,2\},
\]

because the six-cycle aggregate has \(15\) placements.

Similarly,

\[
4\,2\,1^2:\qquad
\frac{420}{15}\{4,2\}=28\{4,2\}.
\]

For the pure connected blocks,

\[
5\,1^3:\quad 56C_5,
\]

\[
6\,1^2:\quad 28C_6,
\]

\[
7\,1:\quad 8C_7.
\]

For the seventh-order joint class,

\[
5\,2\,1:\qquad
\frac{168}{21}\{5,2\}=8\{5,2\}.
\]

## 4. Genuinely new eighth-order classes

The five new objects are exactly

\[
\{2^4\},\qquad
\{4,2,2\},\qquad
\{4,4\},\qquad
\{6,2\},\qquad
C_8.
\]

Hence

\[
\boxed{
m_8
=
\frac{167}{6}
+28\{2,2,2\}
+28\{4,2\}
+56C_5
+8\{5,2\}
+28C_6
+8C_7
+A_8.
}
\]

Using the exact lower-order values

\[
\{2,2,2\}=\frac{131}{420},\quad
\{4,2\}=-\frac{23}{420},\quad
C_5=\frac1{36},\quad
C_6=-\frac1{126},
\]

and the current seventh-order candidates

\[
\{5,2\}=\frac18,\qquad
C_7=-\frac{17}{360},
\]

the inherited part is

\[
\boxed{
m_8=\frac{3329}{90}+A_8.
}
\]

## 5. Exact target interval when \(m_7=862/45\)

If the seventh-moment candidate

\[
m_7=\frac{862}{45}
\]

is eventually certified, the ordinary Hankel moment constraint gives

\[
m_8\ge37.2618657631\ldots,
\]

so

\[
\boxed{
A_8\ge0.2729768742\ldots.
}
\]

The exact degree-four Christoffel target curves give:

### 80%

\[
m_8\le37.4791244743\ldots
\]

or equivalently

\[
\boxed{
A_8\le0.4902355854\ldots.
}
\]

### 80.2%

\[
\boxed{
A_8\le0.4129766383\ldots.
}
\]

### 81%

\[
\boxed{
A_8\le0.3290373911\ldots.
}
\]

### 90%

At this moment order alone the required band would be

\[
\boxed{
A_8\le0.2768820728\ldots,
}
\]

only slightly above the moment-cone floor. This quantifies why 90% will almost certainly require higher moments rather than merely sharpening the eighth-order calculation.

## 6. Research priority

For an 80% theorem, it is not necessary to identify every new eighth-order constant exactly.

It is enough to prove the single aggregate upper bound

\[
\boxed{
A_8\le0.4902355854\ldots.
}
\]

Exact evaluation remains desirable for auditability and future rungs, but a rigorous aggregate upper bound is the shortest route to the headline.
