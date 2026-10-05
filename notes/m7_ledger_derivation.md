# Bell(7) derivation of the seventh-moment ledger

The combinatorial reduction assumes the frozen-singleton identities and 3-block vanishing rule used at m5 and m6. The corrected model value {5,2}=1/8 is exact; C7 and arithmetic transport remain unresolved.

## 1. Bell(7) classes

The block-size signatures and multiplicities are

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

with total \(877\).

## 2. Classes containing a 3-block

The lower-moment local law kills the classes whose nontrivial connected structure contains a 3-block:

\[
3\,1^4,\quad
3\,2\,1^2,\quad
3\,2^2,\quad
3^2\,1,\quad
4\,3.
\]

This is the same mechanism already consumed in the \(m_5,m_6\) ledgers.

## 3. Pair/four-cycle layer

The all-singleton, one-pair, two-pair and one-four-block classes give

\[
1+7+70t_{\rm adj}+35t_{\rm opp}+35\Phi_4.
\]

Using

\[
t_{\rm adj}=\frac7{60},\qquad
t_{\rm opp}=\frac1{30},\qquad
\Phi_4=-\frac1{60},
\]

this is exactly

\[
\boxed{\frac{67}{4}}.
\]

## 4. Frozen-singleton lifts

Deleting singleton blocks leaves overlap ranges unchanged, hence

\[
2^3\,1:\quad 7\{2,2,2\},
\]

\[
4\,2\,1:\quad 7\{4,2\},
\]

\[
5\,1^2:\quad 21C_5,
\]

\[
6\,1:\quad 7C_6.
\]

## 5. New seventh-order classes

Only

\[
\{5,2\},\qquad C_7
\]

remain genuinely seventh-order. Here \(\{5,2\}\) already denotes the sum over all \(21\) placements.

Therefore

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

Substituting the exact lower-order constants gives

\[
\boxed{
m_7
=
\frac{1717}{90}
+\{5,2\}
+C_7.
}
\]

## 6. Corrected model input

The [exact lattice calculation](spectator_52_exact_lattice.md) gives {5,2}=1/8. The ledger therefore becomes

    m7 = 6913/360 + C7.

The numerical C7 candidate -17/360 would give m7=862/45. This remains a conditional scenario until C7 and the arithmetic transports are established.

The Stieltjes floor implied by the supplied moments through m6 is

    m7* = 25866469/1352400.

The explicit-input calculator in `scripts/m7_ledger.py` keeps scenario arithmetic separate from moment certification.
