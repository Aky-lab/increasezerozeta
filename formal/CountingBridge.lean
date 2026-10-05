import Std

/-!
Finite certificate and count accounting for the spectral counting bridge.
The matrix inertia, zero-block decomposition, tail estimate and actual
arithmetic trace cap are inputs outside this formalization.
-/

namespace ZeroZeta

open Lean.Grind Std

/-- Count accounting after an inertia lower bound has been supplied. -/
theorem simple_accounting (s₁ s₂ pairs positive total dim bad : Nat)
    (hinertia : positive ≤ s₁ + s₂ + pairs)
    (hmultiplicity : s₁ + 2*s₂ + 2*pairs ≤ total)
    (hspectral : dim ≤ positive + bad) :
    2*dim ≤ s₁ + total + 2*bad := by omega

/-- Removing an enlarged-interval collar costs two multiplicity units. -/
theorem simple_collar_accounting
    (s₁ s₂ pairs positive dim bad target collar simple : Nat)
    (hinertia : positive ≤ s₁ + s₂ + pairs)
    (hmultiplicity : s₁ + 2*s₂ + 2*pairs ≤ target + collar)
    (hspectral : dim ≤ positive + bad)
    (hremove : s₁ ≤ simple + collar) :
    2*dim ≤ simple + target + 2*bad + 2*collar := by
  have := simple_accounting s₁ s₂ pairs positive (target+collar) dim bad
    hinertia hmultiplicity hspectral
  omega

/-- Distinct blocks need a single collar charge. -/
theorem distinct_collar_accounting
    (s₁ s₂ pairs positive dim bad collar distinct : Nat)
    (hinertia : positive ≤ s₁ + s₂ + pairs)
    (hspectral : dim ≤ positive + bad)
    (hremove : s₁ + s₂ + 2*pairs ≤ distinct + collar) :
    dim ≤ distinct + bad + collar := by omega

/-- The 80% conversion, with finite dimension and collar corrections. -/
theorem eighty_percent_accounting
    (dim bad target collar simple deficit : Nat)
    (hcount : 2*dim ≤ simple + target + 2*bad + 2*collar)
    (hcap : 10*bad ≤ dim)
    (hdeficit : target ≤ dim + deficit) :
    4*target ≤ 5*simple + 10*collar + 9*deficit := by omega

/-- The same cap gives a 90% distinct-count conversion. -/
theorem ninety_percent_distinct_accounting
    (dim bad target collar distinct deficit : Nat)
    (hcount : dim ≤ distinct + bad + collar)
    (hcap : 10*bad ≤ dim)
    (hdeficit : target ≤ dim + deficit) :
    9*target ≤ 10*distinct + 10*collar + 9*deficit := by omega

section OrderedField

variable {R : Type} [Field R] [LE R] [LT R]
  [LawfulOrderLT R] [IsLinearOrder R] [OrderedRing R]
attribute [local instance] Semiring.natCast Ring.intCast

/-- The cubic certificate with its common denominator cleared. -/
def q3Numerator (x : R) : R := 2519-8232*x+7368*x^2-1932*x^3

/-- Scaled denominator of the bounded rational certificate. -/
def certificateDenominator (α x : R) : R := 2519^2*(1+α*x^2)^3

def boundedCertificate (α x : R) : R :=
  (q3Numerator x)^2 / certificateDenominator α x

theorem certificate_denominator_positive (α x : R) (hα : 0 ≤ α) :
    0 < certificateDenominator α x := by
  have hsq := OrderedRing.sq_nonneg (a := x)
  have hprod := OrderedRing.mul_nonneg hα hsq
  have hbase : 0 < (1 : R)+α*x^2 := by grind
  have hbase2 := OrderedRing.mul_pos hbase hbase
  have hbase3 := OrderedRing.mul_pos hbase2 hbase
  unfold certificateDenominator
  grind

theorem bounded_certificate_nonnegative (α x : R) (hα : 0 ≤ α) :
    0 ≤ boundedCertificate α x := by
  have hd := certificate_denominator_positive α x hα
  have hs := OrderedRing.sq_nonneg (a := q3Numerator x)
  have hi := (Field.IsOrdered.inv_nonneg_iff (a := certificateDenominator α x)).mpr
    (Preorder.le_of_lt hd)
  have hp := OrderedRing.mul_nonneg hs hi
  unfold boundedCertificate
  simpa only [Field.div_eq_mul_inv] using hp

/-- Negative-side domination for an explicit admissible parameter family. -/
theorem negative_certificate_numerator (α t : R)
    (hα : 0 ≤ α) (hαone : α ≤ 1)
    (hαcube : 2519^2*α^3 ≤ (1932 : R)^2) (ht : 0 ≤ t) :
    certificateDenominator α (-t) ≤ (q3Numerator (-t))^2 := by
  have ht2 := OrderedRing.mul_nonneg ht ht
  have ht3 := OrderedRing.mul_nonneg ht2 ht
  have ht4 := OrderedRing.mul_nonneg ht3 ht
  have ht5 := OrderedRing.mul_nonneg ht4 ht
  have ht6 := OrderedRing.mul_nonneg ht5 ht
  have hα2 := OrderedRing.mul_le_mul_of_nonneg_left hαone hα
  have hco2 : (0 : R) ≤ 104885808-19036083*α := by grind
  have hco4 : (0 : R) ≤ 86095872-19036083*α^2 := by grind
  have h1 := OrderedRing.mul_nonneg (OrderedRing.ofNat_nonneg (R := R) 41472816) ht
  have h2 := OrderedRing.mul_nonneg hco2 ht2
  have h3 := OrderedRing.mul_nonneg (OrderedRing.ofNat_nonneg (R := R) 131040168) ht3
  have h4 := OrderedRing.mul_nonneg hco4 ht4
  have h5 := OrderedRing.mul_nonneg (OrderedRing.ofNat_nonneg (R := R) 28469952) ht5
  have hco6 : (0 : R) ≤ 1932^2-2519^2*α^3 := by grind
  have h6 := OrderedRing.mul_nonneg hco6 ht6
  have hsum : (0 : R) ≤ 41472816*t+(104885808-19036083*α)*(t*t)
      +131040168*((t*t)*t)+(86095872-19036083*α^2)*(((t*t)*t)*t)
      +28469952*((((t*t)*t)*t)*t)
      +(1932^2-2519^2*α^3)*(((((t*t)*t)*t)*t)*t) := by grind
  have hexpand : (q3Numerator (-t))^2-certificateDenominator α (-t)
      =41472816*t+(104885808-19036083*α)*(t*t)
        +131040168*((t*t)*t)+(86095872-19036083*α^2)*(((t*t)*t)*t)
        +28469952*((((t*t)*t)*t)*t)
        +(1932^2-2519^2*α^3)*(((((t*t)*t)*t)*t)*t) := by
    unfold q3Numerator certificateDenominator
    grind
  have hdiff : 0 ≤ (q3Numerator (-t))^2-certificateDenominator α (-t) := by
    rw [hexpand]
    exact hsum
  exact OrderedAdd.sub_nonneg_iff.mp hdiff

theorem bounded_certificate_negative_majorant (α x : R)
    (hα : 0 ≤ α) (hαone : α ≤ 1)
    (hαcube : 2519^2*α^3 ≤ (1932 : R)^2) (hx : x ≤ 0) :
    1 ≤ boundedCertificate α x := by
  have hneg : 0 ≤ -x := by grind
  have hn := negative_certificate_numerator α (-x) hα hαone hαcube hneg
  have hxx : - -x = x := by grind
  rw [hxx] at hn
  have hd := certificate_denominator_positive α x hα
  have hdiv := (Field.IsOrdered.le_mul_inv_iff_mul_le (1 : R)
    ((q3Numerator x)^2) hd).mpr (by grind)
  unfold boundedCertificate
  simpa only [Field.div_eq_mul_inv] using hdiv

/-- Indicator domination for an arbitrary finite list and nonnegative weights. -/
theorem count_le_certificate_sum {A : Type} (xs : List A)
    (bad : A → Bool) (f : A → R)
    (hnonnegative : ∀ x ∈ xs, 0 ≤ f x)
    (hmajorant : ∀ x ∈ xs, bad x = true → 1 ≤ f x) :
    (xs.countP bad : R) ≤ (xs.map f).sum := by
  induction xs with
  | nil => simp only [List.countP_nil, List.map_nil, List.sum_nil]; grind
  | cons x xs ih =>
    have htail : ∀ y ∈ xs, 0 ≤ f y := by grind
    have htailbad : ∀ y ∈ xs, bad y = true → 1 ≤ f y := by grind
    have hh := ih htail htailbad
    by_cases hx : bad x = true
    · have hhead := hmajorant x (by simp) hx
      simp only [List.countP_cons, hx, ↓reduceIte, List.map_cons, List.sum_cons]
      grind
    · have hhead := hnonnegative x (by simp)
      simp only [List.countP_cons, hx, List.map_cons, List.sum_cons]
      grind

/-- A trace cap of one tenth implies an integer count cap of one tenth. -/
theorem tenth_count_cap {A : Type} (xs : List A)
    (bad : A → Bool) (f : A → R)
    (hnonnegative : ∀ x ∈ xs, 0 ≤ f x)
    (hmajorant : ∀ x ∈ xs, bad x = true → 1 ≤ f x)
    (hcap : 10*(xs.map f).sum ≤ (xs.length : R)) :
    10*xs.countP bad ≤ xs.length := by
  have h := count_le_certificate_sum xs bad f hnonnegative hmajorant
  have hc : ((10*xs.countP bad : Nat) : R) ≤ (xs.length : R) := by grind
  exact OrderedRing.le_of_natCast_le_natCast _ _ hc

variable [DecidableLE R]

/-- The bad count includes equality at the spectral threshold. -/
theorem threshold_count_cap (xs : List R) (α ε : R)
    (hα : 0 ≤ α) (hαone : α ≤ 1)
    (hαcube : 2519^2*α^3 ≤ (1932 : R)^2)
    (hcap : 10*(xs.map (fun x => boundedCertificate α (x-ε))).sum
      ≤ (xs.length : R)) :
    10*xs.countP (fun x => decide (x ≤ ε)) ≤ xs.length := by
  apply tenth_count_cap xs (fun x => decide (x ≤ ε))
    (fun x => boundedCertificate α (x-ε))
  · intro x _
    exact bounded_certificate_nonnegative α (x-ε) hα
  · intro x _ hx
    have hxle : x ≤ ε := of_decide_eq_true hx
    have hshift : x-ε ≤ 0 := by grind
    exact bounded_certificate_negative_majorant α (x-ε) hα hαone hαcube hshift
  · exact hcap

/-- Conditional finite 80% implication, including collar and dimension loss. -/
theorem finite_bounded_certificate_counting
    (xs : List R) (α ε : R)
    (s₁ s₂ pairs positive target collar simple deficit : Nat)
    (hα : 0 ≤ α) (hαone : α ≤ 1)
    (hαcube : 2519^2*α^3 ≤ (1932 : R)^2)
    (hcap : 10*(xs.map (fun x => boundedCertificate α (x-ε))).sum
      ≤ (xs.length : R))
    (hinertia : positive ≤ s₁+s₂+pairs)
    (hmultiplicity : s₁+2*s₂+2*pairs ≤ target+collar)
    (hspectral : xs.length ≤ positive+xs.countP (fun x => decide (x ≤ ε)))
    (hremove : s₁ ≤ simple+collar)
    (hdeficit : target ≤ xs.length+deficit) :
    4*target ≤ 5*simple+10*collar+9*deficit := by
  have hcount := simple_collar_accounting s₁ s₂ pairs positive xs.length
    (xs.countP (fun x => decide (x ≤ ε))) target collar simple
    hinertia hmultiplicity hspectral hremove
  have hbad := threshold_count_cap xs α ε hα hαone hαcube hcap
  exact eighty_percent_accounting xs.length
    (xs.countP (fun x => decide (x ≤ ε))) target collar simple deficit
    hcount hbad hdeficit

end OrderedField

theorem cubic_model_cap : (247 : Rat)/2519 < 1/10 := by decide +kernel
theorem cubic_model_margin : (1 : Rat)/10-247/2519=49/25190 := by decide +kernel
theorem cubic_counting_conversion : (1 : Rat)-2*(247/2519)=2025/2519 := by decide +kernel
theorem cubic_model_trace_identity :
    ((6345361 : Rat)-41472816+104885808*(4/3)-131040168*2
      +86095872*(13/4)-28469952*(101/18)+3732624*(640/63))/6345361
      =247/2519 := by decide +kernel
theorem rational_parameter_admissible :
    (0 : Rat) ≤ 4/5 ∧ (4 : Rat)/5 ≤ 1 ∧
      (2519 : Rat)^2*((4 : Rat)/5)^3 ≤ (1932 : Rat)^2 := by decide +kernel

end ZeroZeta

#print axioms ZeroZeta.simple_accounting
#print axioms ZeroZeta.simple_collar_accounting
#print axioms ZeroZeta.distinct_collar_accounting
#print axioms ZeroZeta.eighty_percent_accounting
#print axioms ZeroZeta.ninety_percent_distinct_accounting
#print axioms ZeroZeta.certificate_denominator_positive
#print axioms ZeroZeta.bounded_certificate_nonnegative
#print axioms ZeroZeta.negative_certificate_numerator
#print axioms ZeroZeta.bounded_certificate_negative_majorant
#print axioms ZeroZeta.count_le_certificate_sum
#print axioms ZeroZeta.tenth_count_cap
#print axioms ZeroZeta.threshold_count_cap
#print axioms ZeroZeta.finite_bounded_certificate_counting
#print axioms ZeroZeta.cubic_model_cap
#print axioms ZeroZeta.cubic_model_margin
#print axioms ZeroZeta.cubic_counting_conversion
#print axioms ZeroZeta.cubic_model_trace_identity
#print axioms ZeroZeta.rational_parameter_admissible
