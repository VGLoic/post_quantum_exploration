# Notes on the articles

These notes are based on the series of articles of Vitalik Buterin about STARK:
- [STARKs, Part I: Proofs with Polynomials](https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html),
- [STARKs, Part II: Thank Goodness It's FRI-day](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html),
- [STARKs, Part 3: Into the Weeds](https://vitalik.eth.limo/general/2018/07/21/starks_part_3.html).

The notes can be considered as an addition to the articles, it is a way for me to spend time on the parts that I find hard or that I don't know.

If you are here, please, read the articles first.

## What is our goal?

As a disclaimer, I will mostly work in the context of `stark` and not `zero-knowledge stark`, meaning I will not care about privacy **at first** but I'll still try to keep an eye on it.

Let us consider that we have some public computation function `f`. A prover wants to generate a proof that it knows `x` such that `f(x) = y`, a verifier wants to take such proof and check whether or not the prover is cheating with high probability.

The process is additionally described by the following properties:
- if `zero-knowledge` is kept, the verifier must not learn anything about `x` during the verification,
- the proof must be **succint**, the verifier must be able to run the verification in a short time compared to the original computation `f`,
- the process must be **complete**: if the prover knows a valid `x` for a `y`, it must be able to generate a proof that will pass the verification of the verifier,
- the process must be **sound**: if the prover cheats, the verifier should reject the proof with high probability.


As a spoiler, most of the hard work that will come later is around *succintness* and *soundness*.

## A disgression on erasure coding 

The mathematics and the idea that will be used is close to the domain of [erasure coding](https://en.wikipedia.org/wiki/Erasure_code):
- given a set of data, one can derive an associated polynomial function and use it in order to oversample the original set,
- the oversampled set of data is sent to the recipient,
- if data is corrupted, for any reason, the recipient is able to detect it. Recovery of valid data is also a main feature but not for interest in our case.


### Down to Earth example

Let us give first the toy example on the Wikipedia page:
Alice wants to send her phone number 123456 to Bob with an atrocious protocol:
- half of the messages disappear,
- messages longer than 5 characters are not allowed,
- communication is very pricey.

Alice can therefore not send directly her phone number. Additionally it would be costly to send it in pieces, it would be at least two messages for the two pieces and then two aknowledgements from Bob, so at the very least four messages.
What she does:
 - she splits the phone number in two pieces:  a = 123 and b = 456,
 - she defines f(x) = a + (b - a)(x - 1). We have f(1) = a and f(2) = b,
 - she derives c = f(3), d = f(4), e = f(5)
 - she sends to Bob the five pieces {a, b, c, d, e}. Each piece should be associated to its input x in order that Bob knows how to reconstruct
 - Bob will receive two pieces, let's say d and e, he can then reconstruct the line:
 - f(x) = - d * (x - 5) + e * (x - 4)
 - Bob can then recover f(1) and f(2) from the line.

 If Bob would have received corrupted data, i.e. the five pieces but with data changed. Bob would be able to check that the pieces do not fit on the same line, hence he will say with certainty that the data is not correct.

 ### A more generic approach with polynomials

 Given a data set of N evaluations, one can build a [N - 1 degree polynomial from it](https://en.wikipedia.org/wiki/Lagrange_polynomial).

 Using this polynomial, we can oversample the original data set, say `N + k`.

 Now, a verifier will take the first N evaluations it receives and build the associated degree N - 1 polynomial. It can then verify if the `k` other points belong to the built polynomial.

 Additionally, if the prover changes one original point and keep the other N - 1 evaluations, this new polynomial is still of order N - 1 but the `k` new evaluations of oversampling will be completely different. In some way, it gives a way to amplify differences from small modifications of the original data.

## Considering the range check problem

Alice wants to prove to Bob that she knows a polynomial `P` such that `0 <= P(x) <= 9, x ∈ [1, 1_000_000]`.

We will add a new hypothesis to the article: P is meant to be a polynomial less than degree 1_000_000.

### The simple one

As a first proof, Alice could send all the evaluations of P to Bob. Bob would then check that all evaluations are indeed between 0 and 9, and he would be convinced because he knows he can build a polynomial from these evaluations. This is not very efficient because Bob would need the 1_000_000 steps to check the statement. Additionally, Alice completely reveals her solution, which is sad.

If Bob would do less than 1_000_000 checks, Alice could start to cheat with a pretty good probability.
For instance if Bob does 750_000 checks, chosen in an uniform way, and Alice cheats on 100 value (999_900 valid values), the probability that the proof is verified would be 
```
1 - (999_000/1_000_000)(998_999/999_999)...((999_000 - 750_000 + 1)/(1_000_000 - 750_000 + 1))
```

### A stronger alternative

Let's try another approach and be closer to the original problem statement `f(x) = y` of the beginning. We introduce a constraint polynomial `C` such that `C(x) = 0, x ∈ [1, 1_000_000]`. This polynomial is straightforwardly defined as
```
C(x) = (x - 1)(x - 2)...(x - 1_000_000)
```
This constraint polynomial is considered a public parameter of the problem, so this is known by the verifier and the prover.

Now our original problem becomes 
```
C(P(x)) = 0, x ∈ [1, 1_000_000]
```

Well, this is nice but we have not actually moved much in terms of nicer proof for this than before.

Let's try to use our erasure coding approach here. Instead of restricting the interval up to 1_000_000, we go to 1_000_000_000.
In this case we can rewrite the problem as
```
C(P(x)) = Z(x) * D(x), x ∈ [1, 1_000_000_000]
```
where `Z(x) = (x - 1)(x - 2)...(x - 1_000_000)`. The `D` polynomial may be computed using a division of `C(P)` by `Z`.
The `Z` polynomial is also considered as a public parameter of the system.

Let us consider the new proof system:
- Alice computes the evaluations of `P` and `D` over the range 1 to 1_000_000_000. Since we become a bit more formal, Alice also takes all the evaluations and put it in a Merkle tree as a way to commit to them,
- Bob selects points at which he wants the evaluations and ask Alice to give him,
- Alice gives the asked evaluations, the Merkle proofs and the root hash of the Merkle tree,
- Bob verifies the Merkle proofs and verifies that the equation `C(P(x)) = Z(x) * D(x)` holds at the selected point.

We blew up the size of the space but in exchange, we know that if changes to `P` or `D` happen, it will have huge impact on the system.

For instance, let us consider 1 bad point in P, e.g. `P(1_000_000) = 10`, in addition to this original evaluation, it will change all the oversampled data: the 990_000_000 evaluations. The probability that Bob picks an invalid point is therefore
```
999_000_001 / 1_000_000_000 > 99%
```

Instead of one points, Bob is encouraged to pick multiple points, e.g. 16, in order to obtain a extremely high probability that such a fake proof would be detected.

The process can be made non interactive using the [Fiat-Shamir heuristic](https://en.wikipedia.org/wiki/Fiat%E2%80%93Shamir_heuristic): the selected points by Bob are pseudo randomly selected by Alice based on the root hash of the evaluations commitments. As a consequence, Alice can generate the proof on her side without the need of Bob. When verifying, Bob must check that the correct points have been selected.

### Completing the soudness check

We saw that the previous version was quite strong against modification of the original data, i.e. if something is changed in the original data, the oversampling will differ and the chance of detecting an error will blow up.

Now, what would happen if the prover modifies directly the oversampled data set. It can easily do so in a way that still respect the equation `C(P(x)) = D(x) * Z(x)`: a prover can simply choose a random number `p_cheat`, set `p_cheat` at the evaluation of `P` for a particular point `x_cheat`, and set `C(p_cheat) / Z(x_cheat)` at the evaluation of `D` for the point `x_cheat`. Our test above would actually not catch it because all we do is check that the relation holds.

What we don't do is ensure that the evaluations of `P` correspond to a polynomial of degree less than 1_000_000. This attack is actually hitting there because by modifying the oversampled data, we are replacing a derived evaluation by an additional constraint, hence actually most probably increasing the degree of P.

The heart of the stark is actually focusing on this new aspect: proving that the polynomial behind the proof has a degree less than a fixed degree.
