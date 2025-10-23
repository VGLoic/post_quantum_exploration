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

## A digression on erasure coding 

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
- Bob verifies the Merkle proofs and verifies that the equation `C(P(x)) = Z(x) * D(x)` holds at the selected points.

We blew up the size of the space but in exchange, we know that if changes to `P` or `D` happen, it will have huge impact on the system.

For instance, let us consider 1 bad point in P, e.g. `P(1_000_000) = 10`, in addition to this original evaluation, it will change all the oversampled data: the 990_000_000 evaluations. The probability that Bob picks an invalid point is therefore
```
999_000_001 / 1_000_000_000 > 99%
```

Instead of one points, Bob is encouraged to pick multiple points, e.g. 16, in order to obtain a extremely high probability that such a fake proof would be detected.

The process can be made non interactive using the [Fiat-Shamir heuristic](https://en.wikipedia.org/wiki/Fiat%E2%80%93Shamir_heuristic): the selected points by Bob are pseudo randomly selected by Alice based on the root hash of the evaluations commitments. As a consequence, Alice can generate the proof on her side without the need of Bob. When verifying, Bob must check that the correct points have been selected.

### Completing the soudness check

We saw that the previous version was quite strong against modification of the original data, i.e. if something is changed in the original data, the oversampling will differ and the chance of detecting an error will blow up.

Now, what would happen if the prover modifies directly the data set without relying on any polynomials. It can easily do so in a way that still respect the equation `C(P(x)) = D(x) * Z(x)`: a prover can simply choose a random number `p_cheat`, set `p_cheat` at the evaluation of `P` for a particular point `x_cheat`, and set `d_cheat = C(p_cheat) / Z(x_cheat)` at the evaluation of `D` for the point `x_cheat`. Our test above would actually not catch it because all we do is check that the relation holds.

What we don't do is ensure that the evaluations of `P` correspond to a polynomial of degree less than 1_000_000. This attack is actually hitting there because by modifying the data directly, we are replacing a supposedly derived evaluation by an additional constraint, hence actually most probably increasing the degree of P (at least the P that would satisfy all the evaluation).

The heart of the stark is actually focusing on this new aspect: proving that the polynomial behind the proof has a degree less than a fixed degree.

## Low degree testing

Let us focus on our new problem: **proving with high probability that a set of evaluations belong to a polynomial of degree less than `D`**.

This problem is called *low degree test* in the litterature. The term `low` is used because the degree `D` is meant to be very small compared to the original data set size.

For instance, in our previous problem, the degree we care about was less than 1_000_000 while the data set size was 1_000_000_000.

In order to be consistent with the process, the proof must be **succint, complete and sound**

### The direct approach

The verifier has access to the whole set of evaluations, meaning it can ask the prover for it.

If the verifier asks for all the evaluations, it can use `D` points in order to build a polynomial of order less than `D`(`<D`). It can then check that the other evaluations are consistent with the interpolated polynomial. In this case, we say that we have an infinite number of *queries*. It is perfectly complete and sound but it is definitely not succint.

Let's try to decrease a bit the number of queries to see how the system behaves. We need in any case `D` evaluations in order to interpolate the polynomial, so at least `D` queries. Let us define `k` the additional number of queries we make: in total we make `D + k` queries.

If we consider a data set with a proportion `p` of points that do not belong to the same polynomial, the probability the test passes with `D + k` queries is `(1 - p)^k`.

For instance, with `p = 0.1`:
- k = 1, P = 0.9,
- k = 5, P = 0.59,
- k = 10, P = 0.348678,
- k = 100, P = 2.656e-5

It is decreasing correctly but unfortunately, we are still above the `D` queries, hence we would never be able to reach succintness with this approach.

This approach is the `direct approach` and will be re-used later on, but not for the entire system.

As the *direct approach* does not work, we will introduce *indirect approaches*, these protocols are called *probabilistically checkable proofs of proximity (PCPP)*. We will look at the FRI protocol in our case, it stands for *Fast RS IOPP* with *RS* for [Reed-Solomon error correction](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction#Constructions) and *IOPP* for *Interactive Oracle Proofs of Proximity*.

### Introducing sublinearity

Let us introduce `N` the number of points we have. We'll stick to our problem instance with `N = 1_000_000_000` and `D = 1_000_000`.

Let us also introduce our polynomial function `P(x)` defined from 1 to N.

From `P`, one can derive the bivariate polynomial `G(x, y), x ∈ [1, 1_000_000], y = z^1_000, z ∈ [1, 1_000_000]` such that `G(x, x^1_000) = P(x)`.
For instance, `P(x) = x^4006 + 4x^3021 + 18x^2000 + 7x^1765 + x^2 + 98x + 8`, then `G(x, y) = y^4 * x^6 + 4x^21 * y^3 + 18y^2 + 7yx^765 + x^2 + 98x + 8`.

The space on which `G` lives is a square, on this square:
- rows are defined with `y` fixed, they are described by polynomial of degree less than 1_000,
- columns are defined with `x` fixed, they are described by polynomial of degree less than 1_000,
- the diagonal is the set of original data, `G(x, x^1000) = P(x)`.

An interactive low degree proof can be made with this:
- the prover evaluates `G` over the whole square,
- the prover commits to the evaluations using a Merkle tree,
- the verifier selects a number of rows and columns,
- the verifier asks `D + k` evaluations on each rows and columns, the evaluations must contain the diagonal one,
- the prover sends the evaluations with the Merkle proofs and the root hash of the commitments,
- the verifier verifies the Merkle proofs, that each row fit in a <1000 degree polynomial, that each column fit in a <1000 degree polynomial.

This provides a statistical proof that
- most rows are described by <1000 degree polynomials,
- most columns are described by <1000 degree polynomials,
- most points on the diagonal are on these polynomials.

As a consequence, most of the diagonal points are described by a <1_000_000 degree polynomial.

The interactive process can be made non interactive using the Fiat-Shamir heuristic.

> [!NOTE]
> This is far from a bullet proof demonstration, it is more an intuition demonstration at this stage.

The process is definitively *complete*. I would not be able to analyze its soudness but I will happily rely on the shoulders of giants for the moment. Let us look at the succintness though.

Let's say for the exercise, `k = 10` and the number of queries rows and columns are 60.

For the prover, one evaluation is `O(D)`, we make evaluations over the square, therefore the total is `O(DN^2)`. It is a good thing we care about the verifier because the complexity on the prover's side is gigantic.

For the verifier, we do 60 interpolations of degree <1000 polynomials and 60 * 1010 Merkle proof verifications. An interpolation using the fast Fourier transform is `O(nlog(n))` and a Merkle proof verification is `O(log(n))`, so we get something like `~60 O(1000log(1000)) + 60 * 1010log(N^2)`. This looks sublinear to me!

This describes our first indirect proof, now let's try to improve it.

### A note about modular arithmetic

I have not said anything about modular arithmetic until now. For the rest of these notes, I will simply assume it applied. I encourage you to read about it, I personally like these [notes of Ben Lynn](https://crypto.stanford.edu/pbc/notes/numbertheory/book.pdf).

I simply want to pin point two important aspects here.

1. Let us define `g` the generator of our prime field, we define the `units` as the successive power of `g`:
    ```
    g, g^2, g^3, ..., g^(p-2), 1
    ```
2. Let us recall the Fermat's little theorem:
    ```
    a^p = a (mod p) for any a in Z_p with p prime.
    ```
    As an application, if `p - 1` is a multiple of integer `k`, then the function `a -> a^k` has only `(p - 1)/k + 1` values, i.e. it has a small image. In terms of units, applying the function map `p - 1` units to `(p - 1) / k` units: `g^4, g^8, ..., g^(p-5), 1`

### Going from the square to a rectangle

The square and the bivariate polynomial G we defined earlier were using the function `x -> x^1_000`, because I know what will come after, let us change this to `x -> x^4`. As a consequence, the rows are now described by polynomials with degree less than 4, the columns are described by polynomials of degree less than `D/4`. The spirit of the previously described proof system is unchanged apart that rows are easier to verify but columns are harder.

A few points related to implementation:
- for the rest of the notes, we will take care of choosing a prime field with a prime `N` such that `N - 1` is divisible by 4 as much as possible,
- the `N` becomes the upper limit of our evaluations range. In practice, if we want something close to 1_000_000_000, we will find a prime number close to 1_000_000_000 (and divisible by 4 a few times),
- the evaluations from 1 to `N` can be ordered using the units, as successive power of our generators, e.g. the first evaluation is for `x = g`, the second is `x = g^2`, etc...

And going back to our square, the column is formed by `y ∈ { units^4 }`. And this is actually 4 times smaller than the number of units, so the column direction of our square just shrank and we now have a rectangle instead of a square.

The previous proof approach still works, it is just that the prover can skip a bunch of evaluations because our space just got smaller.

In terms of succintness, we did not move on the verifier side, however we simplified the work of the prover as the evaluations are now of order `O(D * N * N/4)`.

### An observation on the diagonal

Since the square became a rectangle, the diagonal has interesting properties. It is actually wrapped, it will reach the "bottom" of the rectangle and then re-start at the top.
This wrapping will actually occur 4 times.
As a consequence, the diagonal evaluations of `P` contains 4 points on each row.

For instance, the row defined by `y = g^4` has the following four points:
- `x = g, x^4 = g^4, P(x) = G(g, g^4)`,
- `x = g^(1 + (N-1)/4), x^4 = g^4 * g^(N-1) = g^4, P(x) = G(g^(1 + (N-1)/4), g^4)`,
- `x = g^(1 + 2 * (N-1)/4), x^4 = g^4 * g^2(N-1) = g^4, P(x) = G(g^(1 + 2 * (N-1)/4), g^4)`,
- `x = g^(1 + 3 * (N-1)/4), x^4 = g^4 * g^3(N-1) = g^4, P(x) = G(g^(1 + 3 * (N-1)/4), g^4)`.

### Limiting evaluations to the bare minimum

Instead of evaluating and committing to the whole rectangle and then sending a bunch of rows and columns to the verifier, we will take a drastically less consuming approach:

The prover will:
- evaluate and commit to the diagonal, i.e. `P(x) where x ∈ { units }`,
- evaluate and commit to one column, i.e. `G(x, y) where x = x_c fixed and y ∈ { units^4 }`,
- build the proof:
    - select rows, for each row, the four matching diagonal points are added and the matching column point is added,
    - select enough evaluations of the column to do a **direct** (as in the paragraph `The direct approach`) degree proof of the `G` polynomial on the column,

The verifier will:
- verify that each row is a degree <4 polynomial,
- verify that the column is a degree `<D/4` polynomial,
- verify all the Merkle proofs.

With this, the verifier is convinced:
- most points on the diagonal form rows described by polynomials of degree <4,
- the column is described by a polynomial of degree `<N/4`.

As a consequence, the verifier is convinced that: most diagonal points lie on a `<D` degree polynomial.

> [!NOTE]
> I have a hard time about this part. I understand the fact that we don't need to select rows anymore as we can sample them from the diagonal evaluations. However, I don't understand the fact that we can suddenly work with only one column.

This is an amazing improvement for the prover as it now needs to do evaluations over the set of units for the diagonal (`O(DN)`) and over the column (`O(DN/4)`).

For the verifier, it is an improvement as only a single column needs to be interpolated.

### Going recursive

We started with a polynomial of a degree `<D` and we have a proof about two polynomials, one of degree <4 and one of degree `<D/4`.
The degree <4 is very fast to verify as it implies only 5 evaluations, the degree `<D/4` may still be long to verify.

Now, we have chosen our space in such a way that it is divisible by 4 a few times at least.
This gives us the opportunity to repeat the above process for the polynomial `P_1(x) = G(x_c, x), x ∈ { units^4 }` with now original degree of `<D/4`.
We can consider the rectangle `(x, y), x ∈ { units^4 }, y ∈ { units^8 }` and define the polynomial `G_1(x, y)` such that `G_1(x, x^4) = P_1(x)`. So the column of the last iteration is the diagonal of the new iteration, and at each iteration, the degree of the column polynomial is divided by 4.
The process stops when the diagonal has a degree low enough to allow a direct proof, i.e. simply providing enough evaluations to interpolate and check the degree.

The prover will now: check if the degree to prove is above a fixed threshold
- if no: the prover selects enough, e.g. more than degree threshold, evaluations with their Merkle roots, it adds these to the proof,
- if yes:
    - the prover will select a column point,
    - it evaluates the polynomial over the column and commit to these evaluations,
    - it selects rows, add the associated diagonal evaluations with their Merkle proofs, it also adds the matching evaluation on the column with its Merkle proof, all of this is added to the proof,
    - the column evaluations is used as starting point for the process again.

In terms of verification:
- for each indirect proof,
    - interpolate the rows based on the diagonal evaluations, check that the rows are described by polynomials of degree `<4` and that it correctly interpolates to the committed column evaluation,
    - verify the Merkle proofs of the diagonal evaluations and the column evaluations,
- interpolate the evaluations of the direct proof and verify that it belongs to a polynomial of degree less than the degree threshold.

This recursive process is nice for the prover but is mostly for the verifier. The complexity of the prover is still capped by the diagonal evaluations at `O(DN)`. However, for the verifier, we don't have to do the costly interpolation of the column polynomial so it is a great win, the complexity goes down to a logarithmic.

## The final proof

Now we have our ingredients to define the final proof.

Before the development on the low degree proof, we stopped with a way to prove with good probability that the equation `C(P(x)) = Z(x) * D(x)` holds, i.e. oversampling the original evaluations with a way bigger range, committing to the data and then selects a few evaluations in order to let the verifier verifies that the equation holds at the selected points. We'll keep that part and call it the `spot checks` part.

The weakness in the above proof was that the prover could cheat by modifying directly the evaluations with arbitrary values, though still verifying the relations. We argued that by cheating in this way, the prover was actually increasing the degree of the polynomials underlying the proof.

The low degree proof has been described as a solution to this problem, so we can now complete the proof to defend against this kind of attack. The prover evaluates two polynomials, the `P` and `D` polynomials, hence we will need two low degree proofs:
- one for `P` with maximum degree `D = 1_000_000`,
- one for `D` with maximum degree `10 * D - D = 9 * D = 9_000_000`. The factor 10 is coming from `C` while the `- D` is coming from the division of `C(P(x))` by `Z(x)`.
