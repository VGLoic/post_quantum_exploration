/*
 * Our goal is to implement various STARKs. It will not be the official approach as described in the paper.
 * We will start by following another approach, described in this series of articles [1](https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html), [2](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) and [3](https://vitalik.eth.limo/general/2018/07/21/starks_part_3.html).
 *
 * The first stark we implement is in `range_check.rs`: we prove that we know a polynomial of degree less than 1_000_000 such that 0 <= P(x) <= 9.
 *
 * vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * vvvvvvvvvvvv My notes of the article below vvvvvvvvvvvv
 * vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 * Erasure code toy example
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Let's give the toy example on Wikipedia:
 * Alice wants to send her phone number 123456 to Bob with an atrocious protocol:
 * - half of the messages disappear,
 * - messages longer than 5 characters are not allowed,
 * - communication is very pricey.
 *
 * Alice can therefore not send directly her phone number. Additionally it would be costly to send it in pieces, it would be at least two messages for the two pieces and then two aknowledgements from Bob, so at the very least four messages.
 * What she does:
 *  - she splits the phone number in two pieces:  a = 123 and b = 456,
 *  - she defines f(x) = a + (b - a)(x - 1). We have f(1) = a and f(2) = b,
 *  - she derives c = f(3), d = f(4), e = f(5)
 *  - she sends to Bob the five pieces {a, b, c, d, e}. Each piece should be associated to its input x in order that Bob knows how to reconstruct
 * Bob will receive two pieces, let's say d and e, he can then reconstruct the line:
 * f(x) = - d * (x - 5) + e * (x - 4)
 * Bob can then recover f(1) and f(2) from the line.
 *
 * => The generic idea (in `erasure code`) is to derive a polynomial function with the data at hand, and then oversampling it.
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Simple example with range check
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Alice wants to prove to Bob that she knows a polynomial P such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 *
 * A first proof would be to publish the 1_000_000 values. Bob could then recompute the 1_000_000 values and checking that they are correct.
 * => We want something that is faster than 1_000_000 steps.
 * If the verifier is lazy and checks only a subset of the points, it would be doable (and quite OK) by the prover to provide fake proofs.
 *
 * Let us rewrite the problem by introducing a constraint polynomial C(x) such that C(x) = 0 for 0 <= x <= 9. C is considered a public input and can be built here as x(x - 1)(x - 2)...(x - 9).
 *
 * What Alice wants now is to prove that she knows P such that C(P(x)) = 0 for 1 <= x <= 1_000_000.
 * The values between 1 and 1_000_000 are roots and we can use that to rewrite the equation as follows:
 * ```
 * C(P(x)) = D(x) Z(x) for all x
 * Z(x) = (x - 1)(x - 2)...(x - 1_000_000)
 * ```
 *
 * The problem now becomes: Alice wants to prove that she knows the polynomials P and D verifying the above equations.
 *
 * Let us do this in a three step process:
 * - Alice computes the evaluations of P and D at the points x = 1 to x = 1_000_000_000. Alice is generating data over a way bigger range than the initial one.
 *   She then generates two Merkle trees with these evaluations as a way to commit to these evaluations, one for P and one for D. One can be good with a single Merkle tree.
 *   She sends to Bob the two Merkle roots.
 * - Bob will randomly choose some points and asks Alice that she sends the evaluations of P and D with the Merkle proofs,
 * - Bob verifies the Merkle proofs and verifies that the equation holds at the desired point.
 *
 * In terms of *correctness*, if Alice knows a valid P, she will always be able to generate a proof that will be correctly verified, i.e. the checks of Bob will always pass.
 *
 * In terms of *soundness*, i.e. if Alice has an invalid P, what is the probability that she get caught?
 * C is actually a degree 10 polynomials, C(P) is at least a degree 1_000_000 polynomial.
 * Note: The article claims that P is a 1_000_000 degree polynomial, so that C(P) is a 10_000_000 degree polynomial, I don't understand the saying on the degree of P, maybe it's just an assumed fact for the problem? So I stay with 1_000_000 because this is the degree of Z.
 * Two degree N polynomials will agree on at most N points.
 * Therefore, if Alice forges a different polynomials, she will have at most 1_000_000 similar points, it leaves 999_000_000 different points.
 * The probability that Bob picks a valid point is therefore 1/990 ~ 10^(-2). Now the probability for 16 valid points is 10^(-32).
 *
 * Using [Fiat-Shamir heuristic](https://en.wikipedia.org/wiki/Fiat%E2%80%93Shamir_heuristic), one can turn this into non interactive.
 * The idea is to use the root of the Merkle tree (or two roots, as you prefer), in order to derive which 16 values the verifier "will want to check". Alice can then send in the proof:
 * - the Merkle root (or two),
 * - the evaluation and the proof of the derived points.
 * Bob can then directly verifies all of it.
 *
 * Forging a bad proof would be quite intensive as it would require Alice to generate polynomials with their trees until she finds valid branches. Since the soundness is 1 - 1^(-32), it is considered very safe.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The millionth Fibonnaci number
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * We want to prove that we know the millionth Fibonacci number. In order to prove this, we will prove that we know a polynomial representing the Fibonacci suite.
 *
 * Let us go in a similar route than in the previous case.
 * Let us introduce `C(x_1, x_2, x_3) = x_1 - x_2 - x_3` and `P(x)` the polynomial representing the x-th Fibonacci number.
 * In this case we have `C(P(x + 2), P(x + 1), P(x)) = 0` for all integer x.
 *
 * The problem is then translated to proving that we know P and D such that `C(P(x + 2), P(x + 1), P(x)) = D(x) * Z(x)`, where Z is built using a range of integer.
 *
 * The prover will evaluate P and D over a set of values, overlapping the 1_000_000 target, let's say 1_000_000_000. It will then commit to these evaluations in a Merkle tree.
 * The root of the tree will be used in order to derive which 16 values the prover will give in the proof.
 * For each value, the prover will give the branches of P(v), P(v + 1) and P(v + 2).
 * Additionally, the prover will give the branches for P(0) = 0 and P(1) = 0 as these must respect the Fibonacci suite initial values.
 *
 * In practice, this scheme is not very practical as the numbers in the Fibonacci suite grows fast.
 * => Proving this over a finite field would be a good idea.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * A very important remark about previous soundness check
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Let us assume that the solution must be a polynomial of degree with a degree 1_000_000 or less. Said otherwise, the degree of the solution is known.
 * The assumption of soundness is based on the fact that Alice wants to forge an invalid polynomial P of degree 1_000_000 (or whatever it was first).
 * But what if some values are indeed from a polynomial and some values are not? I.e. the whole set of values are not on the same 1_000_000 degree polynomial.
 *
 * For instance, for every x, Alice generates a random value `p` that she will use to commit for `P(x)` and she will derive `d = C(p) / Z(x)`. The `d` will be used as a commitment value for the `D(x)` value.
 * With these computations, Alice can replace some values and commitments with these freshly computed data.
 * When verifying, Bob will not be able to tell that it was wrong.
 *
 * The protection against this kind of attack is actually quite at the center of the machinery of the STARKs.
 * The engineered tools provide so called `proof of proximity` such that most of the points for P and D match the right kind of polynomial.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Trying to theoritize a bit the problem
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * The problem we care about is:
 * *How do we make sure that given a set of data points, they all belong to the same low-degree polynomial.*
 * Low is to be compared with the number of points, e.g. degree 1_000_000 for 1_000_000_000 points.
 * This problem is called `low-degree testing`.
 *
 * In low-degree testing, the issue goes like this:
 * *Let us consider a set of points. It is claimed that it belongs to a polynomial of degree less than D. How to create the associated succint probabilistic proof that it is true?*
 *
 * If you have the ability to look at every points, a simple solution can be made:
 * - use D points to create the associated degree D-1 polynomial using Lagrange interpolation,
 * - check that the remaining points are on the constructed polynomials.
 * => In this case, we say that we have an infinite number of queries (or at least enough to cover all the cases).
 *
 * Let's try now to decrease the number of queries, while controlling how wrong we can be.
 *
 * *With D points (or less)*, it is not enough as one will always be able to use Lagrange interpolation to build a less than D (at most D-1) degree polynomial with all the points on it.
 *
 * *Let's look at D+1 points*. We can use the first D points to build a <D degree polynomial and then use the last point to check if it belongs to the constructed polynomial.
 * This is not a very strong test as there is a non negligible probability that the remaining points do not belong to this low degree polynomial.
 *
 * *We say that we look at D+k points, or that we do D+k queries*, if a proportion `p` of the points do not belong to the same polynomial, the probability of successfully passing the test is `P_k = (1 - p)^k`.
 * Example with 10% of bad points, p = 0.1:
 * - k = 1, P_1 = 0.9,
 * - k = 5, P_5 = 0.59,
 * - k = 10, P_10 = 0.348678
 * - k = 100, P_100 = 2.656e-5
 *
 *
 * So this is good but unfortunately for us, we care about cases where D can be high. So this protocol can be costly and we would like some proximity proof with less than D queries.
 * => The above direct approach absolutely does not work.
 * => Indirect approaches using auxiliary data have succeed in doing so, they are called `probabilistically checkable proofs of proximity (PCPP)`.
 *
 * We will take a look at the FRI protocol. FRI stands for `Fast RS IOPP` with:
 * - RS: [Reed-Solomon error correction](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction#Constructions),
 * - IOPP: Interactive Oracle Proofs of Proximity.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * A simple example for sublinearity
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Let us consider N, the number of points we have.
 * Let us consider f(x) the <D degree polynomial function where the points lie.
 * We'll consider in this example N = 1_000_000_000 and D = 1_000_000.
 *
 * Let us first define the bivariate polynomial function g(x, y) such that g(x, x^1_000) = f(x).
 * One can derive this function by rewriting f(x) = a_i * x^i = a_i * x^[i % 1_000] * (x^1_000)^[floor(i / 1_000)].
 * E.g. for i = 2345 and a_2345 = 876, we would have 876 * x^2345 = 876 * x^345 * (x^1_000)^2 => 876 * x^345 * y^2.
 * Note that g has at most degree 1_000 in x and 1_000 in y.
 *
 * The prover now commits to the values of g on the square {x, 1 <= x <= N} * {y, 1 <= y <= N} (equivalent to {x, 1 <= x <= N} * {x^1_000, 1 <= x <= N}) by creating a Merkle tree.
 * On the square:
 *  - the diagonal are the values of g(x, x^1_000) = f(x),
 *  - a row contains values derived with y fixed in g(x, y) so it is a polynomial of degree <1000,
 *  - a column contains values derived with x fixed in g(x, y) so it is a polynomial of degree <1000,
 *
 * The verifier will then ask a bunch of rows and columns to the prover, this process can be made non interactive using Fiat-Shamir.
 * For each row or column, the verifier will ask for >1000, e.g. 1010, points, this set MUST contain the point on the diagonal.
 *
 * The prover will give all the evaluations with the merkle proofs.
 *
 * The verifier is able to
 *  - verify the Merkle proofs, i.e. the evaluations are consistent with this commitment,
 *  - verify for each row or column that the values lie on a <1000 degree polynomial.
 *
 * With this the verifier can be convinced that:
 *  - most rows contain values that fit in <1000 degree polynomial,
 *  - most columns contain values that fit in <1000 degree polynomial,
 *  - most of the diagonal points belong to these <1000 degree polynomials.
 *
 * As a consequence, the verifier is convinced that the diagonal points correspond to a degree <1_000_000 polynomial.
 *
 *
 * Some numbers about the protocol:
 * If we pick 30 rows and columns, we get 60 * 1_010 = 60_600 points, this respects the sublinearity requirement as it is less than D = 1_000_000. It is still a lot.
 * In terms of computations it's quite costly as:
 *  - on the verifier side: we need to check the Merkle proofs (log) and interpolate 60 polynomials of degree 1_000 (interpolation is subquadratic with FFT I think) so we get sublinear complexity in D,
 *  - on the prover side: we need to compute g on the square, it is 10^18 values to compute, each evaluation is superlinear so complexity is quite high. Then we need to build the tree but this is only linear.
 *
 * In conclusion, complexity for the prover is quite crazy with this protocol. Let's try to decrease it a bit.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * A note about modular arithmetic
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Let us remember the Fermat's little theorem:
 * ```
 * a^p = a (mod p) for any a in Z_p with p prime.
 * ```
 *
 * Here is an important application:
 * If p - 1 is a multiple of k, then a -> a^k has only (p - 1)/k + 1 values.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Improving our low-degree testing
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * The goal is to improve the past protocol. For instance, the 10^18 evaluations was quite a requirement. Additionally, the actual value of evaluation could go very high.
 *
 * Let us use modular arithmetic on our protocol with a carefully chosen prime of `p = 1_000_005_001`.
 * We have p > N, it was required as we still needed our 1_000_000_000 points.
 * Additionally, p - 1 = 1_000_005_000 = 1_000 * 1_000_005 and we know that x^1_000 has a small image, we can compute the size: [(p-1)/1_000 + 1] = 1_000_006 possible values.
 * The square on which we evaluate g(x, y) is then shrinked from 10^18 to 1_000_000_000 x 1_000_006 = 10^15 values, we gain a factor 1_000 yeay!
 *
 * Now let us do an additional tricky trick observation:
 * We have the original data set of N = 1_000_000_000 points.
 * Let us consider two points with the same image under x -> x^1_000, x_a and x_b, i.e. x_a^1000 = x_b^1000 = m. We have f(a) = g(a, a^1_000) = g(a, m) and f(b) = g(b, b^1_000) = g(b, m).
 * As a consequence f(a) and f(b) are on the same row in the square. *There are actually 1_000 points on each rows from the original data set.*
 *
 * It simplifies the work done by the prover: it will only evaluate and commit one column, i.e. 1_000_0006 values.
 *
 * The verifier will then take a bunch of rows, derive from the original data the 1_000 points that describe them. It can then interpolate to obtain the <1000 degree polynomial of each row.
 * The verifier can then check that the associated value on the column match the interpolation of a row.
 * The verifier also needs to check that the column properly correspond to a <1000 degree polynomial.
 *
 * The article specifies that the prover complexity is now 10^9, I got 10^6 from the single column evaluation, we would need to add the complexity of the evaluation but I don't think it takes it into account in the article. maybe yes actually?
 * Let's think about that.
 *
 * On the verifier side, there is less interpolation because it receives only one column. However it has to do some work on the original data set.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Improving the verifier complexity
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Instead of using the previous g(x, y) such that g(x, x^1000) = f(x), we're gonna use a way smaller degree bound, i.e. 2: g(x, x^2) = f(x).
 * We can use a similar construction than before for the coefficients: 876 * x^2345 = 876 * x^1 * (x^2)^1172 = 876xy^1172.
 * As a consequence:
 *  - the rows are at most degree 2 polynomials,
 *  - the columns are at most N/2 degree polynomials.
 *
 * The previous row check involves now 3 points:
 *  - 2 from the diagonal / original data set to interpolate the line,
 *  - 1 from the committed column.
 *
 * Instead of directly interpolating and check that the column is a <N/2 degree polynomial, we repeat the process:
 * - the prover will take the column as new diagonal, it will generate the commitments. The rows will be at most degree 2 polynomials, the columns will be at most degree N/4 polynomials,
 * - the verifier will be able to do a row check and check that the column is a <N/4 degree polynomial, the recursion continue with N/8,
 * - etc...
 * - until we reach a value of N/m small enough that the verifier can *directly* check the polynomial, i.e. by asking all the values, or at least enough to be safe.
 *
 * By doing so, complexity for the verifier is decreased to logarithmic in N (times log(N) for the Merkle proofs, hence log^2(N)).
 *
 *
 * Note: the real FRI is actually a bit more complex, the finite field is a binary Galois field, exponent for the row is 4 and not 2, and possibly other modifications.
 */

pub mod polynomial;
pub mod range_check;
