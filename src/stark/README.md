# Stark exploration

Our goal is to document and implement the generation and verification of STARK proof for various problems.

We will not follow the official approach as described in the academic papers, we are actually following the three articles ([1](https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html), [2](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) and [3](https://vitalik.eth.limo/general/2018/07/21/starks_part_3.html)) of Vitalik Buterin on the topic.

In this repository, one will find:
- the notes on the articles in the related [document](./ARTICLES_NOTES.md),
- an implementation of a stark proof generation and verification for a range check problem.

> [!WARNING] 
> The code represents my understanding of the Stark, it is not optimized and is not meant for any serious use.

## Scripts

Scripts are available at the root level in order to work with stark:
- `fri_friendly`: check among a list of primes, the ones that are FRI friendly,
- `range_check`: run the proof generation and verification for the range check problem.


## Profiling

[Samply](https://github.com/mstange/samply) is used as a profiling tool. The recommended way to use it is to first build the executables using the `profiling` profile and then use Samply to run it.

For instance for the `range_check` executable:
```bash
cargo build --profile profiling
samply record ./target/profiling/range_check
```

## Improvements to add

### Batch inversion
The current implementation uses the naive approach for computing inverses in a finite field, which is quite costly. We can improve it by using batch inversion. The [Montgomery batch inversion method](https://books.google.fr/books?id=kGu4lTznRdgC&pg=PA54&lpg=PA54&dq=montgomery+batch+inversion&source=bl&ots=tPJcPPOrCe&sig=Z3p_6YYwYloRU-f1K-nnv2D8lGw&hl=en&sa=X&redir_esc=y#v=onepage&q=montgomery%20batch%20inversion&f=false) is a good candidate for this.

### Improve arithmetic operations

Currently, arithmetic operations in the finite field are implemented in a naive way. In particular, reduction modulo a prime is done in every operation, addition or multiplication. There are definitely room for improvement on the various computations, especially the ones involving multiple additions and multiplications.

### Use Goldilocks arithmetic

The Goldilocks prime (2^64 - 2^32 + 1) allows for very efficient arithmetic operations using 64-bit integers. Implementing arithmetic operations using the Goldilocks prime can significantly speed up computations in the finite field.

### Use parallelism

Many parts of the implementation can benefit from parallelism. For instance, fast fourier transforms and Merkle tree computations can be parallelized to take advantage of multi-core processors. It would significantly speed up the proof generation and verification processes.

### Add measurements with various metrics

There are no measurements yet on the implementation. It would be interesting to add measurements on various metrics such as:
- proof generation time,
- proof verification time,
- proof size,
- memory usage,
- CPU usage.


### Add another computation than range check

The current implementation only covers the range check problem. It would be interesting to implement other computations.
