# Stark exploration

Our goal is to document and implement the generation and verification of STARK proof for various problems.

We will not follow the official approach as described in the academic papers, we are actually following the three articles ([1](https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html), [2](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) and [3](https://vitalik.eth.limo/general/2018/07/21/starks_part_3.html)) of Vitalik Buterin on the topic.

In this repository, one will find:
- the notes on the articles in the related [document](./ARTICLES_NOTES.md),
- an implementation of a stark proof generation and verification for a range check problem.


## Improvements to add:
- use fast Fourier transform instead of Lagrange interpolation,
- support modular arithmetic for larger prime fields,
- add measurements on concrete cases,
- add another system than range check.
