# Stark exploration

Our goal is to document and implement the generation and verification of STARK proof for various problems.

We will not follow the official approach as described in the academic papers, we are actually following the three articles ([1](https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html), [2](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) and [3](https://vitalik.eth.limo/general/2018/07/21/starks_part_3.html)) of Vitalik Buterin on the topic.

In this repository, one will find:
- the notes on the articles in the related [document](./ARTICLES_NOTES.md),
- an implementation of a stark proof generation and verification for a range check problem.

> [!WARNING] 
> The code represents my understanding of the Stark, it is not optimized and is not meant for any serious use.

## Improvements to add:
- use fast Fourier transform instead of Lagrange interpolation,
- use [Montgomery batch inversions](https://books.google.fr/books?id=kGu4lTznRdgC&pg=PA54&lpg=PA54&dq=montgomery+batch+inversion&source=bl&ots=tPJcPPOrCe&sig=Z3p_6YYwYloRU-f1K-nnv2D8lGw&hl=en&sa=X&redir_esc=y#v=onepage&q=montgomery%20batch%20inversion&f=false),
- support modular arithmetic for larger prime fields,
- add measurements on concrete cases,
- add another computation than range check.
