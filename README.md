<h1 align="center">ark-msm</h1>

`ark-msm` is a Multi-Scalar Multiplication (MSM) implementation that
incorporates the state-of-the-art MSM optimizations into
[arkworks](https://github.com/arkworks-rs/) with software engineering best
practices in mind. Theories and implementations of the optimization techniques
are documented for further optimizations and development.

**WARNING:** Current implemenation has not received careful code review and is
NOT ready for production use.

## Usage:

* Install [Rust](https://www.rust-lang.org/tools/install):
    ```bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env
    ```

* Run Tests:
    ```bash
    cargo test
    ```

* Benchmarking:
    ```bash
    cargo bench
    ```

## Acknowledgement

The initial implementation of ark-msm was funded by a grant from the [MINA
Foundation](https://minaprotocol.com).

Our algorithms and implementations were heavily based on 
[yrrid software's 2022 zprize submission](https://github.com/yrrid/submission-wasm-msm)
 written in c.

We would like to acknowledge Gregor Mitscha-Baude from O(1) Labs and Niall
Emmart from Yrrid Software for generously taking the time to answer our
technical questions.
