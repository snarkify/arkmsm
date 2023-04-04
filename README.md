<h1 align="center">ArkMSM</h1>

`arkmsm` is a Multi-Scalar Multiplication (MSM) module that incorporates
state-of-the-art MSM optimizations into
[arkworks](https://github.com/arkworks-rs/). This implementation is extensively
documented, enabling developers to quickly learn and experiment with lastest
MSM optimization techniques. Additionally, the code is designed to be modular,
facilitating easy integration with other libraries and projects.

Please note that **the current implementation is intended for 
research and study purposes only and has not been audited for production use. 
Use at your own discretion.**

## Up to 2x Performance Speedup


The table below presents a comparison of the latencies of `arkmsm` and Arkworks
3.0 (baseline) on input sizes ranging from 2^8 to 2^18. Overall, `arkmsm`
achieved up to 2x speedup over the baseline.

| Input Size (ms)  | 2^8      | 2^9      | 2^10     | 2^11     | 2^12     | 2^13     | 2^14     | 2^15     | 2^16     | 2^17     | 2^18     |
|------------------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| Arkworks 3.0     | 13.903   | 24.208   | 37.5     | 67.545   | 121.03   | 204.92   | 375.85   | 693.46   | 1268.6   | 2324.9   | 4391.9   |
| Arkmsm           | 7.665    | 11.982   | 19.514   | 33.197   | 60.593   | 110.68   | 204.33   | 375.17   | 711.6    | 1372.1   | 2742.1   |
| Speedup          | 1.81     | 2.02     | 1.92     | 2.03     | 2.00     | 1.85     | 1.84     | 1.85     | 1.78     | 1.69     | 1.60     |

Performance measured on a AMD EPYC 7282 16-Core Processor.

## Optimizations

The following optimization techniques were used on top of the Pippenger Algorithm used in Arkworks.
1. Batch Addition in Bucket Accumulation (Batch Accumulation)
2. Batch Addition in Bucket Reduction (Batch Reduction)
3. Signed Bucket Indexes (Signed Index)
4. GLV Decomposition (GLV)

Each optimization is implemented in a separate commit so that the impact of each optimization can be accurately measured.

The table below presents a detailed breakdown of the performance improvements achieved with each optimization technique with a 2^12 intput size.

|                    |Latency (ms)| Improvement |
|--------------------|------------|-------------|
| Baseline           | 121.03     |             |
| Batch Accumulation | 92.154     | 31.33%      |
| Signed Index       | 74.429     | 23.81%      |
| Batch Reductio     | 65.023     | 14.47%      |
| GLV                | 60.593     | 7.31%       |
| Overall            |            | 99.74%      |

## Documentation

Detailed documentations about each optimization are available at [https://hackmd.io/@drouyang/msm](https://hackmd.io/@drouyang/msm)


## Getting Started

* Install [Rust](https://www.rust-lang.org/tools/install):
    ```bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env
    ```

* Run Unittests:
    ```bash
    cargo test
    ```

* Run Benchmarks:
    ```bash
    cargo bench
    ```

## Acknowledgement

The `arkmsm` project was funded by a grant from the [MINA
Foundation](https://minaprotocol.com).

Our algorithms and implementations were heavily based on 
the [2022 zPrize submission from Yrrid software](https://github.com/yrrid/submission-wasm-msm)
 written in C.

We would like to acknowledge [Gregor Mitscha-Baude](https://twitter.com/mitschabaude) from [O(1) Labs](https://twitter.com/o1_labs) 
and [Niall Emmart](https://www.linkedin.com/in/niall-emmart-0369384/) from [Yrrid Software](https://www.yrrid.com) 
for generously taking the time to answer our technical questions.

## Questions
For technical questions, please contact `ouyang at snarikify.io`

