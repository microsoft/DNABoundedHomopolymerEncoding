# Bounded Homopolymer Encoding for DNA Storage

Codec for converting binary strings into DNA-like sequences over the alphabet `{A, C, G, T}` while guaranteeing that no symbol repeats more than *k* times consecutively (length *k* homopolymer run). The repository ships a fast C++ encoder/decoder that supports maximum homopolymer run lengths `k ∈ {1, 2, 3, 4, 5}`. See the [writeup](https://github.com/microsoft/DNABoundedHomopolymerEncoding/blob/main/BoundedHomopolymerEncoding.pdf) for how the encoding and decoding algorithms work. This code is written as part of the [DNA Storage Project](https://www.microsoft.com/en-us/research/project/dna-storage/) at Microsoft.

## Performance
- About 50Mbps for encoding and 80 Mbps for decoding on a single core. Can be easily parallelized across multiple cores.
- Max run length 1 is about 3 times faster, because it is basically base 2 to base 3 conversion. 
- Uses GMP multi precision library for maintaining large integers.


## Why bounded homopolymers?

DNA synthesis and sequencing are sensitive to long runs of identical bases (e.g., `AAA` or `CCCC`). By bounding the length of these homopolymer runs, stored strands become more robust to physical errors. The encoder in this repository ensures that every generated strand respects the chosen bound *k* while still allowing you to recover the original binary payload exactly.

## Features

- Deterministic encoding of binary strings into base-4 strings subject to a homopolymer constraint.
- Efficient decoder that returns the original binary message.
- Supports five constraint strengths (`k = 1…5`).
- Includes benchmarking harnesses to measure encoding/decoding throughput with GMP-backed big-integer arithmetic.

## Prerequisites

The C++ sources depend on:

- A C++17-compatible compiler (tested with `g++`).
- [GNU Multiple Precision Arithmetic Library (GMP)](https://gmplib.org/) with C++ bindings (`libgmpxx`).

### Installing GMP

```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install -y libgmp3-dev

# macOS using Homebrew
brew install gmp
```

For Windows, download the prebuilt binaries or build from source following the instructions at [gmplib.org](https://gmplib.org/) or the [NYU GMP guide](https://cs.nyu.edu/~exact/core/gmp/index.html).

## Building the encoder/decoder

Clone this repository and compile the main driver:

```bash
git clone https://github.com/microsoft/DNABoundedHomopolymerEncoding.git
cd DNABoundedHomopolymerEncoding
g++ -O3 BoundedHomopolymerEncoding.cpp -lgmpxx -lgmp -o bounded_homopolymer
```

This produces an executable named `bounded_homopolymer` that links against GMP.

If you also need the rate benchmarking tool, build `BHE_rates.cpp` in the same way:

```bash
g++ -O3 BHE_rates.cpp -lgmpxx -lgmp -o bhe_rates
```

## Usage

### Command line driver

The compiled `bounded_homopolymer` binary encodes random binary messages and checks that the decoder inverts the process.

```bash
./bounded_homopolymer <max_homopolymer_run_length> <encoding_length> <input_data_length> <number_trials>
```

- `max_homopolymer_run_length` (`k`): Choose from `1` to `5`. This bounds the length of any homopolymer run in the encoded strand. For example, `k = 3` guarantees that substrings such as `AAAA` or `GGGG` never appear; the maximum run is `AAA`, `CCC`, `GGG`, or `TTT`.
- `encoding_length`: Desired length (in bases) of the encoded strand.
- `input_data_length`: Length (in bits) of the binary message to encode.
- `number_trials`: Number of random messages to encode/decode for benchmarking.

Example run:

```bash
./bounded_homopolymer 3 25 40 100
```

The program prints the maximum encodable message size for the provided parameters and reports encoding/decoding timings. If any trial fails to round-trip, an error message is emitted.

### Estimating minimal encoding lengths

The companion `bhe_rates` tool leverages the same encoder to estimate the smallest strand length that can host a target number of data bits while honoring the homopolymer constraint. It evaluates each `k ∈ {1,…,5}` and reports the achievable rate (bits per base).

```bash
./bhe_rates <encoding_length>
```

- `encoding_length`: Strand length (in bases) you are willing to allocate.
- Output: For each `k ∈ {1,…,5}`, the tool reports the maximum payload size (in bits) that fits within `encoding_length` bases while respecting the homopolymer constraint, plus the corresponding rate.

Example:

```bash
./bhe_rates 96
```

Sample output:

```
k     max_input_bits    rate(bits/base)
------------------------------------------
1     152              1.583333
2     184              1.916667
3     190              1.979167
4     191              1.989583
5     191              1.989583
```

If a given constraint cannot store any payload at that length, the tool reports `0` bits for that row.

### Using the classes in your own code

Both `BoundedHomopolymerEncoder` and `BoundedHomopolymerDecoder` are self-contained in `BoundedHomopolymerEncoding.cpp`. You can include the file or extract the classes into your project.

```cpp
#include "BoundedHomopolymerEncoding.cpp"

int main() {
	const int k = 3;              // Maximum homopolymer run length
	const int encoded_len = 25;   // Length of the output DNA strand
	const int data_len = 12;      // Length of the binary input

	BoundedHomopolymerEncoder encoder(k, encoded_len, data_len);
	std::string input = "101101001011"; // Provide your own binary string
	std::string encoded = encoder.encode(input);
	// Map digits {0,1,2,3} to {A,C,G,T} (or your preferred ordering) to obtain a DNA strand.

	BoundedHomopolymerDecoder decoder(k, encoded_len, data_len);
	std::string recovered = decoder.decode(encoded);
}
```

## Troubleshooting

- **`Max data bits that can be encoded` warning:** Reduce `input_data_length` or increase `encoding_length` so the message fits the selected constraint.
- **Linker errors about GMP:** Confirm `libgmp` and `libgmpxx` are installed and search paths are visible to your compiler.

## Acknowledgements

The author would like to thank [Jon Schneider](https://jschnei.github.io/) for his insightful comments which motivated the encoding and decoding algorithms.


## Contributing

This project welcomes contributions and suggestions. Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit [Contributor License Agreements](https://cla.opensource.microsoft.com).

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft
trademarks or logos is subject to and must follow
[Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship.
Any use of third-party trademarks or logos are subject to those third-party's policies.
