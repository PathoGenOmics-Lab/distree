# Installation

distree is a single self-contained binary. It has no runtime dependencies, no
interpreter, and no configuration file.

## From Bioconda

=== "conda"

    ```bash
    conda install -c bioconda distree
    ```

=== "mamba"

    ```bash
    mamba install -c bioconda distree
    ```

=== "pixi"

    ```bash
    pixi add --channel bioconda distree
    ```

## From a release binary

Each release publishes binaries for Linux and macOS, on both x86-64 and arm64,
on the [releases page](https://github.com/PathoGenOmics-Lab/distree/releases).
Download the one for your platform, make it executable, and put it on your
`PATH`:

```bash
curl -LO https://github.com/PathoGenOmics-Lab/distree/releases/latest/download/distree-linux-x86_64
chmod +x distree-linux-x86_64
mv distree-linux-x86_64 ~/.local/bin/distree
```

The four assets are `distree-linux-x86_64`, `distree-linux-aarch64`,
`distree-macos-x86_64` and `distree-macos-aarch64`.

!!! note "macOS Gatekeeper"

    The binaries are not notarised, so macOS quarantines a downloaded one. Clear
    the flag with `xattr -d com.apple.quarantine distree` before the first run.

## From source

Any stable Rust toolchain will do; [rustup](https://rustup.rs) is the usual way
to get one.

```bash
git clone https://github.com/PathoGenOmics-Lab/distree.git
cd distree
cargo build --release
```

The binary lands at `target/release/distree`. Copy it wherever you keep local
tools, or run it in place.

## Checking the install

```bash
distree --version
echo '((A:1.0,B:2.0):0.5,C:3.0);' | distree -
```

```
	A	B	C
A	0.0000000000	3.0000000000	4.5000000000
B	3.0000000000	0.0000000000	5.5000000000
C	4.5000000000	5.5000000000	0.0000000000
```

A tab-separated square matrix, tips sorted alphabetically, with an empty cell in
the top-left corner so the header lines up over the data columns.

## Next

[Quick start](quickstart.md) covers the commands most runs actually use.
