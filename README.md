# qrusty
[![License](https://img.shields.io/github/license/Qiskit/qiskit-terra.svg?style=popout-square)](https://opensource.org/licenses/Apache-2.0)

**Qrusty** is a small library for working with quantum-computing in Rust and Python.

It is meant as a sandbox to experiment with code that is meant to be compatible with IBM's Qiskit.

## Installation

```
rustc --version
```
This should show version 1.61, mine prints 

```
>>> rustc 1.61.0 (fe5b13d68 2022-05-18)
```

Next we build using `maturin`

```
# go into pyqrusty dir
cd qrusty/pyqrusty
pip install maturin
maturin develop
```

## License

[Apache License 2.0](LICENSE.txt)
