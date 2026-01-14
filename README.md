# PyO3 + Maturin: Build Rust, Install as a Python Module

This README describes a minimal and reproducible workflow to compile a Rust library into a Python-importable extension module using **PyO3** and **maturin**, and then import it from Python.

---

## 1. Prerequisites

### System tools

- Rust toolchain (stable)
- Python 3.8 or newer (a virtual environment is recommended)
- A C/C++ build toolchain available on your OS  
  - macOS: Xcode Command Line Tools  
  - Windows: Visual Studio Build Tools (MSVC)  
  - Linux: GCC or Clang  

### Python tooling

Install maturin in your active Python environment:
```bash
python -m pip install -U pip maturin
```

Verify installations:
```bash
maturin --version
rustc --version
python --version
```

---

2. Create a New Project

Generate a PyO3 project scaffold:

```bash
maturin new -b pyo3 my_rust_py
cd my_rust_py
```

This creates:
	•	Cargo.toml
	•	src/lib.rs
	•	pyproject.toml

---

3. Configure Cargo.toml

Ensure the library is built as a dynamic library and that PyO3 is enabled.

[package]
name = "my_rust_py"
version = "0.1.0"
edition = "2021"

[lib]
name = "my_rust_py"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.21", features = ["extension-module"] }

Notes
	•	cdylib is required for Python extension modules.
	•	extension-module ensures correct linking against Python.

---

4. Configure pyproject.toml

A typical pyproject.toml for maturin:

[build-system]
requires = ["maturin>=1.5,<2.0"]
build-backend = "maturin"

[project]
name = "my-rust-py"
version = "0.1.0"
requires-python = ">=3.8"

[tool.maturin]
bindings = "pyo3"
module-name = "my_rust_py"

Important distinctions
	•	project.name is the Python package name.
	•	module-name is the name you import in Python.

---

5. Implement the Rust Module (src/lib.rs)

use pyo3::prelude::*;

/// Add two integers.
#[pyfunction]
fn add(a: i64, b: i64) -> i64 {
    a + b
}

/// Python module definition.
#[pymodule]
fn my_rust_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    Ok(())
}

Key points
	•	#[pyfunction] exposes a Rust function to Python.
	•	#[pymodule] defines the Python module and registers functions.

---

6. Build and Install (Development Mode)

Activate your virtual environment, then run:

```bash
maturin develop
```
maturin develop
```
If multiple Python interpreters are installed:

maturin develop -i python3.11

This compiles the extension and installs it into the active environment.

---

7. Test in Python

```bash
python -c "import my_rust_py; print(my_rust_py.add(2, 3))"
```

Expected output:

5


---

8. Build a Wheel for Distribution

```bash
maturin build --release
```

Install the wheel:

```bash
python -m pip install target/wheels/*.whl
```


---

9. Common Issues

Import name mismatch

If ModuleNotFoundError occurs, check:
	•	[tool.maturin].module-name in pyproject.toml
	•	The function name in #[pymodule] fn my_rust_py(...)

Wrong Python environment

Confirm the interpreter:

which python
python -c "import sys; print(sys.executable)"
python -m pip -V

macOS build errors

Install command line tools:

xcode-select --install

Windows build errors

Ensure MSVC Build Tools are installed and available.

---

10. Minimal Project Structure

my_rust_py/
  Cargo.toml
  pyproject.toml
  src/
    lib.rs


---

11. Command Summary

# Create project
maturin new -b pyo3 my_rust_py
cd my_rust_py

# Install into current Python environment
maturin develop

# Build wheels
maturin build --release

# Test import
python -c "import my_rust_py; print(my_rust_py.add(2, 3))"

If you would like, I can also provide:
- A version with PyO3 class bindings (`#[pyclass]`) instead of only functions.  
- A template for mixed Rust–Python workflows using `pyo3 + matplotlib` via callbacks.  
- A variant adapted for publishing to PyPI.