# Gates' Error-Correcting Codes in Supersymmetry

This repository implements the groundbreaking connection between supersymmetry representations and error-correcting codes discovered by Gates et al. The implementation demonstrates how Adinkra graphs (used to visualize supersymmetry field relationships) naturally encode error-correcting code structures.

## Overview

In supersymmetry theory, **Adinkra graphs** are bipartite graphs that visualize the relationships between bosonic and fermionic fields through colored and dashed edges representing different supersymmetry operators. Gates and collaborators discovered that these structures exhibit remarkable similarities to linear block codes and other error-correcting codes.

### Key Concepts

- **Adinkra Graphs**: Bipartite graphs where nodes represent fields (bosonic/fermionic) and colored edges represent supersymmetry transformations
- **Linear Block Codes**: Error-correcting codes that map k-bit messages to n-bit codewords with systematic error detection/correction
- **Code Mapping**: The structural correspondence between Adinkra edge patterns and binary codewords

## Features

### Core Implementations

1. **Linear Block Codes** (`linear_block_codes.py`)
   - General linear block code class with encoding/decoding
   - Hamming codes, repetition codes, parity check codes
   - Syndrome computation and error detection
   - Minimum distance calculation

2. **Adinkra Graphs** (`adinkra_graphs.py`)
   - Adinkra graph construction and manipulation
   - Node types (bosonic/fermionic fields)
   - Colored and dashed edges
   - Adjacency matrix generation
   - Code mapping functionality

3. **Demonstration** (`gates_ecc_demo.py`)
   - Comprehensive examples showing the theory
   - Visual representation of mappings
   - Analysis of code properties
   - Physical interpretation

## Installation

```bash
pip install -r requirements.txt
```

### Requirements
- Python 3.7+
- NumPy >= 1.21.0
- Matplotlib >= 3.5.0

## Usage

### Basic Example

```python
from linear_block_codes import hamming_code
from adinkra_graphs import create_simple_adinkra

# Create a Hamming code
code = hamming_code(3)  # (7,4,3) Hamming code
print(f"Code parameters: {code.get_info()}")

# Create an Adinkra graph
adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
print(f"Adinkra has {len(adinkra.nodes)} nodes and {len(adinkra.edges)} edges")

# Convert Adinkra to error-correcting code
ecc_code = adinkra.to_linear_code()
print(f"Resulting code: {ecc_code}")
```

### Run Complete Demonstration

```bash
python gates_ecc_demo.py
```

This runs a comprehensive demonstration showing:
- Basic linear block codes
- Adinkra graph construction
- Mapping between Adinkra structures and codewords
- Hypercube-based Adinkra graphs
- Analysis of error-correcting properties
- Physical interpretation of the connections

## Mathematical Background

### Adinkra Graph Structure

An Adinkra graph is a bipartite graph G = (V_B ∪ V_F, E) where:
- V_B: Bosonic field vertices
- V_F: Fermionic field vertices  
- E: Colored and dashed edges representing SUSY generators

### Code Mapping

The mapping from Adinkra graphs to linear codes works as follows:

1. **Edge Structure → Binary Matrix**: Each color's adjacency matrix becomes a binary representation
2. **Dashing Pattern → Sign Structure**: Dashed edges encode sign information
3. **SUSY Constraints → Parity Checks**: Supersymmetry relations become code constraints

### Physical Significance

| Adinkra Element | Code Element | Physical Meaning |
|-----------------|--------------|------------------|
| Bosonic/Fermionic fields | Information/Parity bits | Field content |
| Edge colors | Generator polynomials | SUSY operators |
| Edge dashing | Sign patterns | Transformation signs |
| SUSY constraints | Parity equations | Physical consistency |

## Examples

### Simple 2×2 Adinkra

```python
# Create 2 bosonic + 2 fermionic fields with 2 SUSY generators
adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)

# Fields: φ₀, φ₁ (bosonic) | ψ₀, ψ₁ (fermionic)  
# Generators: Q₀ (RED), Q₁ (GREEN)

# Convert to code
code = adinkra.to_linear_code()
print(f"Code parameters: (n={code.n}, k={code.k}, d={code.d})")
```

### Hypercube Adinkra

```python
# Create hypercube-based Adinkra for higher dimensions
hypercube = create_hypercube_adinkra(dimension=3)

# 2³ = 8 vertices arranged as hypercube
# Alternating bosonic/fermionic based on vertex parity
print(f"Hypercube Adinkra: {len(hypercube.nodes)} nodes")
```

## Theory and Background

This implementation is based on the theoretical work showing that:

1. **Supersymmetry representations** naturally encode error-correcting structures
2. **Adinkra graphs** can be systematically mapped to linear block codes
3. **Code properties** reflect underlying supersymmetry constraints and stability
4. **Error correction** corresponds to constraint satisfaction in SUSY theories

### Key Insights

- Error correction capability ↔ Stability of field configurations
- Code rate ↔ Efficiency of field representation  
- Syndrome computation ↔ Detection of SUSY violations
- Minimum distance ↔ Robustness of supersymmetric structure

## Research Applications

This implementation enables research into:

- **Quantum Error Correction**: Connections between SUSY and quantum codes
- **AdS/CFT Correspondence**: Error correction in holographic theories
- **Supersymmetric Field Theory**: Computational tools for SUSY analysis
- **Information Theory**: Novel code constructions from physics

## File Structure

```
error-correcting-codes/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── linear_block_codes.py     # Linear block code implementations
├── adinkra_graphs.py         # Adinkra graph classes and functions
└── gates_ecc_demo.py         # Comprehensive demonstration script
```

## Contributing

Contributions are welcome! Areas for enhancement:

- Additional code families (BCH, Reed-Solomon, etc.)
- More sophisticated Adinkra constructions
- Visualization tools for graphs and codes
- Performance optimizations
- Extended physical interpretations

## References

- Gates, S.J. et al. "Adinkras and the Dynamics of Superspace Prepotentials" 
- "Error-correcting Codes and Supersymmetry Representations" (research papers)
- [WINLAB Research](https://www.winlab.rutgers.edu/) on supersymmetry and codes

## License

This project is open source. See individual files for specific licensing information.

---

*This implementation demonstrates the profound connections between supersymmetry theory and information theory, revealing how fundamental physics naturally incorporates error-correction structures.*