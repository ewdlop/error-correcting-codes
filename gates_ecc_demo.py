"""
Gates' Error-Correcting Codes Demonstration

This script demonstrates the connection between supersymmetry representations
(via Adinkra graphs) and error-correcting codes, following the work of 
Gates et al.

The demonstration shows how:
1. Adinkra graphs encode supersymmetry field relationships
2. These structures can be mapped to binary codewords
3. The resulting codewords form error-correcting codes
4. The code properties reflect the underlying supersymmetry constraints
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List
from linear_block_codes import LinearBlockCode, hamming_code, repetition_code
from adinkra_graphs import (
    AdinkraGraph, create_simple_adinkra, create_hypercube_adinkra,
    EdgeColor, FieldType
)


def demonstrate_basic_concepts():
    """Demonstrate basic error-correcting code concepts."""
    print("=" * 60)
    print("GATES' ERROR-CORRECTING CODES IN SUPERSYMMETRY")
    print("=" * 60)
    print()
    
    print("1. BASIC LINEAR BLOCK CODES")
    print("-" * 30)
    
    # Create a simple Hamming code
    hamming = hamming_code(3)
    print(f"Hamming Code: {hamming}")
    print(f"Generator Matrix:\n{hamming.G}")
    print(f"Code info: {hamming.get_info()}")
    print()
    
    # Demonstrate encoding
    message = [1, 0, 1, 1]
    codeword = hamming.encode(message)
    print(f"Message: {message}")
    print(f"Encoded: {codeword}")
    print(f"Is valid codeword: {hamming.is_codeword(codeword)}")
    print()


def demonstrate_adinkra_basics():
    """Demonstrate basic Adinkra graph concepts."""
    print("2. ADINKRA GRAPHS IN SUPERSYMMETRY")
    print("-" * 35)
    
    # Create a simple 2x2 Adinkra graph
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    print(f"Adinkra Graph (dimension {adinkra.dimension}):")
    print(f"  Bosonic nodes: {len(adinkra.bosonic_nodes)}")
    print(f"  Fermionic nodes: {len(adinkra.fermionic_nodes)}")
    print(f"  Total edges: {len(adinkra.edges)}")
    print(f"  Available colors: {[c.name for c in adinkra.available_colors]}")
    print()
    
    # Show node details
    print("Nodes:")
    for node in sorted(adinkra.nodes, key=lambda x: x.id):
        print(f"  {node.id}: {node.field_type.value} field at {node.position}")
    print()
    
    # Show edge details
    print("Edges (first 6):")
    for i, edge in enumerate(sorted(adinkra.edges, key=lambda x: (x.color.value, x.source.id))):
        if i >= 6:
            break
        dash_str = "dashed" if edge.dashed else "solid"
        print(f"  {edge.source.field_name} →[{edge.color.name}, {dash_str}]→ {edge.target.field_name}")
    print()
    
    # Check Adinkra rules
    rules = adinkra.check_adinkra_rules()
    print("Adinkra rule validation:")
    for rule, passed in rules.items():
        status = "✓" if passed else "✗"
        print(f"  {rule}: {status}")
    print()


def demonstrate_adinkra_to_code_mapping():
    """Demonstrate mapping from Adinkra graphs to error-correcting codes."""
    print("3. ADINKRA → ERROR-CORRECTING CODE MAPPING")
    print("-" * 45)
    
    # Create Adinkra graph
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    # Convert to codeword representation
    codewords = adinkra.to_codeword_representation()
    
    print("Codeword representation by color:")
    for color, codeword in codewords.items():
        print(f"  {color.name}: {codeword}")
    print()
    
    # Create corresponding linear code
    try:
        linear_code = adinkra.to_linear_code()
        print(f"Corresponding Linear Code: {linear_code}")
        print(f"Generator Matrix:")
        print(linear_code.G)
        print(f"Code parameters: {linear_code.get_info()}")
        print()
        
        # Show all codewords
        all_codewords = linear_code.get_all_codewords()
        print(f"All {len(all_codewords)} codewords:")
        for i, cw in enumerate(all_codewords):
            print(f"  {i:2d}: {cw}")
        print()
        
    except Exception as e:
        print(f"Error creating linear code: {e}")
        print()


def demonstrate_hypercube_adinkra():
    """Demonstrate hypercube-based Adinkra graphs."""
    print("4. HYPERCUBE ADINKRA GRAPHS")
    print("-" * 30)
    
    # Create hypercube Adinkra for different dimensions
    for dim in [2, 3]:
        print(f"Dimension {dim} Hypercube Adinkra:")
        
        hypercube = create_hypercube_adinkra(dim)
        print(f"  Nodes: {len(hypercube.nodes)} (2^{dim} hypercube vertices)")
        print(f"  Bosonic: {len(hypercube.bosonic_nodes)}")
        print(f"  Fermionic: {len(hypercube.fermionic_nodes)}")
        print(f"  Edges: {len(hypercube.edges)}")
        
        # Show field arrangement
        print("  Field arrangement:")
        for node in sorted(hypercube.nodes, key=lambda x: x.position):
            print(f"    {node.position}: {node.field_type.value} ({node.field_name})")
        
        # Get code parameters
        n, k = hypercube.get_code_parameters()
        print(f"  Corresponding code parameters: n={n}, k={k}")
        
        print()


def demonstrate_adjacency_matrices():
    """Show adjacency matrices for different colors."""
    print("5. ADJACENCY MATRICES BY COLOR")
    print("-" * 35)
    
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    for color in adinkra.available_colors:
        adj_matrix = adinkra.get_adjacency_matrix(color)
        print(f"{color.name} adjacency matrix:")
        print(adj_matrix)
        print(f"  Non-zero entries: {np.count_nonzero(adj_matrix)}")
        print(f"  Positive entries: {np.sum(adj_matrix > 0)}")
        print(f"  Negative entries: {np.sum(adj_matrix < 0)}")
        print()


def analyze_code_properties():
    """Analyze the error-correcting properties of Adinkra-derived codes."""
    print("6. ERROR-CORRECTING PROPERTIES ANALYSIS")
    print("-" * 42)
    
    # Compare different Adinkra configurations
    configs = [
        (2, 2, 2, "Simple 2x2, dim=2"),
        (3, 3, 2, "Simple 3x3, dim=2"), 
        (2, 2, 3, "Simple 2x2, dim=3"),
    ]
    
    print("Configuration comparison:")
    print("B×F×D  | Code    | Rate  | Min Dist | Error Corr | Description")
    print("-------|---------|-------|----------|------------|------------------")
    
    for n_bos, n_fer, dim, desc in configs:
        try:
            adinkra = create_simple_adinkra(n_bos, n_fer, dim)
            code = adinkra.to_linear_code()
            info = code.get_info()
            
            print(f"{n_bos}×{n_fer}×{dim}  | "
                  f"({info['n']},{info['k']},{info['d']}) | "
                  f"{info['rate']:.3f} | "
                  f"{info['d']:8d} | "
                  f"{info['error_correction_capability']:10d} | "
                  f"{desc}")
        except Exception as e:
            print(f"{n_bos}×{n_fer}×{dim}  | Error: {str(e)[:30]}...")
    
    print()


def demonstrate_physical_interpretation():
    """Explain the physical interpretation of the code structure."""
    print("7. PHYSICAL INTERPRETATION")
    print("-" * 28)
    print("""
The connection between Adinkra graphs and error-correcting codes reveals
deep structural relationships in supersymmetry:

ADINKRA GRAPH ELEMENTS → CODE ELEMENTS:
• Bosonic/Fermionic fields → Information/Parity bits
• Edge colors → Different generator polynomials
• Edge dashing → Sign patterns in codewords  
• SUSY constraints → Parity check equations
• Field transformations → Syndrome computations

PHYSICAL SIGNIFICANCE:
• Error correction ↔ Constraint satisfaction in SUSY
• Minimum distance ↔ Stability of field configurations
• Code rate ↔ Efficiency of field representation
• Syndrome ↔ Violation of supersymmetry relations

This mapping suggests that:
1. SUSY representations have inherent error-correcting structure
2. Quantum error correction may be fundamental to SUSY
3. Code properties reflect physical stability and constraints
""")


def visualize_simple_example():
    """Create a visual representation of a simple example."""
    print("8. SIMPLE VISUAL EXAMPLE")
    print("-" * 26)
    
    # Create minimal example
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    print("Consider a 2D supersymmetry with 2 bosonic + 2 fermionic fields:")
    print()
    print("Fields: φ₀(B), φ₁(B) | ψ₀(F), ψ₁(F)")
    print("SUSY generators: Q₀ (RED), Q₁ (GREEN)")
    print()
    
    # Show the structure
    print("Connections (showing dashing pattern):")
    for edge in sorted(adinkra.edges, key=lambda x: (x.color.value, x.source.id, x.target.id)):
        dash_symbol = "- -" if edge.dashed else "───"
        color_symbol = "R" if edge.color == EdgeColor.RED else "G"
        print(f"  {edge.source.field_name} {dash_symbol}[{color_symbol}]{dash_symbol} {edge.target.field_name}")
    print()
    
    # Show the resulting code
    try:
        code = adinkra.to_linear_code()
        print(f"Resulting error-correcting code: {code}")
        print("This code can detect/correct errors in the field relationships!")
    except Exception as e:
        print(f"Code generation: {e}")
    
    print()


def main():
    """Run the complete demonstration."""
    try:
        demonstrate_basic_concepts()
        demonstrate_adinkra_basics()
        demonstrate_adinkra_to_code_mapping()
        demonstrate_hypercube_adinkra()
        demonstrate_adjacency_matrices()
        analyze_code_properties()
        demonstrate_physical_interpretation()
        visualize_simple_example()
        
        print("=" * 60)
        print("DEMONSTRATION COMPLETE")
        print("=" * 60)
        print()
        print("This demonstrates Gates et al.'s discovery that supersymmetry")
        print("representations naturally encode error-correcting code structures,")
        print("revealing deep connections between physics and information theory.")
        
    except Exception as e:
        print(f"Error in demonstration: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()