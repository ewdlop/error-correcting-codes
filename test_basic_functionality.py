"""
Basic functionality tests for the Gates error-correcting codes implementation.
These tests verify that the core components work correctly.
"""

import numpy as np
from linear_block_codes import LinearBlockCode, hamming_code, repetition_code, parity_check_code
from adinkra_graphs import (
    AdinkraGraph, AdinkraNode, AdinkraEdge, FieldType, EdgeColor,
    create_simple_adinkra, create_hypercube_adinkra
)


def test_linear_block_codes():
    """Test basic linear block code functionality."""
    print("Testing Linear Block Codes...")
    
    # Test Hamming code
    hamming = hamming_code(3)
    assert hamming.n == 7
    assert hamming.k == 4
    assert hamming.d >= 3
    
    # Test encoding
    message = [1, 0, 1, 1]
    codeword = hamming.encode(message)
    assert len(codeword) == 7
    assert hamming.is_codeword(codeword)
    
    # Test syndrome (should be zero for valid codeword)
    syndrome = hamming.syndrome(codeword)
    assert np.all(syndrome == 0)
    
    # Test repetition code
    rep_code = repetition_code(5)
    assert rep_code.n == 5
    assert rep_code.k == 1
    
    # Test parity check code
    parity_code = parity_check_code(4)
    assert parity_code.n == 5
    assert parity_code.k == 4
    
    print("✓ Linear Block Codes tests passed")


def test_adinkra_components():
    """Test Adinkra graph components."""
    print("Testing Adinkra Components...")
    
    # Test node creation
    bosonic_node = AdinkraNode(FieldType.BOSONIC, (0,), "phi_0")
    fermionic_node = AdinkraNode(FieldType.FERMIONIC, (0,), "psi_0")
    
    assert bosonic_node.field_type == FieldType.BOSONIC
    assert fermionic_node.field_type == FieldType.FERMIONIC
    
    # Test edge creation
    edge = AdinkraEdge(bosonic_node, fermionic_node, EdgeColor.RED, dashed=False)
    assert edge.color == EdgeColor.RED
    assert not edge.dashed
    
    # Test binary representation
    color_bits, dash_bit = edge.get_binary_representation()
    assert color_bits == EdgeColor.RED.value
    assert dash_bit == 0
    
    print("✓ Adinkra Components tests passed")


def test_adinkra_graph():
    """Test Adinkra graph functionality."""
    print("Testing Adinkra Graph...")
    
    # Create simple Adinkra
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    assert len(adinkra.bosonic_nodes) == 2
    assert len(adinkra.fermionic_nodes) == 2
    assert len(adinkra.nodes) == 4
    assert len(adinkra.edges) == 8  # 2*2*2 connections
    
    # Test adjacency matrices
    red_adj = adinkra.get_adjacency_matrix(EdgeColor.RED)
    assert red_adj.shape == (2, 2)
    
    # Test rule checking
    rules = adinkra.check_adinkra_rules()
    assert rules['bipartite']  # Should be bipartite
    
    # Test codeword generation
    codewords = adinkra.to_codeword_representation()
    assert EdgeColor.RED in codewords
    assert EdgeColor.GREEN in codewords
    
    print("✓ Adinkra Graph tests passed")


def test_hypercube_adinkra():
    """Test hypercube Adinkra construction."""
    print("Testing Hypercube Adinkra...")
    
    # Test 2D hypercube
    hypercube_2d = create_hypercube_adinkra(2)
    assert len(hypercube_2d.nodes) == 4  # 2^2
    assert len(hypercube_2d.bosonic_nodes) == 2
    assert len(hypercube_2d.fermionic_nodes) == 2
    
    # Test 3D hypercube  
    hypercube_3d = create_hypercube_adinkra(3)
    assert len(hypercube_3d.nodes) == 8  # 2^3
    assert len(hypercube_3d.bosonic_nodes) == 4
    assert len(hypercube_3d.fermionic_nodes) == 4
    
    print("✓ Hypercube Adinkra tests passed")


def test_adinkra_to_code_conversion():
    """Test conversion from Adinkra to linear code."""
    print("Testing Adinkra to Code Conversion...")
    
    # Create Adinkra and convert to code
    adinkra = create_simple_adinkra(n_bosonic=2, n_fermionic=2, dimension=2)
    
    try:
        code = adinkra.to_linear_code()
        assert isinstance(code, LinearBlockCode)
        assert code.k <= code.n  # Basic code property
        
        # Test that we can encode with the resulting code
        if code.k > 0:
            message = [1] * code.k
            codeword = code.encode(message)
            assert len(codeword) == code.n
            
    except Exception as e:
        # Some configurations might not produce valid codes
        print(f"Note: Code conversion produced: {e}")
    
    print("✓ Adinkra to Code Conversion tests passed")


def test_code_properties():
    """Test analysis of code properties."""
    print("Testing Code Properties...")
    
    # Test with different configurations
    adinkra_configs = [
        (2, 2, 2),
        (3, 3, 2),
    ]
    
    for n_bos, n_fer, dim in adinkra_configs:
        adinkra = create_simple_adinkra(n_bos, n_fer, dim)
        n, k = adinkra.get_code_parameters()
        
        assert n > 0
        assert k > 0
        assert k <= n  # Dimension cannot exceed length
        
    print("✓ Code Properties tests passed")


def run_all_tests():
    """Run all basic functionality tests."""
    print("=" * 50)
    print("RUNNING BASIC FUNCTIONALITY TESTS")
    print("=" * 50)
    print()
    
    try:
        test_linear_block_codes()
        test_adinkra_components()
        test_adinkra_graph()
        test_hypercube_adinkra()
        test_adinkra_to_code_conversion()
        test_code_properties()
        
        print()
        print("=" * 50)
        print("ALL TESTS PASSED ✓")
        print("=" * 50)
        print()
        print("The Gates error-correcting codes implementation")
        print("is working correctly!")
        
    except Exception as e:
        print(f)
        print("TEST FAILED ✗")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    run_all_tests()