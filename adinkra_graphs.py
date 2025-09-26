"""
Adinkra Graphs Implementation for Gates' Error-Correcting Codes

This module implements Adinkra graphs used in supersymmetry representations
and demonstrates how they can be mapped to error-correcting code structures,
following the work of Gates et al.

Adinkra graphs visualize the relationship between bosonic and fermionic fields
through colored and dashed edges representing different supersymmetry operators.
"""

import numpy as np
from typing import Dict, List, Tuple, Set, Optional, Union
from dataclasses import dataclass
from enum import Enum
import itertools
from linear_block_codes import LinearBlockCode


class FieldType(Enum):
    """Types of fields in supersymmetry."""
    BOSONIC = "B"
    FERMIONIC = "F"


class EdgeColor(Enum):
    """Colors for Adinkra graph edges representing different SUSY generators."""
    RED = 0
    GREEN = 1  
    BLUE = 2
    YELLOW = 3
    # Can extend for higher dimensions


@dataclass
class AdinkraNode:
    """Represents a node in an Adinkra graph."""
    field_type: FieldType
    position: Tuple[int, ...]  # Position in the lattice/grid
    field_name: str
    
    def __post_init__(self):
        self.id = f"{self.field_type.value}_{self.field_name}_{self.position}"
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        return isinstance(other, AdinkraNode) and self.id == other.id


@dataclass
class AdinkraEdge:
    """Represents an edge in an Adinkra graph."""
    source: AdinkraNode
    target: AdinkraNode
    color: EdgeColor
    dashed: bool = False  # Dashing represents sign flip/negative
    
    def __post_init__(self):
        # Ensure edges connect different field types
        if self.source.field_type == self.target.field_type:
            raise ValueError("Adinkra edges must connect bosonic and fermionic fields")
        # Create a unique ID for hashing
        self.id = f"{self.source.id}_{self.target.id}_{self.color.name}_{self.dashed}"
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        return isinstance(other, AdinkraEdge) and self.id == other.id
    
    def get_binary_representation(self) -> Tuple[int, int]:
        """
        Get binary representation of the edge for code mapping.
        
        Returns:
            (color_bits, dash_bit) where color_bits encode the color
            and dash_bit encodes the dashing state
        """
        color_value = self.color.value
        dash_bit = 1 if self.dashed else 0
        return (color_value, dash_bit)


class AdinkraGraph:
    """
    Represents an Adinkra graph for supersymmetry visualization.
    
    An Adinkra graph is a bipartite graph where:
    - One set of vertices represents bosonic fields
    - Another set represents fermionic fields  
    - Edges are colored and may be dashed
    - Each color represents a different supersymmetry generator
    """
    
    def __init__(self, dimension: int = 2):
        """
        Initialize an Adinkra graph.
        
        Args:
            dimension: Number of supersymmetry generators (colors)
        """
        self.dimension = dimension
        self.nodes: Set[AdinkraNode] = set()
        self.edges: Set[AdinkraEdge] = set()
        self.bosonic_nodes: Set[AdinkraNode] = set()
        self.fermionic_nodes: Set[AdinkraNode] = set()
        
        # Available colors based on dimension
        self.available_colors = list(EdgeColor)[:dimension]
    
    def add_node(self, node: AdinkraNode) -> None:
        """Add a node to the graph."""
        self.nodes.add(node)
        if node.field_type == FieldType.BOSONIC:
            self.bosonic_nodes.add(node)
        else:
            self.fermionic_nodes.add(node)
    
    def add_edge(self, edge: AdinkraEdge) -> None:
        """Add an edge to the graph."""
        if edge.color not in self.available_colors:
            raise ValueError(f"Color {edge.color} not available for dimension {self.dimension}")
        
        # Ensure both nodes are in the graph
        self.add_node(edge.source)
        self.add_node(edge.target)
        self.edges.add(edge)
    
    def get_adjacency_matrix(self, color: EdgeColor) -> np.ndarray:
        """
        Get adjacency matrix for a specific color.
        
        Args:
            color: Edge color to consider
            
        Returns:
            Adjacency matrix where entry (i,j) represents edge from 
            bosonic node i to fermionic node j
        """
        bosonic_list = sorted(list(self.bosonic_nodes), key=lambda x: x.id)
        fermionic_list = sorted(list(self.fermionic_nodes), key=lambda x: x.id)
        
        n_bosonic = len(bosonic_list)
        n_fermionic = len(fermionic_list)
        
        adj_matrix = np.zeros((n_bosonic, n_fermionic), dtype=int)
        
        for edge in self.edges:
            if edge.color != color:
                continue
                
            try:
                i = bosonic_list.index(edge.source if edge.source.field_type == FieldType.BOSONIC else edge.target)
                j = fermionic_list.index(edge.target if edge.target.field_type == FieldType.FERMIONIC else edge.source)
                
                # Use dashing to determine sign (1 for solid, -1 for dashed in matrix representation)
                value = -1 if edge.dashed else 1
                adj_matrix[i, j] = value
                
            except ValueError:
                continue
        
        return adj_matrix
    
    def to_codeword_representation(self) -> Dict[EdgeColor, np.ndarray]:
        """
        Convert Adinkra graph structure to codeword representation.
        
        This maps the edge structure to binary vectors that can be interpreted
        as codewords in an error-correcting code.
        
        Returns:
            Dictionary mapping each color to its binary representation
        """
        codewords = {}
        
        for color in self.available_colors:
            # Get adjacency matrix for this color
            adj_matrix = self.get_adjacency_matrix(color)
            
            # Flatten and convert to binary (taking absolute value for structure, 
            # separate tracking of dashing)
            structure_bits = (np.abs(adj_matrix).flatten() > 0).astype(int)
            
            # Get dashing pattern
            dashing_bits = []
            for edge in self.edges:
                if edge.color == color:
                    dashing_bits.append(1 if edge.dashed else 0)
            
            # Pad dashing_bits to match structure if needed
            while len(dashing_bits) < len(structure_bits):
                dashing_bits.append(0)
            
            dashing_array = np.array(dashing_bits[:len(structure_bits)])
            
            # Combine structure and dashing information
            combined = np.concatenate([structure_bits, dashing_array])
            codewords[color] = combined
            
        return codewords
    
    def check_adinkra_rules(self) -> Dict[str, bool]:
        """
        Check if the graph satisfies Adinkra graph rules.
        
        Returns:
            Dictionary of rule checks and their results
        """
        results = {}
        
        # Rule 1: Bipartite structure (bosonic-fermionic connections only)
        bipartite_valid = True
        for edge in self.edges:
            if edge.source.field_type == edge.target.field_type:
                bipartite_valid = False
                break
        results['bipartite'] = bipartite_valid
        
        # Rule 2: Each node should have edges of all available colors
        color_completeness = True
        for node in self.nodes:
            node_colors = set()
            for edge in self.edges:
                if edge.source == node or edge.target == node:
                    node_colors.add(edge.color)
            if len(node_colors) != len(self.available_colors):
                color_completeness = False
                break
        results['color_completeness'] = color_completeness
        
        # Rule 3: Dashing pattern should satisfy certain constraints
        # (This is a simplified check - actual constraints depend on specific representation)
        dashing_consistent = True
        results['dashing_consistent'] = dashing_consistent
        
        return results
    
    def get_code_parameters(self) -> Tuple[int, int]:
        """
        Get parameters for the corresponding linear code.
        
        Returns:
            (n, k) where n is codeword length and k is dimension
        """
        n_bosonic = len(self.bosonic_nodes)
        n_fermionic = len(self.fermionic_nodes)
        
        # Codeword length includes structure and dashing information
        n = n_bosonic * n_fermionic * 2  # Factor of 2 for structure + dashing
        
        # Dimension related to number of independent constraints
        k = len(self.available_colors)
        
        return n, k
    
    def to_linear_code(self) -> LinearBlockCode:
        """
        Convert the Adinkra graph to a linear block code.
        
        This creates a generator matrix where each row corresponds to 
        a color's codeword representation.
        """
        codewords = self.to_codeword_representation()
        
        if not codewords:
            raise ValueError("No codewords generated from empty graph")
        
        # Stack codewords as rows of generator matrix
        generator_rows = []
        for color in sorted(self.available_colors, key=lambda x: x.value):
            if color in codewords:
                generator_rows.append(codewords[color])
        
        if not generator_rows:
            raise ValueError("No valid generator rows created")
        
        # Ensure all rows have the same length
        max_length = max(len(row) for row in generator_rows)
        padded_rows = []
        for row in generator_rows:
            padded_row = np.zeros(max_length, dtype=int)
            padded_row[:len(row)] = row
            padded_rows.append(padded_row)
        
        generator_matrix = np.array(padded_rows)
        return LinearBlockCode(generator_matrix)


def create_simple_adinkra(n_bosonic: int = 2, n_fermionic: int = 2, dimension: int = 2) -> AdinkraGraph:
    """
    Create a simple Adinkra graph for demonstration.
    
    Args:
        n_bosonic: Number of bosonic fields
        n_fermionic: Number of fermionic fields  
        dimension: Number of supersymmetry generators
        
    Returns:
        AdinkraGraph with basic structure
    """
    graph = AdinkraGraph(dimension)
    
    # Create bosonic nodes
    bosonic_nodes = []
    for i in range(n_bosonic):
        node = AdinkraNode(FieldType.BOSONIC, (i,), f"phi_{i}")
        bosonic_nodes.append(node)
        graph.add_node(node)
    
    # Create fermionic nodes  
    fermionic_nodes = []
    for i in range(n_fermionic):
        node = AdinkraNode(FieldType.FERMIONIC, (i,), f"psi_{i}")
        fermionic_nodes.append(node)
        graph.add_node(node)
    
    # Add edges with different colors
    available_colors = graph.available_colors
    
    for i, bosonic in enumerate(bosonic_nodes):
        for j, fermionic in enumerate(fermionic_nodes):
            for k, color in enumerate(available_colors):
                # Create edge with some pattern for dashing
                dashed = (i + j + k) % 2 == 1  # Simple pattern
                edge = AdinkraEdge(bosonic, fermionic, color, dashed)
                graph.add_edge(edge)
    
    return graph


def create_hypercube_adinkra(dimension: int) -> AdinkraGraph:
    """
    Create an Adinkra graph based on hypercube structure.
    
    This represents a common structure in supersymmetry where fields
    are arranged on a hypercube with edges representing SUSY transformations.
    
    Args:
        dimension: Dimension of the hypercube (number of SUSY generators)
        
    Returns:
        AdinkraGraph representing hypercube structure
    """
    graph = AdinkraGraph(dimension)
    
    # Generate all vertices of the hypercube (2^dimension vertices)
    vertices = list(itertools.product([0, 1], repeat=dimension))
    
    # Create nodes - alternate between bosonic and fermionic based on parity
    nodes = {}
    for vertex in vertices:
        parity = sum(vertex) % 2
        field_type = FieldType.BOSONIC if parity == 0 else FieldType.FERMIONIC
        field_name = ''.join(map(str, vertex))
        node = AdinkraNode(field_type, vertex, field_name)
        nodes[vertex] = node
        graph.add_node(node)
    
    # Add edges along hypercube edges
    available_colors = graph.available_colors
    
    for vertex in vertices:
        for i in range(dimension):
            # Create adjacent vertex by flipping i-th bit
            adjacent = list(vertex)
            adjacent[i] = 1 - adjacent[i]
            adjacent = tuple(adjacent)
            
            if adjacent in nodes:
                source = nodes[vertex] 
                target = nodes[adjacent]
                color = available_colors[i]
                
                # Dashing pattern based on specific rules
                # (simplified - in practice depends on representation theory)
                dashed = vertex[i] == 1  
                
                edge = AdinkraEdge(source, target, color, dashed)
                graph.add_edge(edge)
    
    return graph