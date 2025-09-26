"""
Linear Block Codes Implementation for Gates' Error-Correcting Codes Concept

This module implements basic linear block codes that serve as the foundation
for understanding the mapping between supersymmetry representations and 
error-correcting code structures as discovered by Gates et al.
"""

import numpy as np
from typing import List, Tuple, Optional
import itertools


class LinearBlockCode:
    """
    Represents a linear block code (n, k, d) where:
    - n: codeword length
    - k: message length (dimension)
    - d: minimum distance
    """
    
    def __init__(self, generator_matrix: np.ndarray):
        """
        Initialize with generator matrix G.
        
        Args:
            generator_matrix: k x n binary matrix where each row is a basis codeword
        """
        self.G = np.array(generator_matrix, dtype=int) % 2
        self.k, self.n = self.G.shape
        self.d = self._compute_minimum_distance()
        
        # Compute parity check matrix H
        self.H = self._compute_parity_check_matrix()
    
    def _compute_parity_check_matrix(self) -> np.ndarray:
        """Compute the parity check matrix H such that GH^T = 0."""
        # For systematic codes, if G = [I_k | P], then H = [-P^T | I_{n-k}]
        # For general case, we need to find the null space of G
        
        # Create augmented matrix and perform Gaussian elimination
        identity = np.eye(self.n, dtype=int)
        augmented = np.vstack([self.G, identity])
        
        # Gaussian elimination over GF(2)
        for i in range(self.k):
            # Find pivot
            pivot_row = i
            for j in range(i + 1, self.k):
                if augmented[j, i] == 1:
                    pivot_row = j
                    break
            
            if augmented[pivot_row, i] == 0:
                continue
                
            # Swap rows if needed
            if pivot_row != i:
                augmented[[i, pivot_row]] = augmented[[pivot_row, i]]
            
            # Eliminate
            for j in range(self.k):
                if j != i and augmented[j, i] == 1:
                    augmented[j] = (augmented[j] + augmented[i]) % 2
        
        # Extract the parity check matrix
        null_space_dim = self.n - self.k
        if null_space_dim > 0:
            H_rows = []
            for i in range(self.k, self.n):
                row = np.zeros(self.n, dtype=int)
                for j in range(self.n):
                    if j < self.k:
                        row[j] = augmented[j, i] if j < len(augmented) else 0
                    elif j == i:
                        row[j] = 1
                H_rows.append(row)
            
            if H_rows:
                return np.array(H_rows)
        
        # Fallback: create a basic parity check matrix
        if self.n > self.k:
            H = np.zeros((self.n - self.k, self.n), dtype=int)
            # Simple parity check: sum of all bits = 0
            if self.n - self.k >= 1:
                H[0, :] = 1
            return H
        else:
            return np.zeros((0, self.n), dtype=int)
    
    def _compute_minimum_distance(self) -> int:
        """Compute the minimum Hamming distance between distinct codewords."""
        min_dist = self.n + 1
        codewords = self.get_all_codewords()
        
        for i in range(len(codewords)):
            for j in range(i + 1, len(codewords)):
                dist = np.sum((codewords[i] + codewords[j]) % 2)  # XOR gives Hamming distance
                if dist < min_dist:
                    min_dist = dist
        
        return min_dist
    
    def encode(self, message: List[int]) -> np.ndarray:
        """
        Encode a message using the generator matrix.
        
        Args:
            message: k-bit message vector
            
        Returns:
            n-bit codeword
        """
        msg = np.array(message, dtype=int) % 2
        if len(msg) != self.k:
            raise ValueError(f"Message length must be {self.k}, got {len(msg)}")
        
        return (msg @ self.G) % 2
    
    def syndrome(self, received: List[int]) -> np.ndarray:
        """
        Compute syndrome of received word.
        
        Args:
            received: n-bit received vector
            
        Returns:
            (n-k)-bit syndrome vector
        """
        recv = np.array(received, dtype=int) % 2
        if len(recv) != self.n:
            raise ValueError(f"Received word length must be {self.n}, got {len(recv)}")
        
        return (recv @ self.H.T) % 2
    
    def get_all_codewords(self) -> List[np.ndarray]:
        """Generate all possible codewords."""
        codewords = []
        for message_bits in itertools.product([0, 1], repeat=self.k):
            codeword = self.encode(list(message_bits))
            codewords.append(codeword)
        return codewords
    
    def is_codeword(self, word: List[int]) -> bool:
        """Check if a word is a valid codeword."""
        syndrome = self.syndrome(word)
        return np.all(syndrome == 0)
    
    def get_info(self) -> dict:
        """Get information about this code."""
        return {
            'n': self.n,
            'k': self.k, 
            'd': self.d,
            'rate': self.k / self.n,
            'error_correction_capability': (self.d - 1) // 2,
            'error_detection_capability': self.d - 1
        }
    
    def __str__(self) -> str:
        return f"LinearBlockCode({self.n}, {self.k}, {self.d})"
    
    def __repr__(self) -> str:
        return self.__str__()


def hamming_code(r: int) -> LinearBlockCode:
    """
    Create a Hamming code with parameter r.
    
    Args:
        r: Number of parity bits (r >= 2)
        
    Returns:
        Hamming(2^r - 1, 2^r - r - 1, 3) code
    """
    if r < 2:
        raise ValueError("r must be at least 2")
    
    n = 2**r - 1
    k = n - r
    
    # Create parity check matrix H
    H = []
    for i in range(1, n + 1):
        binary_repr = format(i, f'0{r}b')
        H.append([int(b) for b in binary_repr])
    
    H = np.array(H).T
    
    # Create generator matrix G from H
    # Rearrange H to standard form [A | I]
    identity_cols = []
    for i in range(r):
        col = np.zeros(r, dtype=int)
        col[i] = 1
        for j in range(n):
            if np.array_equal(H[:, j], col):
                identity_cols.append(j)
                break
    
    # Create systematic generator matrix
    G = np.zeros((k, n), dtype=int)
    info_positions = [i for i in range(n) if i not in identity_cols]
    
    for i in range(k):
        G[i, info_positions[i]] = 1  # Information bits
        
    # Add parity bits
    for i in range(k):
        info_pos = info_positions[i]
        for j, parity_pos in enumerate(identity_cols):
            G[i, parity_pos] = H[j, info_pos]
    
    return LinearBlockCode(G)


def repetition_code(n: int) -> LinearBlockCode:
    """
    Create a repetition code that repeats each bit n times.
    
    Args:
        n: Number of repetitions (must be odd for error correction)
        
    Returns:
        Repetition code (n, 1, n)
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    
    G = np.ones((1, n), dtype=int)
    return LinearBlockCode(G)


def parity_check_code(k: int) -> LinearBlockCode:
    """
    Create a single parity check code.
    
    Args:
        k: Number of information bits
        
    Returns:
        Parity check code (k+1, k, 2)
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    
    G = np.zeros((k, k + 1), dtype=int)
    
    # Information bits
    for i in range(k):
        G[i, i] = 1
    
    # Parity bit (sum of all information bits)
    for i in range(k):
        G[i, k] = 1
    
    return LinearBlockCode(G)