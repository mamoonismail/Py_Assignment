#!/usr/bin/env python3

from seq_Complexity import count_kmers
import pytest

# To avoid repetition 
@pytest.fixture 

def seq():	
	seq = 'AGGATGAATGG'
	return seq

#test k = 2 produces 16 kmers 
def test_count_kmers(seq):
	counts = count_kmers(seq, 2)
	assert len(counts) == 6 
