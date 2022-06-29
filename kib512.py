import calc_consts_kib512
import numpy as np

""" pre-processing"""
# user input as np.uint8 matrix

# multiple of 8 by 8? 512-bit input. 

# how will the multi-block processing work? xor them, add them?

# flag indicating end of data

# flag indicating end of block

# endiennes of message

""" compression word matrix function """
# define size of single block and how to compress each block

# 4 by 4 matrix for each.
# top and bottom row and second top and second below of the input matrix

""" compression function """
# how many rounds for compression? 

# matrix multipication with bitmasking

# have mutiple mini compression functions for the main one

# value depending on size? second parameter for size of hash? this would mean the output isn't only 512-bits. Check BLAKE2
