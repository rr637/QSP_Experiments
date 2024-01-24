import numpy as np

def str_replace(string, char, index):
    assert(index >= 0 and index < len(string))
    output = string[0:index]
    output += char
    output += string[index + 1:]
    return output

#Requires [src] and [dest] to be integers that differ by exactly 1 bit,
# returns the index of that bit.
def find_target(src, dest):
    xor = src ^ dest
    output = xor.bit_length() - 1
    assert(1 << output == xor)
    return output

#Requires [src] and [dest] to be integers that differ by exactly 1 bit,
# returns (control, target) such that [src] and [dest] differ on the [target]
# bit and both src and dest have a 1 at the [control] bit. Also guarantees that
# control and target differ by exactly 1, and requires that such a
# (control, target) pair can be found.
def find_control_target(src, dest):
    target = find_target(src, dest)
    if (src >> (target + 1)) % 2: return target + 1, target
    return target - 1, target

#Takes [number], flips the [index]th bit, returns the result
def bit_flip(number, index):
    return number ^ (1 << index)

def get_bit(number, index):
    return (number >> index) % 2

def max_index(array):
    return np.where(array == max(array))[0][0]

def norm(v):
    return np.vdot(v, v) ** 0.5
