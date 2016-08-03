"""
Collection of function used to aid in text field processing in BAC and related
scripts
"""
from __future__ import print_function


def remove_digits(input_str):
    """
    Remove all digits from input string.

    Args:
        input_str (str): Input string

    Returns:
        str: Input text stripped of digits

    """

    return ''.join([i for i in input_str if not i.isdigit()])