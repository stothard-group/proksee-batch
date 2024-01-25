from typing import List
from typing import Tuple


def remove_covered_features(features: List[Tuple[int, int, int]]) -> List[bool]:
    """
    Identifies features that are completely covered by another feature with an
    equal or higher score.

    Parameters:
    features (list of tuples): A list of tuples, each containing a start position,
                               an end position, and a score.

    Returns:
    list: A list of boolean values indicating if the corresponding feature is not covered.
    """

    # Initialize an array to keep track of whether each feature is covered
    is_not_covered = [True] * len(features)

    # Iterate through each feature
    for i, (start, end, score) in enumerate(features):
        # Check for any other feature that covers this feature
        for j, (other_start, other_end, other_score) in enumerate(features):
            if i != j and other_score >= score:
                # Check if the other feature completely covers the current feature
                if other_start <= start and other_end >= end:
                    is_not_covered[i] = False
                    break

    return is_not_covered
