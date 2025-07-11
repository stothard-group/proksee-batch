from typing import List
from typing import Tuple


def remove_covered_features(features: List[Tuple[int, int, float]]) -> List[bool]:
    """
    Identifies features that are completely covered by another feature with an
    equal or higher score.

    Parameters:
    features (list of tuples): A list of tuples, each containing a start position,
                               an end position, and a score.

    Returns:
    list: A list of boolean values indicating if the corresponding feature is not covered.
    """
    # Count the number of input features.
    num_features = len(features)

    assert num_features > 0, "The input list of features is empty."

    # Validate input types
    for start, end, score in features:
        # Type checking is already enforced by type hints
        # Ensure score is numeric (this is redundant with type hints but kept for runtime safety)
        pass

    # Initialize an array to keep track of whether each feature is covered
    is_not_covered = [True] * len(features)

    # Iterate through each feature
    for i, (start, end, score) in enumerate(features):
        # Check for any other feature that covers this feature
        for j, (other_start, other_end, other_score) in enumerate(features):
            if i != j:
                if other_score > score:
                    # Check if the other feature completely covers the current feature
                    if other_start <= start and other_end >= end:
                        is_not_covered[i] = False
                        break
                elif other_score == score:
                    if other_start <= start and other_end >= end:
                        # Check if the other feature covers the current feature and is larger.
                        if other_start < start or other_end > end:
                            is_not_covered[i] = False
                            break
                        # Check if the other feature covers the current feature and is the same size.
                        elif other_start == start and other_end == end and i > j:
                            is_not_covered[i] = False
                            break

    assert (
        len(is_not_covered) == num_features
    ), "The number of output values does not match the number of input features."

    # Check that there is at least one feature that is not covered.
    assert any(
        is_not_covered
    ), "All features are covered by other features. This can't be right."

    return is_not_covered
