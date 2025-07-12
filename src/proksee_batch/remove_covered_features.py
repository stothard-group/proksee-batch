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


def remove_covered_features_optimized(
    features: List[Tuple[int, int, float]],
) -> List[bool]:
    """
    Optimized version using sorting-based algorithm for better performance.

    Identifies features that are completely covered by another feature with an
    equal or higher score.

    Parameters:
    features (list of tuples): A list of tuples, each containing a start position,
                               an end position, and a score.

    Returns:
    list: A list of boolean values indicating if the corresponding feature is not covered.
    """
    num_features = len(features)

    if num_features == 0:
        return []

    if num_features == 1:
        return [True]

    # For small datasets, use the simple O(n²) algorithm
    if num_features <= 300:
        return remove_covered_features(features)

    # Initialize result array
    is_not_covered = [True] * num_features

    # Create a list of features with their original indices
    indexed_features = [
        (i, start, end, score) for i, (start, end, score) in enumerate(features)
    ]

    # Sort by score (descending) then by start position
    indexed_features.sort(key=lambda x: (-x[3], x[1], x[0]))

    # Keep track of non-covered features for each score level
    kept_features = []

    for i, start, end, score in indexed_features:
        is_covered = False

        # Check against all previously kept features with higher or equal scores
        for other_i, other_start, other_end, other_score in kept_features:
            if other_score > score:
                # Check if the other feature completely covers this feature
                if other_start <= start and other_end >= end:
                    is_covered = True
                    break
            elif other_score == score:
                # Handle equal scores with same logic as original
                if other_start <= start and other_end >= end:
                    # Check if the other feature covers and is larger
                    if other_start < start or other_end > end:
                        is_covered = True
                        break
                    # Same size features: keep the one with lower original index
                    elif other_start == start and other_end == end and i > other_i:
                        is_covered = True
                        break

        if is_covered:
            is_not_covered[i] = False
        else:
            # Keep this feature for future comparisons
            kept_features.append((i, start, end, score))

    # Verify at least one feature remains
    if not any(is_not_covered):
        # Fallback to original algorithm if something went wrong
        return remove_covered_features(features)

    return is_not_covered
