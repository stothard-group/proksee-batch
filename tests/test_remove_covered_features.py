import random
import time
from typing import List
from typing import Tuple

from proksee_batch.remove_covered_features import remove_covered_features
from proksee_batch.remove_covered_features import remove_covered_features_optimized


def test_remove_covered_features() -> None:
    """Test the remove_covered_features function. It should take a list of
    tuples. Each tuple should have three elements: a start position, an end
    position, and a score. The function should return a list of boolean values.
    The boolean values should indicate whether the corresponding tuple is
    completely covered by another tuple with an equal or higher score. True
    means it is not covered, False means it is covered. Positions are 1-based,
    inclusive.
    """
    # Test case 1: No overlap
    # 1-10, 11-20, 21-30
    # 1-10 is not covered
    # 11-20 is not covered
    # 21-30 is not covered
    # Expected result: [True, True, True]
    test_list = [(1, 10, 1.0), (11, 20, 2.0), (21, 30, 3.0)]
    expected_result = [True, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result
    test_list = [(1, 10, 3.0), (11, 20, 2.0), (21, 30, 1.0)]
    expected_result = [True, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 2: Overlap
    # 1-10, 1-15, 20-30
    # 1-10 is covered
    # 5-15 is not covered
    # 20-30 is not covered
    # Expected result: [False, True, True]
    test_list = [(1, 10, 1.0), (1, 15, 2.0), (20, 30, 1.0)]
    expected_result = [False, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 3: Overlap
    # 1-10, 1-15, 1-20
    # 1-10 is covered
    # 1-15 is covered
    # 1-20 is not covered
    # Expected result: [False, False, True]
    test_list = [(1, 10, 1.0), (1, 15, 2.0), (1, 20, 3.0)]
    expected_result = [False, False, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 4: Overlap
    # 1-20, 5-20, 10-20
    # 1-20 is not covered
    # 5-20 is covered
    # 10-20 is covered
    # Expected result: [True, False, False]
    test_list = [(1, 20, 3.0), (5, 20, 2.0), (10, 20, 1.0)]
    expected_result = [True, False, False]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 5: Overlap
    # 1-20, 5-20, 10-20
    # 1-20 is not covered
    # 5-20 is not covered
    # 10-20 is not covered
    # Expected result: [True, True, True]
    test_list = [(1, 20, 1.0), (5, 20, 2.0), (10, 20, 3.0)]
    expected_result = [True, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 6: Multiple nested overlaps
    # 1-20, 2-19, 3-18, 4-17, 5-16, 6-15, 7-14, 8-13, 9-12, 10-11
    # 1-20 is not covered
    # 10-11 is not covered
    # All other ranges are covered by 1-20
    # Expected result: [True, False, False, False, False, False, False, False, False, True]
    test_list = [
        (1, 20, 10.0),
        (2, 19, 2.0),
        (3, 18, 3.0),
        (4, 17, 4.0),
        (5, 16, 5.0),
        (6, 15, 6.0),
        (7, 14, 7.0),
        (8, 13, 8.0),
        (9, 12, 9.0),
        (10, 11, 11.0),
    ]
    expected_result = [
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
    ]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 7: Multiple features with the same locations but different scores
    # 1-10, 1-10, 1-10
    # One feature is not covered
    # Expected result: [False, False, True]
    test_list = [(1, 10, 1.0), (1, 10, 2.0), (1, 10, 3.0)]
    expected_result = [False, False, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 8: Multiple features with the same locations and scores
    # 1-10, 1-10, 1-10
    # One feature is not covered
    # Expected result: [True, False, False]
    test_list = [(1, 10, 1.0), (1, 10, 1.0), (1, 10, 1.0)]
    expected_result = [True, False, False]
    result = remove_covered_features(test_list)
    assert result == expected_result


def generate_random_features(
    n: int, max_pos: int = 10000, max_score: float = 100.0
) -> List[Tuple[int, int, float]]:
    """Generate random features for performance testing."""
    features = []
    for _ in range(n):
        start = random.randint(1, max_pos - 100)
        end = random.randint(start + 1, min(start + 500, max_pos))
        score = random.uniform(0.1, max_score)
        features.append((start, end, score))
    return features


def generate_nested_features(n: int) -> List[Tuple[int, int, float]]:
    """Generate features with known nesting patterns for correctness testing."""
    features = []
    base_start = 1
    base_end = 1000

    # Add a large encompassing feature with highest score
    features.append((base_start, base_end, 100.0))

    # Add nested features with decreasing scores
    for i in range(1, n):
        start = base_start + i * 10
        end = base_end - i * 10
        score = 100.0 - i
        if start < end:
            features.append((start, end, score))

    return features


def test_remove_covered_features_performance() -> None:
    """Test performance characteristics and ensure reasonable execution time."""
    # Test with progressively larger datasets
    sizes = [10, 50, 100, 200]

    for size in sizes:
        # Generate random features
        features = generate_random_features(size)

        # Measure execution time
        start_time = time.time()
        result = remove_covered_features(features)
        end_time = time.time()

        execution_time = end_time - start_time

        # Basic correctness checks
        assert len(result) == len(features), f"Result length mismatch for size {size}"
        assert any(result), f"All features marked as covered for size {size}"
        assert all(
            isinstance(x, bool) for x in result
        ), f"Non-boolean results for size {size}"

        # Performance check - should complete reasonably quickly
        # For O(n²) algorithm, even 200 features should complete in under 1 second
        assert (
            execution_time < 1.0
        ), f"Performance too slow for size {size}: {execution_time:.3f}s"

        print(
            f"Size {size}: {execution_time:.4f}s, kept {sum(result)}/{len(features)} features"
        )


def test_remove_covered_features_nested_correctness() -> None:
    """Test correctness with known nested patterns."""
    # Generate features with predictable nesting
    features = generate_nested_features(10)
    result = remove_covered_features(features)

    # First feature (largest, highest score) should never be covered
    assert result[0] is True, "Largest feature with highest score was marked as covered"

    # Most nested features should be covered by the encompassing one
    covered_count = sum(1 for x in result if not x)
    assert covered_count > 0, "No features were marked as covered in nested test"

    # Last feature (smallest, lowest score) should be covered
    assert result[-1] is False, "Smallest nested feature was not marked as covered"


def test_remove_covered_features_edge_cases() -> None:
    """Test edge cases that might cause performance issues."""

    # All identical features
    identical_features = [(100, 200, 50.0)] * 20
    result = remove_covered_features(identical_features)
    # Only first should remain (due to index-based tie breaking)
    assert result[0] is True
    assert all(not x for x in result[1:])

    # Features with identical positions but different scores
    score_variant_features = [(100, 200, float(i)) for i in range(20)]
    result = remove_covered_features(score_variant_features)
    # Only highest score should remain
    assert result[-1] is True  # Highest score (19.0)
    assert all(not x for x in result[:-1])

    # Linear chain of overlapping features
    chain_features = [(i * 10, i * 10 + 50, float(i)) for i in range(20)]
    result = remove_covered_features(chain_features)
    # All should remain as none completely covers another
    assert all(result), "Chain features incorrectly marked as covered"


def test_remove_covered_features_deterministic() -> None:
    """Test that results are deterministic for the same input."""
    features = generate_random_features(50, max_pos=1000)

    # Run multiple times and ensure same result
    result1 = remove_covered_features(features)
    result2 = remove_covered_features(features)
    result3 = remove_covered_features(features)

    assert result1 == result2 == result3, "Results not deterministic"


def test_optimized_vs_original_correctness() -> None:
    """Test that optimized version produces identical results to original."""
    test_cases = [
        # All original test cases
        [(1, 10, 1.0), (11, 20, 2.0), (21, 30, 3.0)],
        [(1, 10, 1.0), (1, 15, 2.0), (20, 30, 1.0)],
        [(1, 10, 1.0), (1, 15, 2.0), (1, 20, 3.0)],
        [(1, 20, 3.0), (5, 20, 2.0), (10, 20, 1.0)],
        [(1, 20, 1.0), (5, 20, 2.0), (10, 20, 3.0)],
        [(1, 10, 1.0), (1, 10, 2.0), (1, 10, 3.0)],
        [(1, 10, 1.0), (1, 10, 1.0), (1, 10, 1.0)],
        # Random test cases
        generate_random_features(25),
        generate_random_features(75),
        generate_random_features(150),
        # Nested features
        generate_nested_features(20),
        generate_nested_features(60),
    ]

    for i, features in enumerate(test_cases):
        original_result = remove_covered_features(features)
        optimized_result = remove_covered_features_optimized(features)

        assert (
            original_result == optimized_result
        ), f"Results differ for test case {i}: original={original_result}, optimized={optimized_result}"


def test_optimized_performance() -> None:
    """Test that optimized version is faster for larger datasets."""
    sizes = [300, 500, 1000, 2000]

    for size in sizes:
        features = generate_random_features(size)

        # Time original algorithm
        start_time = time.time()
        original_result = remove_covered_features(features)
        original_time = time.time() - start_time

        # Time optimized algorithm
        start_time = time.time()
        optimized_result = remove_covered_features_optimized(features)
        optimized_time = time.time() - start_time

        # Ensure results are identical
        assert original_result == optimized_result, f"Results differ for size {size}"

        # For larger datasets, optimized should be faster
        speedup = original_time / optimized_time if optimized_time > 0 else float("inf")
        print(
            f"Size {size}: Original={original_time:.4f}s, Optimized={optimized_time:.4f}s, Speedup={speedup:.1f}x"
        )

        # Optimized should be at least as fast for large datasets
        if size >= 500:
            assert (
                optimized_time <= original_time * 1.2
            ), f"Optimized version significantly slower for size {size}"


def test_optimized_edge_cases() -> None:
    """Test optimized version with edge cases."""

    # Empty list
    assert remove_covered_features_optimized([]) == []

    # Single feature
    assert remove_covered_features_optimized([(1, 10, 1.0)]) == [True]

    # Small datasets should use original algorithm
    small_features = [(1, 10, 1.0), (5, 15, 2.0)]
    original_result = remove_covered_features(small_features)
    optimized_result = remove_covered_features_optimized(small_features)
    assert original_result == optimized_result
