from proksee_batch.remove_covered_features import remove_covered_features


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
    test_list = [(1, 10, 1), (11, 20, 2), (21, 30, 3)]
    expected_result = [True, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result
    test_list = [(1, 10, 3), (11, 20, 2), (21, 30, 1)]
    expected_result = [True, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 2: Overlap
    # 1-10, 1-15, 20-30
    # 1-10 is covered
    # 5-15 is not covered
    # 20-30 is not covered
    # Expected result: [False, True, True]
    test_list = [(1, 10, 1), (1, 15, 2), (20, 30, 1)]
    expected_result = [False, True, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 3: Overlap
    # 1-10, 1-15, 1-20
    # 1-10 is covered
    # 1-15 is covered
    # 1-20 is not covered
    # Expected result: [False, False, True]
    test_list = [(1, 10, 1), (1, 15, 2), (1, 20, 3)]
    expected_result = [False, False, True]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 4: Overlap
    # 1-20, 5-20, 10-20
    # 1-20 is not covered
    # 5-20 is covered
    # 10-20 is covered
    # Expected result: [True, False, False]
    test_list = [(1, 20, 3), (5, 20, 2), (10, 20, 1)]
    expected_result = [True, False, False]
    result = remove_covered_features(test_list)
    assert result == expected_result

    # Test case 5: Overlap
    # 1-20, 5-20, 10-20
    # 1-20 is not covered
    # 5-20 is not covered
    # 10-20 is not covered
    # Expected result: [True, True, True]
    test_list = [(1, 20, 1), (5, 20, 2), (10, 20, 3)]
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
        (1, 20, 10),
        (2, 19, 2),
        (3, 18, 3),
        (4, 17, 4),
        (5, 16, 5),
        (6, 15, 6),
        (7, 14, 7),
        (8, 13, 8),
        (9, 12, 9),
        (10, 11, 11),
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
