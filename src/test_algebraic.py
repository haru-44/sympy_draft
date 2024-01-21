from algebraic_ntheory import class_number_of_quadratic_field, QuadraticFormClassGroup, _class_number_neg_naive

import pytest

def test_quadratic_form_class_group_reduce():
    qf = QuadraticFormClassGroup(-20)
    assert qf.reduction(1, 2, 6) == (1, 0, 5)
    assert qf.reduction(18, -14, 3) == (2, 2, 3)

def test_quadratic_form_class_group_composition():
     qf = QuadraticFormClassGroup(-15) # = Z_2
     assert qf.composition(1, 1, 4, 1, 1, 4) == (1, 1, 4)
     assert qf.composition(2, 1, 2, 1, 1, 4) == (2, 1, 2)
     assert qf.composition(2, 1, 2, 2, 1, 2) == (1, 1, 4)

     qf = QuadraticFormClassGroup(-44) # = Z_3
     assert qf.composition(1, 0, 11, 1, 0, 11) == (1, 0, 11)
     assert qf.composition(1, 0, 11, 3, 2, 4) == (3, 2, 4)
     assert qf.composition(3, 2, 4, 3, 2, 4) == (3, -2, 4)
     assert qf.composition(3, 2, 4, 3, -2, 4) == (1, 0, 11)

def test_quadratic_form_class_group_composition_double():
    qf = QuadraticFormClassGroup(-32) # = Z_2
    assert qf.composition_double(1, 0, 8) == (1, 0, 8)
    assert qf.composition_double(3, 2, 3) == (1, 0, 8)

    qf = QuadraticFormClassGroup(-63) # = Z_4
    assert qf.composition_double(2, 1, 8) == (4, 1, 4)
    assert qf.composition_double(4, 1, 4) == (1, 1, 16)
    assert qf.composition_double(2, -1, 8) == (4, 1, 4)

def test_quadratic_form_class_group_n_times():
    qf = QuadraticFormClassGroup(-83) # Z_3
    assert qf.n_times(3, 1, 7, 2) == (3, -1, 7)
    assert qf.n_times(3, 1, 7, 3) == (1, 1, 21)

def test_quadratic_form_class_group_counting_naive():
    assert QuadraticFormClassGroup(-3).counting_naive() == 1
    assert QuadraticFormClassGroup(-4).counting_naive() == 1
    A014599 = [1, 1, 1, 2, 1, 3, 1, 3, 2, 4, 1, 5, 2, 4,
               3, 4, 1, 7, 2, 5, 3, 6, 2, 8, 2, 5, 3, 8,
               2, 10, 2, 5, 5, 6, 3, 10, 2, 7, 4, 10, 1,
               11, 4, 6, 5, 8, 2, 13, 4, 9, 4, 6, 3, 14,
               4, 7, 5, 12, 2, 15, 3, 6, 7, 12, 4, 13, 2]
    disc = -3
    for h in A014599:
        assert QuadraticFormClassGroup(disc).counting_naive() == h
        disc -= 4
    A000003 = [1, 1, 1, 1, 2, 2, 1, 2, 2, 2, 3, 2, 2, 4,
               2, 2, 4, 2, 3, 4, 4, 2, 3, 4, 2, 6, 3, 2,
               6, 4, 3, 4, 4, 4, 6, 4, 2, 6, 4, 4, 8, 4,
               3, 6, 4, 4, 5, 4, 4, 6, 6, 4, 6, 6, 4, 8,
               4, 2, 9, 4, 6, 8, 4, 4, 8, 8, 3, 8, 8, 4]
    disc = -4
    for h in A000003:
        assert QuadraticFormClassGroup(disc).counting_naive() == h
        disc -= 4

def test_quadratic_form_class_group_counting_order():
    for disc in range(3, 1000, 4):
        qf = QuadraticFormClassGroup(-disc)
        assert qf.counting_order() == qf.counting_naive()

    for disc in range(4, 1000, 4):
        qf = QuadraticFormClassGroup(-disc)
        assert qf.counting_order() == qf.counting_naive()

@pytest.mark.skipif(True, reason="slow")
def test_class_number_of_quadratic_field():
    # imaginary
    ds = [1, 2, 3, 7, 11, 19, 43, 67, 163] # A003173
    assert all(class_number_of_quadratic_field(d=-d) == 1 for d in ds)

    discs_list = [(1, [3, 4, 7, 8, 11, 19, 43, 67, 163]), # A014602
                  (3, [23, 31, 59, 83, 107, 139, 211, 283,
                       307, 331, 379, 499, 547, 643, 883,
                       907]), # A006203
                  (23, [647, 1039, 1103, 1279, 1447, 1471,
                        1811, 1979, 2411, 2671, 3491, 3539,
                        3847, 3923, 4211, 4783, 5387, 5507,
                        5531, 6563, 6659, 6703, 7043, 9587,
                        9931, 10867, 10883, 12203, 12739,
                        13099, 13187, 15307, 15451, 16267,
                        17203, 17851, 18379, 20323]), # A046020
                 ]
    for h, discs in discs_list:
        assert all(class_number_of_quadratic_field(disc=-disc) == h for disc in discs)

    # real
    ds_list = [(1, [2, 3, 5, 6, 7, 11, 13, 14, 17, 19,
                    21, 22, 23, 29, 31, 33, 37]), # A003172
               (2, [10, 15, 26, 30, 34, 35, 39, 42, 51]), # A029702
               (3, [79, 142, 223, 229, 254, 257, 321, 326,
                    359, 443, 469, 473, 659, 733, 761, 839]), # A029703
               (4, [82, 130, 145, 170, 195, 210, 219, 231,
                     255, 274, 290, 291, 322, 323, 330, 370,
                     390, 410, 434, 435]), # A029704
               (5, [401, 439, 499, 727, 817, 982, 1093, 1126,
                    1327, 1393, 1429, 1486, 1641, 1766, 1897,
                    2027, 2081, 2153]), # A029705
               (6, [235, 346, 427, 506, 574, 697, 785, 786,
                    842, 874, 894, 895, 898, 899, 906, 985,
                    1086, 1191, 1211]), # A218038
               (7, [577, 1009, 1087, 1294, 1601, 1761, 1934,
                    2029, 2251, 2302, 2467, 2913, 4139, 4229,
                    4702, 5039, 5273]), # A218039
               (8, [226, 399, 442, 646, 799, 870, 910, 994,
                    1023, 1122, 1155, 1239, 1290, 1299, 1351,
                    1443, 1446, 1590]), # A218040
               (9, [1129, 1654, 3137, 3719, 4409, 4534, 5521,
                    5623, 5878, 6809, 7573, 7873]), # A218041
               (10, [1111, 1226, 2031, 2335, 2362, 2602, 2986,
                     3129, 3246, 3379, 3585, 3598, 3599, 3722,
                     3782, 3966, 4097]) # A218042
               ]
    for h, ds in ds_list:
        assert all(class_number_of_quadratic_field(d=d) == h for d in ds)

