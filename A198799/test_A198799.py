import sympy
import unittest

import A198799

class TestA198799(unittest.TestCase):
    def test_count_circle_with_brute(self):
        for n in range(1, 10000, 13):
            assert A198799.count_brute(n) == A198799.count_circle(n)

    def test_count_circle(self):
        assert A198799.count_circle(49) == 2
        assert A198799.count_circle(13 * 13) == 2
        assert A198799.count_circle(637) == 3
        assert A198799.count_circle(1729) == 4
        assert A198799.count_circle(8281) == 5
        assert A198799.count_circle(12103) == 6

        # test adding factors of 3,5^2,11^2,17^2,23^2 don't change

        assert A198799.count_circle(49 * 3) == 2
        assert A198799.count_circle(49 * 3 ** 2) == 2
        assert A198799.count_circle(49 * 3 ** 3) == 2

        assert A198799.count_circle(49 * 5) == 0
        assert A198799.count_circle(49 * 5 ** 2) == 2
        assert A198799.count_circle(49 * 5 ** 3) == 0
        assert A198799.count_circle(49 * 5 ** 4) == 2

        assert A198799.count_circle(49 * 3 * 5 ** 2 * 23 ** 2) == 2


    def test_gen_small(self):
        primes = [2, 3, 5, 7, 11, 13]
        smallest = A198799.gen_small((1,1,1), count=10, primes=primes)
        for small in smallest:
            assert sorted(sympy.factorint(small).values()) == [1,1,1], small

if __name__ == '__main__':
    unittest.main()
