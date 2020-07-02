import sympy
import unittest

import A198799

class TestA198799(unittest.TestCase):
    def setUp(self):
        self.valid_primes = [p for p in sympy.primerange(2, 400) if p % 6 == 1]

    def n_to_sig(self, n):
        sig = []
        factored = sympy.factorint(n).items()
        for p_i, e_i in factored:
            if p_i == 3: continue
            if p_i % 6 == 1:
                sig.append(e_i)
            else:
                if e_i % 2 == 1:
                    assert False, factored
                    # Maybe return None
        return sig


    def test_count_circle_with_brute(self):
        for n in range(2, 10000, 13):
            assert A198799.count_brute(n) == A198799.count_circle(n)

    def test_count_circle(self):
        assert A198799.count_circle(49) == 2
        assert A198799.count_circle(13 * 13) == 2
        assert A198799.count_circle(637) == 3
        assert A198799.count_circle(1729) == 4
        assert A198799.count_circle(8281) == 5
        assert A198799.count_circle(12103) == 6

        # test adding any number of 3 doesn't change
        assert A198799.count_circle(49 * 3) == 2
        assert A198799.count_circle(49 * 3 ** 2) == 2
        assert A198799.count_circle(49 * 3 ** 3) == 2

        # test adding even counts of 5,11,17,23 doesn't change
        assert A198799.count_circle(49 * 5) == 0
        assert A198799.count_circle(49 * 5 ** 2) == 2
        assert A198799.count_circle(49 * 5 ** 3) == 0
        assert A198799.count_circle(49 * 5 ** 4) == 2

        assert A198799.count_circle(49 * 3 * 5 ** 2 * 23 ** 2) == 2

    def test_gen_all_signatures_small(self):
        results = list(A198799.gen_all_signatures(2000, self.valid_primes))
        assert results == [
            ((1,), 7),
            ((2,), 7 ** 2),
            ((1,1), 7 * 13),
            ((3,), 7 ** 3),
            ((2,1), 7 ** 2 * 13),
            ((1,1,1), 7 * 13 * 19),
        ], results

    def test_gen_all_signatures(self):
        N = 10000
        expected = {}
        for n in range(2, N):
            factored = sympy.factorint(n)
            if all(prime in self.valid_primes for prime in factored):
                sig = tuple(sorted(factored.values(), reverse=True))
                if sig not in expected:
                    expected[sig] = n
        # Python dict is stable by insertion order
        expected = list(expected.items())

        results = list(A198799.gen_all_signatures(N, self.valid_primes))
        assert results == expected, (results, expected)

    def test_count_signature_with_brute(self):
        for sig, small in A198799.gen_all_signatures(1000, self.valid_primes):
            count_brute = A198799.count_brute(small)
            count_sig = A198799.count_signature(sig)
            assert count_brute == count_sig, (small, sig)

    def test_count_signature_with_circle(self):
        for sig, small in A198799.gen_all_signatures(100000, self.valid_primes):
            count_circle = A198799.count_circle(small)
            count_sig = A198799.count_signature(sig)
            assert count_circle == count_sig, (small, sig)

    def test_gen_small(self):
        primes = [2, 3, 5, 7, 11, 13]
        smallest = A198799.gen_small((1,1,1), count=10, primes=primes)
        assert len(smallest) == 10
        for small in smallest:
            assert sorted(sympy.factorint(small).values()) == [1,1,1], small

    def test_A198774(self):
        # All of these should return m=3
        n = [
            637, 931, 1183, 1519, 1813, 1911, 2107, 2401, 2527, 2548, 2793,
            9751, 18228, 35017, 69727, 105889, 142492, 179977,
        ]

        for an in n:
            assert A198799.count_circle(an) == 3
            sig = self.n_to_sig(an)
            assert A198799.count_signature(sig) == 3

    def test_A198775(self):
        # All of these should return m=4
        n = [
            1729, 2821, 3367, 3913, 4123, 4459, 4921, 5187, 5551, 5719, 6097,
            22204, 24583, 25753, 27937, 30121, 32116, 34333, 37387, 44548,
        ]

        for an in n:
            assert A198799.count_circle(an) == 4
            sig = self.n_to_sig(an)
            assert A198799.count_signature(sig) == 4

    def test_previous_A198799(self):
        A = [
            49, 637, 1729, 8281, 12103, 1529437, 53599, 157339, 593047,
            19882681, 375193, 68574961, 2989441, 7709611, 1983163,
            47738317081, 4877509, 21169376772835837, 18384457, 377770939,
        ]

        for m, an in enumerate(A, 2):
            if an < 1e12:
                assert A198799.count_circle(an) == m
            sig = self.n_to_sig(an)
            assert A198799.count_signature(sig) == m


if __name__ == '__main__':
    unittest.main()
