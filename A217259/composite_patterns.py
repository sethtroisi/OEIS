import math
import sympy

# Search for another pattern like in A063680

"""
let t = (p^n-d) / (p-1)
and t be a prime

with p = 3
sigma((p-1) * t) = sigma(2) * sigma(t) = 3 * t

with d = 7
sigma((p-1) * t + d) = sigma(p^n) = (p^(n+1) - 1) / (p-1)

After a bunch of complicated expansion
(3^(n+1) - 1)/(3-1) - d  =  (1 + (p-1)) * (1 + t)
(3^(n+1) - 1)/(3-1) - d  =  p * t + p
(3^(n+1) - 1)/(3-1) =  p * t + p + d
(3^(n+1) - 1)/(3-1) = p*(p^n-d)/(p-1) + p + d
(3^(n+1) - 1)/(3-1) = p*(p^n+1-(d+1))/(p-1) + p + d
...
0 = 0!

When I try with other setups (p=5, d=21) (p=7, d=43)
This fails because need DivisorSum[p-1] == p
which only happens for p = 3
"""

COMPOSITES = {
        2: [434, 8575, 8825],
        4: [305635357],
        6: [#104, 147, 596, 1415,
            4850, 5337, 370047, 1630622, 35020303],
            #120221396, 3954451796, 742514284703],
        7: [74, 531434, 387420482],
        8: [#27,1615,1885,
            218984,4218475,312016315,746314601,1125845307],
            #    1132343549,1296114929,9016730984,
            #    303419868239,1197056419121,2065971192041,
            #    2948269852109,4562970154601
        10: [195556, 1152136225],
        12: [   #65,170,209,1394,3393,4407,4556,
                11009,13736,27674,
                38009,38845,47402,76994,157994,162393,184740,
                186686,209294,680609,825359,954521,1243574],
             #   2205209,3515609,4347209,5968502,6539102,6916241,
             #   8165294,10352294,10595009,10786814
}


for D, composites in COMPOSITES.items():
    print(f"Searching {D=} {composites = }")
    for m in range(1, 9):
        if m == 0: continue
        for d in range(-16, 17):
            if m > 1 and math.gcd(m, d) == m:
                # Turns into m * (x + d/m)
                continue

            # (3^n-7)/2 is prime -> 3^n - 7 is a solution
            # solution + 7 -> 3^n

            transformed = [m * x + d for x in composites if m * x + d > 2]
            factored = list(map(sympy.factorint, transformed))
            max_exp = [max(factors.values()) for factors in factored]
            num_factors = [len(factors) for factors in factored]

            # This finds the rule for D=7
            if max(max_exp) >= 8 and min(max_exp) > 2:
                print("\t", m, d, "BIG EXP", factored)

            if D == 4:
                # Only 1 term so hard to generalize
                continue

            # Would hit a lot of primes so only makes sense for 6,8,12
            if max(num_factors) == 1 and len(composites) >= 4:
                print("\t\t", m, d, "RULE 2", factored)

            all_factors = set(f for factors in factored for f in factors.keys())
            common_factor = [k for k in all_factors if all([k in factors for factors in factored])
                    and m % k != 0]
            if common_factor and max(common_factor) > 3:
                print("\t\t", m, d, "COMMON FACTOR", common_factor, factored)
