import gmpy2


a = 2783000011521200439604138286698962565968882984749634760
b = 6302545829970210767062208359546775789414210954519059753

strA = str(a)
strB = str(b)
assert len(strA) == len(strB)

t = strA + strB
assert gmpy2.is_prime(int(t))


strT = str(t)
for i in range(1, len(strA)+1):
    q = int(strA[:i] + strB[::-1][:i][::-1])

    pairs = [int(strA[i] + strB[::-1][i]) for i in range(i)]
    print ()
    print (len(pairs), len(set(pairs)), pairs)
    # sorted(set(range(0,100)) - set(pairs)))

    print (i, len(str(q)), q)
    assert len(pairs) == len(set(pairs))
    assert gmpy2.is_prime(q)
