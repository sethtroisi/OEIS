import math

def v(x, y):
    return 3 * x * x + 4 * y * y

def count_pairs(N):
    max_x = math.isqrt(N // 3)
    assert v(max_x, 0) <= N < v(max_x + 1, 0);

    y = math.isqrt(N // 4)
    assert v(0, y) <= N < v(0, y + 1)

    # (0, 0) doesn't count
    count = -1
    for x in range(0, max_x+1):
        # reduce y till valid
        while v(x, y) > N:
            y -= 1

        assert y >= 0
        count += (y + 1)

    return count


for bits in range(25, 51):
    count = count_pairs(2 ** bits)
    print(f"| {bits:2d} | ???  | {count} |")
