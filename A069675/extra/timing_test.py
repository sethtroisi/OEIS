import time
import gmpy2

def timeit(a, d, b):
   n = gmpy2.mpz(a * 10 ** d + b)
   is_simple = False
   for i in range(2, 100000):
     if n % i == 0:
       is_simple = True
       break

   t0 = time.time()
   res = gmpy2.is_prime(n)
   t1 = time.time()
   print (res, is_simple, t1 - t0)


timeit(8, 100, 9)   # 0
timeit(6, 855, 7)   # 0.29
timeit(9, 2914, 7)  # 7
timeit(3, 4992, 7)  # 27
timeit(6, 20812, 1) # 1360
       
