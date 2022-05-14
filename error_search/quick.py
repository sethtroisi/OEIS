#!/usr/bin/env python

"""
Search for sequences with all prime/composite EXCEPT for a single term (in the middle)
"""

import pprint
import gzip
import gmpy2
import itertools
import os.path

IGNORE = {
    # First pass with all but one index are prime
    'A029911', 'A029913', 'A035209', 'A039786', 'A043298', 'A046414',
    'A047811', 'A054799', 'A054964', 'A055114', 'A056955', 'A059814',
    'A064699', 'A068871', 'A070159', 'A073319', 'A074647', 'A076033',
    'A078447', 'A078860', 'A082595', 'A082700', 'A082843', 'A082869',
    'A089472', 'A095841', 'A099194', 'A106317', 'A107290', 'A109853',
    'A110074', 'A113910', 'A116974', 'A121837', 'A123976', 'A124112',
    'A125635', 'A126788', 'A133323', 'A137189', 'A138000', 'A141250',
    'A141414', 'A143386', 'A144671', 'A145985', 'A160227', 'A167236',
    'A167806', 'A175189', 'A175524', 'A177378', 'A179538', 'A181659',
    'A185720', 'A186884', 'A191611', 'A192141', 'A192229', 'A195258',
    'A200325', 'A203304', 'A208494', 'A209930', 'A212605', 'A212912',
    'A214032', 'A214788', 'A214791', 'A214792', 'A214794', 'A216059',
    'A216538', 'A216827', 'A217789', 'A226977', 'A227627', 'A227757',
    'A228811', 'A228851', 'A236213', 'A237579', 'A239058', 'A239278',
    'A240102', 'A240103', 'A240915', 'A240960', 'A248526', 'A250291',
    'A250292', 'A254713', 'A255420', 'A256170', 'A257110', 'A257771',
    'A259255', 'A260408', 'A260940', 'A262086', 'A263874', 'A263925',
    'A265653', 'A269983', 'A270449', 'A270834', 'A272069', 'A272649',
    'A273064', 'A273377', 'A273460', 'A275274', 'A276435', 'A276495',
    'A279192', 'A280945', 'A292510', 'A293712', 'A294717', 'A295425',
    'A303604', 'A306661', 'A307512', 'A320641', 'A321866', 'A322394',
    'A324256', 'A326715', 'A327914', 'A327915', 'A330225', 'A331948',
    'A334814', 'A335883', 'A335976', 'A336144', 'A340418', 'A342566',
    'A343829', 'A348883', 'A350402', 'A351398',

    # 2nd pass with all but one index prime
    'A003655', 'A015915', 'A015919', 'A021016', 'A029910', 'A060566',
    'A064555', 'A077199', 'A078859', 'A117081', 'A138715', 'A161896',
    'A172514', 'A185656', 'A260406', 'A268465', 'A279124', 'A283354',
    'A292996', 'A307474', 'A340281',

    # 3rd pass with last index decreasing
    'A009856', 'A010957', 'A010958', 'A012610', 'A025590', 'A027534',
    'A027535', 'A036110', 'A036510', 'A036535', 'A060303', 'A072005',
    'A093235', 'A093456', 'A105340', 'A114831', 'A122373', 'A127352',
    'A131316', 'A136154', 'A170921', 'A191829', 'A193546', 'A195603',
    'A198445', 'A215524', 'A246637', 'A280036', 'A285102', 'A288303',
    'A294826', 'A294832', 'A322760', 'A335629', 'A339313', 'A348830',
    'A349599',

    # 4th pass with last index decreasing > 160 characters
    'A003280', 'A003298', 'A011574', 'A011575', 'A011576', 'A011577',
    'A012827', 'A014364', 'A017703', 'A017704', 'A049025', 'A049026',
    'A063317', 'A068846', 'A068918', 'A069053', 'A072240', 'A072323',
    'A073545', 'A076090', 'A080326', 'A087331', 'A087332', 'A090208',
    'A090754', 'A090755', 'A094276', 'A103716', 'A103717', 'A112909',
    'A114095', 'A115184', 'A123257', 'A156339', 'A161117', 'A174782',
    'A291948', 'A292013', 'A310194', 'A311093', 'A311163', 'A311403',
    'A311774', 'A311823', 'A311858', 'A311867', 'A312016', 'A312057',
    'A312141', 'A312528', 'A313030', 'A313189', 'A313196', 'A313600',
    'A314406', 'A315881', 'A349599',

    # Looking at more non-prime indexes
    'A099954', 'A126225', 'A161427', 'A177892', 'A209296', 'A237359',

}


def split_on_condition(seq, condition):
    l1, l2 = itertools.tee((condition(item), item) for item in seq)
    return [i for p, i in l1 if p], [i for p, i in l2 if not p]


def test(seq, terms):
    if len(terms) < 6:
        return False

    def all_but_one(terms, condition):
        """Check if all but one term in a match, return b of single item"""
        B1, B2 = split_on_condition(terms, condition)
        return (len(B1) == 1, B1) if len(B1) <= len(B2) else (len(B2) == 1, B2)

    is_same_primality, not_same_prime = all_but_one(terms[1:], gmpy2.is_prime)
    # Skip first term which can be weird (e.g. 0,1,2...)
    not_same_prime = sorted(set(not_same_prime) - set(range(-1,5)))
    is_same_primality = is_same_primality and len(terms) >= 10 and len(not_same_prime) > 0
    if is_same_primality and gmpy2.is_prime(min(not_same_prime)):
        # Manually inspected > 300 sequences with the error as the last term and found no obvious mistakes.
        # Can ignore only one prime term sequences by uncommenting here
        is_same_primality = False

    is_non_decreasing = [b >= a for a, b in zip(terms, terms[1:])]
    is_last_only_decreasing = all(is_non_decreasing[:-1]) and not is_non_decreasing[-1]
    sum_length = sum(len(str(a)) + 1 for a in terms)
    maybe_truncated = sum_length > 200 and is_last_only_decreasing

    for cond, wrong in [
        (is_same_primality, not_same_prime),
        (maybe_truncated, terms[-2:]),
    ]:
        if cond:
            url = f"https://oeis.org/{seq}"
            if False:
                print(url)
            else:
                print(f"Inspect {seq}: https://oeis.org/{seq}")
                print("\t:", terms)
                print("\t:", wrong)
                print("*", url, "->", ', '.join(map(str, wrong)), "seem wrong")
                print()
            return True

    return False



def parse():
    temp = []

    # TODO load name file as well
    with gzip.open("stripped.gz") as f:
        for line in f:
            line = line.decode()
            if line.startswith("#"):
                continue
            line = line.strip("\n, ").split(",")
            seq, *terms = line
            seq = seq.strip()
            if seq in IGNORE:
                continue

            terms = list(map(int, terms))
            assert terms, (seq, line)
            if test(seq, terms):
                temp.append(seq)

    print(f"Found {len(temp)} Sequences to inspect")
    pprint.pprint(temp, indent=4, compact=True, width=75)


if __name__ == "__main__":
    parse()
