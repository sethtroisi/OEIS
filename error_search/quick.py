"""
Search for sequences with all prime/composite EXCEPT for a single term (in the middle)
"""
import pprint
import gzip
import gmpy2


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

}


def test(seq, terms):
    # TODO might be nice to add cache infront of terms to handle repeat terms
    is_prime = list(map(gmpy2.is_prime, terms[1:]))

    not_prime = [t for t, is_p in zip(terms[1:], is_prime) if not is_p]
    is_mostly_prime = len(terms) > 10 and is_prime.count(False) == 1 and (set(not_prime) - {-1, 0, 1, 4})

#    not_composite = [t for t, is_p in zip(terms[1:], is_prime) if is_p]
#    is_mostly_composite = len(terms) > 10 and is_prime.count(True) == 1 and (set(not_composite) - {-2, 2, 3, 5})

    for cond, wrong in [
        (is_mostly_prime, not_prime),
#        (is_mostly_composite, not_composite)
    ]:
        if cond:
            url = f"https://oeis.org/{seq}"
            if False:
                print(url)
            else:
                print(f"Inspect {seq}: https://oeis.org/{seq}")
                print("\t:", terms)
                print("\t:", is_prime)
                print("\t:", wrong)
                print("*", url, "->", ', '.join(map(str, wrong)), "seem wrong")
                print()
            return True

    return False



def parse():
    temp = []

    # See https://davidbieber.com/snippets/2020-06-28-oeis-download/
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

    pprint.pprint(temp, indent=4)


if __name__ == "__main__":
    parse()
