# flake8: noqa: E501
import unittest
import time
import random
from itertools import combinations
from binascii import hexlify

from bls_py.aggregation_info import AggregationInfo
from bls_py.bls import BLS
from bls_py.ec import (default_ec, default_ec_twist, generator_Fq, generator_Fq2,
                 hash_to_point_Fq, hash_to_point_Fq2, sw_encode, twist, untwist,
                 y_for_x)
from bls_py.fields import Fq, Fq2, Fq6, Fq12, bls12381_q as Q
from bls_py.keys import PrivateKey, PublicKey, ExtendedPrivateKey
from bls_py.signature import Signature
from bls_py.threshold import Threshold
from bls_py.util import hash256
from bls_py.tdata import (fq_list, fq_res_list, fq2_list, fq2_res_list,
                          fq6_list, fq6_res_list, fq12_list, fq12_res_list)


def rand_scalar(ec=default_ec):
    return random.randrange(1, ec.n)


class TestBLS(unittest.TestCase):

    def test_fields(self):
        a = Fq(Q, 30)
        b = Fq(Q, -18)
        c = Fq2(Q, a, b)
        d = Fq2(Q, a + a, -5)
        e = c * d
        f = e * d
        assert(f != e)
        e_sq = e*e
        e_sqrt = e_sq.modsqrt()
        assert(pow(e_sqrt, 2) == e_sq)

        a2 = Fq(Q,
                3012492130751239573498573249085723940848571098237509182375)
        b2 = Fq(Q, 3432984572394572309458723045723849)
        c2 = Fq2(Q, a2, b2)
        assert(b2 != c2)

        g = Fq6(Q, c, d, d*d*c)
        h = Fq6(Q, a + a*c, c*b*a, b*b*d*21)
        i = Fq12(Q, g, h)
        assert(~(~i) == i)
        assert((~(i.root)) * i.root == Fq6.one(Q))
        x = Fq12(Q, Fq6.zero(Q), i.root)
        assert((~x) * x == Fq12.one(Q))

        j = Fq6(Q, a + a*c, Fq2.zero(Q), Fq2.zero(Q))
        j2 = Fq6(Q, a + a*c, Fq2.zero(Q), Fq2.one(Q))
        assert(j == (a + a*c))
        assert(j2 != (a + a*c))
        assert(j != j2)

        # Test frob_coeffs
        one = Fq(default_ec.q, 1)
        two = one + one
        a = Fq2(default_ec.q, two, two)
        b = Fq6(default_ec.q, a, a, a)
        c = Fq12(default_ec.q, b, b)
        for base in (a, b, c):
            for expo in range(1, base.extension):
                assert base.qi_power(expo) == pow(base, pow(default_ec.q, expo))


    def test_ec(self):
        g = generator_Fq(default_ec)

        assert(g.is_on_curve())
        assert(2*g == g + g)
        assert((3*g).is_on_curve())
        assert(3*g == g + g + g)
        P = hash_to_point_Fq(bytes([]))
        assert(P.is_on_curve())
        assert(P.serialize() == bytes.fromhex("12fc5ad5a2fbe9d4b6eb0bc16d530e5f263b6d59cbaf26c3f2831962924aa588ab84d46cc80d3a433ce064adb307f256"))

        g2 = generator_Fq2(default_ec_twist)
        assert(g2.x * (2 * g2.y) == 2*(g2.x * g2.y))
        assert(g2.is_on_curve())
        s = g2 + g2
        assert(untwist(twist(s)) == s)
        assert(untwist(5 * twist(s)) == 5 * s)
        assert(5 * twist(s) == twist(5 * s))
        assert(s.is_on_curve())
        assert(g2.is_on_curve())
        assert(g2 + g2 == 2 * g2)
        assert(g2 * 5 == (g2 * 2) + (2 * g2) + g2)
        y = y_for_x(g2.x, default_ec_twist, Fq2)
        assert(y[0] == g2.y or y[1] == g2.y)
        assert(hash_to_point_Fq2("chia") == hash_to_point_Fq2("chia"))

        g_j = generator_Fq(default_ec_twist).to_jacobian()
        g2_j = generator_Fq2(default_ec_twist).to_jacobian()
        g2_j2 = (generator_Fq2(default_ec_twist) * 2).to_jacobian()
        assert(g.to_jacobian().to_affine() == g)
        assert((g_j * 2).to_affine() == g * 2)
        assert((g2_j + g2_j2).to_affine() == g2 * 3)

        assert(sw_encode(Fq(default_ec.q, 0)).infinity)
        assert(sw_encode(Fq(default_ec.q, 1)) == sw_encode(Fq(default_ec.q, -1)).negate())
        assert(sw_encode(Fq(default_ec.q, 0x019cfaba0c258165d092f6bca9a081871e62a126c499340dc71c0e9527f923f3b299592a7a9503066cc5362484d96dd7)) == generator_Fq())
        assert(sw_encode(Fq(default_ec.q, 0x186417302d5a65347a88b0f999ab2b504614aa5e2eebdeb1a014c40bceb7d2306c12a6d436befcf94d39c9db7b263cd4)) == generator_Fq().negate())


    def test_vectors(self):
        sk1 = PrivateKey.from_seed(bytes([1, 2, 3, 4, 5]))
        pk1 = sk1.get_public_key()
        sig1 = sk1.sign(bytes([7, 8, 9]))

        sk2 = PrivateKey.from_seed(bytes([1, 2, 3, 4, 5, 6]))
        pk2 = sk2.get_public_key()
        sig2 = sk2.sign(bytes([7, 8, 9]))
        assert(sk1.serialize() == bytes.fromhex("022fb42c08c12de3a6af053880199806532e79515f94e83461612101f9412f9e"))
        assert(pk1.get_fingerprint() == 0x26d53247)
        assert(pk2.get_fingerprint() == 0x289bb56e)
        assert(sig1.serialize() == bytes.fromhex("93eb2e1cb5efcfb31f2c08b235e8203a67265bc6a13d9f0ab77727293b74a357ff0459ac210dc851fcb8a60cb7d393a419915cfcf83908ddbeac32039aaa3e8fea82efcb3ba4f740f20c76df5e97109b57370ae32d9b70d256a98942e5806065"))
        assert(sig2.serialize() == bytes.fromhex("975b5daa64b915be19b5ac6d47bc1c2fc832d2fb8ca3e95c4805d8216f95cf2bdbb36cc23645f52040e381550727db420b523b57d494959e0e8c0c6060c46cf173872897f14d43b2ac2aec52fc7b46c02c5699ff7a10beba24d3ced4e89c821e"))

        agg_sig = BLS.aggregate_sigs([sig1, sig2])
        agg_pk = BLS.aggregate_pub_keys([pk1, pk2], True)
        agg_sk = BLS.aggregate_priv_keys([sk1, sk2], [pk1, pk2], True)
        assert(agg_sig.serialize() == bytes.fromhex("0a638495c1403b25be391ed44c0ab013390026b5892c796a85ede46310ff7d0e0671f86ebe0e8f56bee80f28eb6d999c0a418c5fc52debac8fc338784cd32b76338d629dc2b4045a5833a357809795ef55ee3e9bee532edfc1d9c443bf5bc658"))
        assert(agg_sk.sign(bytes([7, 8, 9])).serialize() == agg_sig.serialize())


        assert(BLS.verify(sig1))
        assert(BLS.verify(agg_sig))

        agg_sig.set_aggregation_info(AggregationInfo.from_msg(agg_pk, bytes([7, 8, 9])))
        assert(BLS.verify(agg_sig))

        sig1.set_aggregation_info(sig2.aggregation_info)
        assert(not BLS.verify(sig1))

        sig3 = sk1.sign(bytes([1, 2, 3]))
        sig4 = sk1.sign(bytes([1, 2, 3, 4]))
        sig5 = sk2.sign(bytes([1, 2]))


        agg_sig2 = BLS.aggregate_sigs([sig3, sig4, sig5])
        assert(BLS.verify(agg_sig2))
        assert(agg_sig2.serialize() == bytes.fromhex("8b11daf73cd05f2fe27809b74a7b4c65b1bb79cc1066bdf839d96b97e073c1a635d2ec048e0801b4a208118fdbbb63a516bab8755cc8d850862eeaa099540cd83621ff9db97b4ada857ef54c50715486217bd2ecb4517e05ab49380c041e159b"))


    def test_vectors2(self):
        m1 = bytes([1, 2, 3, 40])
        m2 = bytes([5, 6, 70, 201])
        m3 = bytes([9, 10, 11, 12, 13])
        m4 = bytes([15, 63, 244, 92, 0, 1])

        sk1 = PrivateKey.from_seed(bytes([1, 2, 3, 4, 5]))
        sk2 = PrivateKey.from_seed(bytes([1, 2, 3, 4, 5, 6]))

        sig1 = sk1.sign(m1)
        sig2 = sk2.sign(m2)
        sig3 = sk2.sign(m1)
        sig4 = sk1.sign(m3)
        sig5 = sk1.sign(m1)
        sig6 = sk1.sign(m4)

        sig_L = BLS.aggregate_sigs([sig1, sig2])
        sig_R = BLS.aggregate_sigs([sig3, sig4, sig5])
        assert(BLS.verify(sig_L))
        assert(BLS.verify(sig_R))

        sig_final = BLS.aggregate_sigs([sig_L, sig_R, sig6])
        assert(sig_final.serialize() == bytes.fromhex("07969958fbf82e65bd13ba0749990764cac81cf10d923af9fdd2723f1e3910c3fdb874a67f9d511bb7e4920f8c01232b12e2fb5e64a7c2d177a475dab5c3729ca1f580301ccdef809c57a8846890265d195b694fa414a2a3aa55c32837fddd80"))
        assert(BLS.verify(sig_final))
        quotient = sig_final.divide_by([sig2, sig5, sig6])
        assert(BLS.verify(quotient))
        assert(BLS.verify(sig_final))
        assert(quotient.serialize() == bytes.fromhex("8ebc8a73a2291e689ce51769ff87e517be6089fd0627b2ce3cd2f0ee1ce134b39c4da40928954175014e9bbe623d845d0bdba8bfd2a85af9507ddf145579480132b676f027381314d983a63842fcc7bf5c8c088461e3ebb04dcf86b431d6238f"))
        assert(quotient.divide_by([]) == quotient)
        try:
            quotient.divide_by([sig6])
            assert(False)  # Should fail due to not subset
        except:
            pass
        sig_final.divide_by([sig1]) # Should not throw
        try:
            sig_final.divide_by([sig_L]) # Should throw due to not unique
            assert(False)  # Should fail due to not unique
        except:
            pass

        # Divide by aggregate
        sig7 = sk2.sign(m3)
        sig8 = sk2.sign(m4)
        sig_R2 = BLS.aggregate_sigs([sig7, sig8])
        sig_final2 = BLS.aggregate_sigs([sig_final, sig_R2])
        quotient2 = sig_final2.divide_by([sig_R2])
        assert(BLS.verify(quotient2))
        assert(quotient2.serialize() == bytes.fromhex("06af6930bd06838f2e4b00b62911fb290245cce503ccf5bfc2901459897731dd08fc4c56dbde75a11677ccfbfa61ab8b14735fddc66a02b7aeebb54ab9a41488f89f641d83d4515c4dd20dfcf28cbbccb1472c327f0780be3a90c005c58a47d3"))


    def test_vectors3(self):
        seed = bytes([1, 50, 6, 244, 24, 199, 1, 25])
        esk =  ExtendedPrivateKey.from_seed(seed)
        assert(esk.private_key.get_public_key().get_fingerprint() == 0xa4700b27)
        assert(hexlify(esk.chain_code).decode('ascii') == "d8b12555b4cc5578951e4a7c80031e22019cc0dce168b3ed88115311b8feb1e3")
        esk77 = esk.private_child(77 + 2**31)
        assert(hexlify(esk77.chain_code).decode('ascii') == "f2c8e4269bb3e54f8179a5c6976d92ca14c3260dd729981e9d15f53049fd698b")
        assert(esk77.private_key.get_public_key().get_fingerprint() == 0xa8063dcf)

        assert(esk.private_child(3)
                  .private_child(17)
                  .private_key
                  .get_public_key()
                  .get_fingerprint() == 0xff26a31f)

        assert(esk.get_extended_public_key()
                  .public_child(3)
                  .public_child(17)
                  .get_public_key()
                  .get_fingerprint() == 0xff26a31f)


    def test1(self):
        seed = bytes([0, 50, 6, 244, 24, 199, 1, 25, 52, 88, 192,
                      19, 18, 12, 89, 6, 220, 18, 102, 58, 209,
                      82, 12, 62, 89, 110, 182, 9, 44, 20, 254, 22])
        sk = PrivateKey.from_seed(seed)
        pk = sk.get_public_key()

        msg = bytes([100, 2, 254, 88, 90, 45, 23])

        sig = sk.sign(msg)

        sk_bytes = sk.serialize()
        pk_bytes = pk.serialize()
        sig_bytes = sig.serialize()

        sk = PrivateKey.from_bytes(sk_bytes)
        pk = PublicKey.from_bytes(pk_bytes)
        sig = Signature.from_bytes(sig_bytes)

        sig.set_aggregation_info(AggregationInfo.from_msg(pk, msg))
        assert(BLS.verify(sig))

        seed = bytes([1]) + seed[1:]
        sk1 = PrivateKey.from_seed(seed)
        seed = bytes([2]) + seed[1:]
        sk2 = PrivateKey.from_seed(seed)

        pk1 = sk1.get_public_key()
        sig1 = sk1.sign(msg)

        pk2 = sk2.get_public_key()
        sig2 = sk2.sign(msg)

        agg_sig = BLS.aggregate_sigs([sig1, sig2])
        agg_pubkey = BLS.aggregate_pub_keys([pk1, pk2], True)

        agg_sig.set_aggregation_info(AggregationInfo.from_msg(agg_pubkey, msg))
        assert(BLS.verify(agg_sig))

        seed = bytes([3]) + seed[1:]
        sk3 = PrivateKey.from_seed(seed)
        pk3 = sk3.get_public_key()
        msg2 = bytes([100, 2, 254, 88, 90, 45, 23])

        sig1 = sk1.sign(msg)
        sig2 = sk2.sign(msg)
        sig3 = sk3.sign(msg2)
        agg_sig_l = BLS.aggregate_sigs([sig1, sig2])
        agg_sig_final = BLS.aggregate_sigs([agg_sig_l, sig3])

        sig_bytes = agg_sig_final.serialize()

        agg_sig_final = Signature.from_bytes(sig_bytes)
        a1 = AggregationInfo.from_msg(pk1, msg)
        a2 = AggregationInfo.from_msg(pk2, msg)
        a3 = AggregationInfo.from_msg(pk3, msg2)
        a1a2 = AggregationInfo.merge_infos([a1, a2])
        a_final = AggregationInfo.merge_infos([a1a2, a3])
        print(a_final)
        agg_sig_final.set_aggregation_info(a_final)
        assert(BLS.verify(agg_sig_final))

        assert(BLS.verify(agg_sig_l))
        agg_sig_final = agg_sig_final.divide_by([agg_sig_l])

        assert(BLS.verify(agg_sig_final))

        agg_sk = BLS.aggregate_priv_keys([sk1, sk2], [pk1, pk2], True)
        agg_sk.sign(msg)

        seed = bytes([1, 50, 6, 244, 24, 199, 1, 25, 52, 88, 192,
                      19, 18, 12, 89, 6, 220, 18, 102, 58, 209,
                      82, 12, 62, 89, 110, 182, 9, 44, 20, 254, 22])

        esk = ExtendedPrivateKey.from_seed(seed)
        epk = esk.get_extended_public_key()

        sk_child = esk.private_child(0).private_child(5)
        pk_child = epk.public_child(0).public_child(5)

        assert(sk_child.get_extended_public_key() == pk_child)


    def test2(self):
        seed = bytes([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        seed2 = bytes([1, 20, 102, 229, 1, 157])

        sk = PrivateKey.from_seed(seed)
        sk_cp = PrivateKey.from_seed(seed)
        sk2 = PrivateKey.from_seed(seed2)
        pk = sk.get_public_key()
        pk2 = sk2.get_public_key()
        assert(sk == sk_cp)
        assert(sk != sk2)
        assert(pk.get_fingerprint() == 0xddad59bb)

        pk2_ser = pk2.serialize()
        pk2_copy = PublicKey.from_bytes(pk2_ser)
        assert(pk2 == pk2_copy)
        assert(pk != pk2)
        assert(pk2.size() == 48)
        assert(sk2.size() == 32)

        message = bytes("this is the message", "utf-8")
        sig = sk.sign(message)
        sig_ser = sig.serialize()
        sig_cp = Signature.from_bytes(sig_ser)
        a1 = AggregationInfo.from_msg(pk, message)
        sig_cp.set_aggregation_info(a1)
        a2 = sig_cp.get_aggregation_info()
        assert(a1 == a2)
        sig2 = sk2.sign(message)

        assert(sig.size() == 96)
        assert(sig != sig2)
        assert(sig == sig_cp)

        sig_agg = BLS.aggregate_sigs([sig, sig2])

        result = BLS.verify(sig_cp)
        result2 = BLS.verify(sig2)
        result3 = BLS.verify(sig_agg)
        assert(result)
        assert(result2)
        assert(result3)

        sk2 = sk


    def _test_threshold_instance(self, T, N):
        commitments = []
        # fragments[i][j] = fragment held by player i,
        #                   received from player j
        fragments = [[None] * N for _ in range(N)]
        secrets = []

        # Step 1 : PrivateKey.new_threshold
        for player in range(N):
            secret_key, commi, frags = PrivateKey.new_threshold(T, N)
            for target, frag in enumerate(frags):
                fragments[target][player] = frag
            commitments.append(commi)
            secrets.append(secret_key)

        # Step 2 : Threshold.verify_secret_fragment
        for player_source in range(1, N+1):
            for player_target in range(1, N+1):
                assert Threshold.verify_secret_fragment(
                    T, fragments[player_target - 1][player_source - 1],
                    player_target, commitments[player_source - 1])

        # Step 3 : master_pubkey = BLS.aggregate_pub_keys(...)
        #          secret_share = BLS.aggregate_priv_keys(...)
        master_pubkey = BLS.aggregate_pub_keys(
               [PublicKey.from_g1(cpoly[0].to_jacobian())
                for cpoly in commitments],
               False)

        secret_shares = [BLS.aggregate_priv_keys(map(PrivateKey, row), None, False)
                         for row in fragments]

        master_privkey = BLS.aggregate_priv_keys(secrets, None, False)
        msg = 'Test'
        signature_actual = master_privkey.sign(msg)

        # Step 4 : sig_share = secret_share.sign_threshold(...)
        # Check every combination of T players
        for X in combinations(range(1, N+1), T):
            # X: a list of T indices like [1, 2, 5]

            # Check underlying secret key is correct
            r = Threshold.interpolate_at_zero(X,
                    [secret_shares[x-1].value for x in X])
            secret_cand = PrivateKey(r)
            assert secret_cand == master_privkey

            # Check signatures
            signature_shares = [secret_shares[x-1].sign_threshold(msg, x, X)
                                for x in X]
            signature_cand = BLS.aggregate_sigs_simple(signature_shares)
            assert signature_cand == signature_actual

        # Check that the signature actually verifies the message
        agg_info = AggregationInfo.from_msg(master_pubkey, msg)
        signature_actual.set_aggregation_info(agg_info)
        assert BLS.verify(signature_actual)

        # Step 4b : Alternatively, we can add the lagrange coefficients
        # to 'unit' signatures.
        for X in combinations(range(1, N+1), T):
            # X: a list of T indices like [1, 2, 5]

            # Check signatures
            signature_shares = [secret_shares[x-1].sign(msg) for x in X]
            signature_cand = Threshold.aggregate_unit_sigs(signature_shares, X, T)
            assert signature_cand == signature_actual


    def test_threshold(self):
        self._test_threshold_instance(1, 1)
        self._test_threshold_instance(1, 2)
        self._test_threshold_instance(2, 2)
        for T in range(1, 6):
            self._test_threshold_instance(T, 5)


class TestFields(unittest.TestCase):

    def print_fqx(self, fqx):
        print(f'    {fqx},')

    def test_fq(self):
        a = int(fq_list[0]); b = int(fq_list[1]);
        c = int(fq_list[2]); d = int(fq_list[3]);
        fqa = fq_list[0]; fqb = fq_list[1];
        fqc = fq_list[2]; fqd = fq_list[3];

        # fq op fq
        assert fqa+fqb == fq_res_list[0]
        assert fqa+fqc == fq_res_list[1]
        assert fqa+fqd == fq_res_list[2]
        assert fqb+fqc == fq_res_list[3]
        assert fqb+fqd == fq_res_list[4]
        assert fqc+fqd == fq_res_list[5]

        assert fqa*fqb == fq_res_list[6]
        assert fqa*fqc == fq_res_list[7]
        assert fqa*fqd == fq_res_list[8]
        assert fqb*fqc == fq_res_list[9]
        assert fqb*fqd == fq_res_list[10]
        assert fqc*fqd == fq_res_list[11]

        assert fqa-fqb == fq_res_list[12]
        assert fqa-fqc == fq_res_list[13]
        assert fqa-fqd == fq_res_list[14]
        assert fqb-fqc == fq_res_list[15]
        assert fqb-fqd == fq_res_list[16]
        assert fqc-fqd == fq_res_list[17]

        # int op fq
        assert a+fqb == fq_res_list[0]
        assert a+fqc == fq_res_list[1]
        assert a+fqd == fq_res_list[2]
        assert b+fqc == fq_res_list[3]
        assert b+fqd == fq_res_list[4]
        assert c+fqd == fq_res_list[5]

        assert a*fqb == fq_res_list[6]
        assert a*fqc == fq_res_list[7]
        assert a*fqd == fq_res_list[8]
        assert b*fqc == fq_res_list[9]
        assert b*fqd == fq_res_list[10]
        assert c*fqd == fq_res_list[11]

        assert a-fqb == fq_res_list[12]
        assert a-fqc == fq_res_list[13]
        assert a-fqd == fq_res_list[14]
        assert b-fqc == fq_res_list[15]
        assert b-fqd == fq_res_list[16]
        assert c-fqd == fq_res_list[17]

        # fq op int
        assert fqa+b == fq_res_list[0]
        assert fqa+c == fq_res_list[1]
        assert fqa+d == fq_res_list[2]
        assert fqb+c == fq_res_list[3]
        assert fqb+d == fq_res_list[4]
        assert fqc+d == fq_res_list[5]

        assert fqa*b == fq_res_list[6]
        assert fqa*c == fq_res_list[7]
        assert fqa*d == fq_res_list[8]
        assert fqb*c == fq_res_list[9]
        assert fqb*d == fq_res_list[10]
        assert fqc*d == fq_res_list[11]

        assert fqa-b == fq_res_list[12]
        assert fqa-c == fq_res_list[13]
        assert fqa-d == fq_res_list[14]
        assert fqb-c == fq_res_list[15]
        assert fqb-d == fq_res_list[16]
        assert fqc-d == fq_res_list[17]

        # op fq
        assert -fqa == fq_res_list[18]
        assert -(-fqa) == fqa
        assert -fqb == fq_res_list[19]
        assert -(-fqb) == fqb
        assert ~fqa == fq_res_list[20]
        assert ~(~fqa) == fqa
        assert ~fqb == fq_res_list[21]
        assert ~(~fqb) == fqb

    def test_fq2(self):
        a = int(fq_list[0]); b = int(fq_list[1]);
        c = int(fq_list[2]); d = int(fq_list[3]);
        fqa = fq_list[0]; fqb = fq_list[1];
        fqc = fq_list[2]; fqd = fq_list[3];
        fq2a = fq2_list[0]; fq2b = fq2_list[1]
        fq2c = fq2_list[2]; fq2d = fq2_list[3];

        # fq2 op fq2
        assert fq2a+fq2b == fq2_res_list[0]
        assert fq2a+fq2c == fq2_res_list[1]
        assert fq2a+fq2d == fq2_res_list[2]
        assert fq2b+fq2c == fq2_res_list[3]
        assert fq2b+fq2d == fq2_res_list[4]
        assert fq2c+fq2d == fq2_res_list[5]

        assert fq2a*fq2b == fq2_res_list[6]
        assert fq2a*fq2c == fq2_res_list[7]
        assert fq2a*fq2d == fq2_res_list[8]
        assert fq2b*fq2c == fq2_res_list[9]
        assert fq2b*fq2d == fq2_res_list[10]
        assert fq2c*fq2d == fq2_res_list[11]

        assert fq2a-fq2b == fq2_res_list[12]
        assert fq2a-fq2c == fq2_res_list[13]
        assert fq2a-fq2d == fq2_res_list[14]
        assert fq2b-fq2c == fq2_res_list[15]
        assert fq2b-fq2d == fq2_res_list[16]
        assert fq2c-fq2d == fq2_res_list[17]

        # fq op fq2
        assert fqa+fq2b == fq2_res_list[18]
        assert fqa+fq2c == fq2_res_list[19]
        assert fqa+fq2d == fq2_res_list[20]
        assert fqb+fq2c == fq2_res_list[21]
        assert fqb+fq2d == fq2_res_list[22]
        assert fqc+fq2d == fq2_res_list[23]

        assert fqa*fq2b == fq2_res_list[24]
        assert fqa*fq2c == fq2_res_list[25]
        assert fqa*fq2d == fq2_res_list[26]
        assert fqb*fq2c == fq2_res_list[27]
        assert fqb*fq2d == fq2_res_list[28]
        assert fqc*fq2d == fq2_res_list[29]

        assert fqa-fq2b == fq2_res_list[30]
        assert fqa-fq2c == fq2_res_list[31]
        assert fqa-fq2d == fq2_res_list[32]
        assert fqb-fq2c == fq2_res_list[33]
        assert fqb-fq2d == fq2_res_list[34]
        assert fqc-fq2d == fq2_res_list[35]

        # fq2 op fq
        assert fq2a+fqb == fq2_res_list[36]
        assert fq2a+fqc == fq2_res_list[37]
        assert fq2a+fqd == fq2_res_list[38]
        assert fq2b+fqc == fq2_res_list[39]
        assert fq2b+fqd == fq2_res_list[40]
        assert fq2c+fqd == fq2_res_list[41]

        assert fq2a*fqb == fq2_res_list[42]
        assert fq2a*fqc == fq2_res_list[43]
        assert fq2a*fqd == fq2_res_list[44]
        assert fq2b*fqc == fq2_res_list[45]
        assert fq2b*fqd == fq2_res_list[46]
        assert fq2c*fqd == fq2_res_list[47]

        assert fq2a-fqb == fq2_res_list[48]
        assert fq2a-fqc == fq2_res_list[49]
        assert fq2a-fqd == fq2_res_list[50]
        assert fq2b-fqc == fq2_res_list[51]
        assert fq2b-fqd == fq2_res_list[52]
        assert fq2c-fqd == fq2_res_list[53]

        # int op fq2
        assert a+fq2b == fq2_res_list[18]
        assert a+fq2c == fq2_res_list[19]
        assert a+fq2d == fq2_res_list[20]
        assert b+fq2c == fq2_res_list[21]
        assert b+fq2d == fq2_res_list[22]
        assert c+fq2d == fq2_res_list[23]

        assert a*fq2b == fq2_res_list[24]
        assert a*fq2c == fq2_res_list[25]
        assert a*fq2d == fq2_res_list[26]
        assert b*fq2c == fq2_res_list[27]
        assert b*fq2d == fq2_res_list[28]
        assert c*fq2d == fq2_res_list[29]

        assert a-fq2b == fq2_res_list[30]
        assert a-fq2c == fq2_res_list[31]
        assert a-fq2d == fq2_res_list[32]
        assert b-fq2c == fq2_res_list[33]
        assert b-fq2d == fq2_res_list[34]
        assert c-fq2d == fq2_res_list[35]

        # fq2 op int
        assert fq2a+b == fq2_res_list[36]
        assert fq2a+c == fq2_res_list[37]
        assert fq2a+d == fq2_res_list[38]
        assert fq2b+c == fq2_res_list[39]
        assert fq2b+d == fq2_res_list[40]
        assert fq2c+d == fq2_res_list[41]

        assert fq2a*b == fq2_res_list[42]
        assert fq2a*c == fq2_res_list[43]
        assert fq2a*d == fq2_res_list[44]
        assert fq2b*c == fq2_res_list[45]
        assert fq2b*d == fq2_res_list[46]
        assert fq2c*d == fq2_res_list[47]

        assert fq2a-b == fq2_res_list[48]
        assert fq2a-c == fq2_res_list[49]
        assert fq2a-d == fq2_res_list[50]
        assert fq2b-c == fq2_res_list[51]
        assert fq2b-d == fq2_res_list[52]
        assert fq2c-d == fq2_res_list[53]

        # op fq2
        assert -fq2a == fq2_res_list[54]
        assert -(-fq2a) == fq2a
        assert -fq2b == fq2_res_list[55]
        assert -(-fq2b) == fq2b
        assert ~fq2a == fq2_res_list[56]
        assert ~(~fq2a) == fq2a
        assert ~fq2b == fq2_res_list[57]
        assert ~(~fq2b) == fq2b

    def test_fq6(self):
        a = int(fq_list[0]); b = int(fq_list[1]);
        c = int(fq_list[2]); d = int(fq_list[3]);
        fqa = fq_list[0]; fqb = fq_list[1];
        fqc = fq_list[2]; fqd = fq_list[3];
        fq2a = fq2_list[0]; fq2b = fq2_list[1]
        fq2c = fq2_list[2]; fq2d = fq2_list[3];
        fq6a = fq6_list[0]; fq6b = fq6_list[1]
        fq6c = fq6_list[2]; fq6d = fq6_list[3];

        # fq6 op fq6
        assert fq6a+fq6b == fq6_res_list[0]
        assert fq6a+fq6c == fq6_res_list[1]
        assert fq6a+fq6d == fq6_res_list[2]
        assert fq6b+fq6c == fq6_res_list[3]
        assert fq6b+fq6d == fq6_res_list[4]
        assert fq6c+fq6d == fq6_res_list[5]

        assert fq6a*fq6b == fq6_res_list[6]
        assert fq6a*fq6c == fq6_res_list[7]
        assert fq6a*fq6d == fq6_res_list[8]
        assert fq6b*fq6c == fq6_res_list[9]
        assert fq6b*fq6d == fq6_res_list[10]
        assert fq6c*fq6d == fq6_res_list[11]

        assert fq6a-fq6b == fq6_res_list[12]
        assert fq6a-fq6c == fq6_res_list[13]
        assert fq6a-fq6d == fq6_res_list[14]
        assert fq6b-fq6c == fq6_res_list[15]
        assert fq6b-fq6d == fq6_res_list[16]
        assert fq6c-fq6d == fq6_res_list[17]

        # fq2 op fq6
        assert fq2a+fq6b == fq6_res_list[18]
        assert fq2a+fq6c == fq6_res_list[19]
        assert fq2a+fq6d == fq6_res_list[20]
        assert fq2b+fq6c == fq6_res_list[21]
        assert fq2b+fq6d == fq6_res_list[22]
        assert fq2c+fq6d == fq6_res_list[23]

        assert fq2a*fq6b == fq6_res_list[24]
        assert fq2a*fq6c == fq6_res_list[25]
        assert fq2a*fq6d == fq6_res_list[26]
        assert fq2b*fq6c == fq6_res_list[27]
        assert fq2b*fq6d == fq6_res_list[28]
        assert fq2c*fq6d == fq6_res_list[29]

        assert fq2a-fq6b == fq6_res_list[30]
        assert fq2a-fq6c == fq6_res_list[31]
        assert fq2a-fq6d == fq6_res_list[32]
        assert fq2b-fq6c == fq6_res_list[33]
        assert fq2b-fq6d == fq6_res_list[34]
        assert fq2c-fq6d == fq6_res_list[35]

        # fq6 op fq2
        assert fq6a+fq2b == fq6_res_list[36]
        assert fq6a+fq2c == fq6_res_list[37]
        assert fq6a+fq2d == fq6_res_list[38]
        assert fq6b+fq2c == fq6_res_list[39]
        assert fq6b+fq2d == fq6_res_list[40]
        assert fq6c+fq2d == fq6_res_list[41]

        assert fq6a*fq2b == fq6_res_list[42]
        assert fq6a*fq2c == fq6_res_list[43]
        assert fq6a*fq2d == fq6_res_list[44]
        assert fq6b*fq2c == fq6_res_list[45]
        assert fq6b*fq2d == fq6_res_list[46]
        assert fq6c*fq2d == fq6_res_list[47]

        assert fq6a-fq2b == fq6_res_list[48]
        assert fq6a-fq2c == fq6_res_list[49]
        assert fq6a-fq2d == fq6_res_list[50]
        assert fq6b-fq2c == fq6_res_list[51]
        assert fq6b-fq2d == fq6_res_list[52]
        assert fq6c-fq2d == fq6_res_list[53]

        # fq op fq6
        assert fqa+fq6b == fq6_res_list[54]
        assert fqa+fq6c == fq6_res_list[55]
        assert fqa+fq6d == fq6_res_list[56]
        assert fqb+fq6c == fq6_res_list[57]
        assert fqb+fq6d == fq6_res_list[58]
        assert fqc+fq6d == fq6_res_list[59]

        assert fqa*fq6b == fq6_res_list[60]
        assert fqa*fq6c == fq6_res_list[61]
        assert fqa*fq6d == fq6_res_list[62]
        assert fqb*fq6c == fq6_res_list[63]
        assert fqb*fq6d == fq6_res_list[64]
        assert fqc*fq6d == fq6_res_list[65]

        assert fqa-fq6b == fq6_res_list[66]
        assert fqa-fq6c == fq6_res_list[67]
        assert fqa-fq6d == fq6_res_list[68]
        assert fqb-fq6c == fq6_res_list[69]
        assert fqb-fq6d == fq6_res_list[70]
        assert fqc-fq6d == fq6_res_list[71]

        # fq6 op fq
        assert fq6a+fqb == fq6_res_list[72]
        assert fq6a+fqc == fq6_res_list[73]
        assert fq6a+fqd == fq6_res_list[74]
        assert fq6b+fqc == fq6_res_list[75]
        assert fq6b+fqd == fq6_res_list[76]
        assert fq6c+fqd == fq6_res_list[77]

        assert fq6a*fqb == fq6_res_list[78]
        assert fq6a*fqc == fq6_res_list[79]
        assert fq6a*fqd == fq6_res_list[80]
        assert fq6b*fqc == fq6_res_list[81]
        assert fq6b*fqd == fq6_res_list[82]
        assert fq6c*fqd == fq6_res_list[83]

        assert fq6a-fqb == fq6_res_list[84]
        assert fq6a-fqc == fq6_res_list[85]
        assert fq6a-fqd == fq6_res_list[86]
        assert fq6b-fqc == fq6_res_list[87]
        assert fq6b-fqd == fq6_res_list[88]
        assert fq6c-fqd == fq6_res_list[89]

        # int op fq6
        assert a+fq6b == fq6_res_list[54]
        assert a+fq6c == fq6_res_list[55]
        assert a+fq6d == fq6_res_list[56]
        assert b+fq6c == fq6_res_list[57]
        assert b+fq6d == fq6_res_list[58]
        assert c+fq6d == fq6_res_list[59]

        assert a*fq6b == fq6_res_list[60]
        assert a*fq6c == fq6_res_list[61]
        assert a*fq6d == fq6_res_list[62]
        assert b*fq6c == fq6_res_list[63]
        assert b*fq6d == fq6_res_list[64]
        assert c*fq6d == fq6_res_list[65]

        assert a-fq6b == fq6_res_list[66]
        assert a-fq6c == fq6_res_list[67]
        assert a-fq6d == fq6_res_list[68]
        assert b-fq6c == fq6_res_list[69]
        assert b-fq6d == fq6_res_list[70]
        assert c-fq6d == fq6_res_list[71]

        # fq6 op int
        assert fq6a+b == fq6_res_list[72]
        assert fq6a+c == fq6_res_list[73]
        assert fq6a+d == fq6_res_list[74]
        assert fq6b+c == fq6_res_list[75]
        assert fq6b+d == fq6_res_list[76]
        assert fq6c+d == fq6_res_list[77]

        assert fq6a*b == fq6_res_list[78]
        assert fq6a*c == fq6_res_list[79]
        assert fq6a*d == fq6_res_list[80]
        assert fq6b*c == fq6_res_list[81]
        assert fq6b*d == fq6_res_list[82]
        assert fq6c*d == fq6_res_list[83]

        assert fq6a-b == fq6_res_list[84]
        assert fq6a-c == fq6_res_list[85]
        assert fq6a-d == fq6_res_list[86]
        assert fq6b-c == fq6_res_list[87]
        assert fq6b-d == fq6_res_list[88]
        assert fq6c-d == fq6_res_list[89]

        # op fq6
        assert -fq6a == fq6_res_list[90]
        assert -(-fq6a) == fq6a
        assert -fq6b == fq6_res_list[91]
        assert -(-fq6b) == fq6b
        assert ~fq6a == fq6_res_list[92]
        assert ~(~fq6a) == fq6a
        assert ~fq6b == fq6_res_list[93]
        assert ~(~fq6b) == fq6b

    def test_fq12(self):
        a = int(fq_list[0]); b = int(fq_list[1]);
        c = int(fq_list[2]); d = int(fq_list[3]);
        fqa = fq_list[0]; fqb = fq_list[1];
        fqc = fq_list[2]; fqd = fq_list[3];
        fq2a = fq2_list[0]; fq2b = fq2_list[1]
        fq2c = fq2_list[2]; fq2d = fq2_list[3];
        fq6a = fq6_list[0]; fq6b = fq6_list[1]
        fq6c = fq6_list[2]; fq6d = fq6_list[3];
        fq12a = fq12_list[0]; fq12b = fq12_list[1]
        fq12c = fq12_list[2]; fq12d = fq12_list[3];

        # fq12 op fq12
        assert fq12a+fq12b == fq12_res_list[0]
        assert fq12a+fq12c == fq12_res_list[1]
        assert fq12a+fq12d == fq12_res_list[2]
        assert fq12b+fq12c == fq12_res_list[3]
        assert fq12b+fq12d == fq12_res_list[4]
        assert fq12c+fq12d == fq12_res_list[5]

        assert fq12a*fq12b == fq12_res_list[6]
        assert fq12a*fq12c == fq12_res_list[7]
        assert fq12a*fq12d == fq12_res_list[8]
        assert fq12b*fq12c == fq12_res_list[9]
        assert fq12b*fq12d == fq12_res_list[10]
        assert fq12c*fq12d == fq12_res_list[11]

        assert fq12a-fq12b == fq12_res_list[12]
        assert fq12a-fq12c == fq12_res_list[13]
        assert fq12a-fq12d == fq12_res_list[14]
        assert fq12b-fq12c == fq12_res_list[15]
        assert fq12b-fq12d == fq12_res_list[16]
        assert fq12c-fq12d == fq12_res_list[17]

        # fq6 op fq12
        assert fq6a+fq12b == fq12_res_list[18]
        assert fq6a+fq12c == fq12_res_list[19]
        assert fq6a+fq12d == fq12_res_list[20]
        assert fq6b+fq12c == fq12_res_list[21]
        assert fq6b+fq12d == fq12_res_list[22]
        assert fq6c+fq12d == fq12_res_list[23]

        assert fq6a*fq12b == fq12_res_list[24]
        assert fq6a*fq12c == fq12_res_list[25]
        assert fq6a*fq12d == fq12_res_list[26]
        assert fq6b*fq12c == fq12_res_list[27]
        assert fq6b*fq12d == fq12_res_list[28]
        assert fq6c*fq12d == fq12_res_list[29]

        assert fq6a-fq12b == fq12_res_list[30]
        assert fq6a-fq12c == fq12_res_list[31]
        assert fq6a-fq12d == fq12_res_list[32]
        assert fq6b-fq12c == fq12_res_list[33]
        assert fq6b-fq12d == fq12_res_list[34]
        assert fq6c-fq12d == fq12_res_list[35]

        # fq12 op fq6
        assert fq12a+fq6b == fq12_res_list[36]
        assert fq12a+fq6c == fq12_res_list[37]
        assert fq12a+fq6d == fq12_res_list[38]
        assert fq12b+fq6c == fq12_res_list[39]
        assert fq12b+fq6d == fq12_res_list[40]
        assert fq12c+fq6d == fq12_res_list[41]

        assert fq12a*fq6b == fq12_res_list[42]
        assert fq12a*fq6c == fq12_res_list[43]
        assert fq12a*fq6d == fq12_res_list[44]
        assert fq12b*fq6c == fq12_res_list[45]
        assert fq12b*fq6d == fq12_res_list[46]
        assert fq12c*fq6d == fq12_res_list[47]

        assert fq12a-fq6b == fq12_res_list[48]
        assert fq12a-fq6c == fq12_res_list[49]
        assert fq12a-fq6d == fq12_res_list[50]
        assert fq12b-fq6c == fq12_res_list[51]
        assert fq12b-fq6d == fq12_res_list[52]
        assert fq12c-fq6d == fq12_res_list[53]

        # fq2 op fq12
        assert fq2a+fq12b == fq12_res_list[54]
        assert fq2a+fq12c == fq12_res_list[55]
        assert fq2a+fq12d == fq12_res_list[56]
        assert fq2b+fq12c == fq12_res_list[57]
        assert fq2b+fq12d == fq12_res_list[58]
        assert fq2c+fq12d == fq12_res_list[59]

        assert fq2a*fq12b == fq12_res_list[60]
        assert fq2a*fq12c == fq12_res_list[61]
        assert fq2a*fq12d == fq12_res_list[62]
        assert fq2b*fq12c == fq12_res_list[63]
        assert fq2b*fq12d == fq12_res_list[64]
        assert fq2c*fq12d == fq12_res_list[65]

        assert fq2a-fq12b == fq12_res_list[66]
        assert fq2a-fq12c == fq12_res_list[67]
        assert fq2a-fq12d == fq12_res_list[68]
        assert fq2b-fq12c == fq12_res_list[69]
        assert fq2b-fq12d == fq12_res_list[70]
        assert fq2c-fq12d == fq12_res_list[71]

        # fq12 op fq2
        assert fq12a+fq2b == fq12_res_list[72]
        assert fq12a+fq2c == fq12_res_list[73]
        assert fq12a+fq2d == fq12_res_list[74]
        assert fq12b+fq2c == fq12_res_list[75]
        assert fq12b+fq2d == fq12_res_list[76]
        assert fq12c+fq2d == fq12_res_list[77]

        assert fq12a*fq2b == fq12_res_list[78]
        assert fq12a*fq2c == fq12_res_list[79]
        assert fq12a*fq2d == fq12_res_list[80]
        assert fq12b*fq2c == fq12_res_list[81]
        assert fq12b*fq2d == fq12_res_list[82]
        assert fq12c*fq2d == fq12_res_list[83]

        assert fq12a-fq2b == fq12_res_list[84]
        assert fq12a-fq2c == fq12_res_list[85]
        assert fq12a-fq2d == fq12_res_list[86]
        assert fq12b-fq2c == fq12_res_list[87]
        assert fq12b-fq2d == fq12_res_list[88]
        assert fq12c-fq2d == fq12_res_list[89]

        # fq op fq12
        assert fqa+fq12b == fq12_res_list[90]
        assert fqa+fq12c == fq12_res_list[91]
        assert fqa+fq12d == fq12_res_list[92]
        assert fqb+fq12c == fq12_res_list[93]
        assert fqb+fq12d == fq12_res_list[94]
        assert fqc+fq12d == fq12_res_list[95]

        assert fqa*fq12b == fq12_res_list[96]
        assert fqa*fq12c == fq12_res_list[97]
        assert fqa*fq12d == fq12_res_list[98]
        assert fqb*fq12c == fq12_res_list[99]
        assert fqb*fq12d == fq12_res_list[100]
        assert fqc*fq12d == fq12_res_list[101]

        assert fqa-fq12b == fq12_res_list[102]
        assert fqa-fq12c == fq12_res_list[103]
        assert fqa-fq12d == fq12_res_list[104]
        assert fqb-fq12c == fq12_res_list[105]
        assert fqb-fq12d == fq12_res_list[106]
        assert fqc-fq12d == fq12_res_list[107]

        # fq12 op fq
        assert fq12a+fqb == fq12_res_list[108]
        assert fq12a+fqc == fq12_res_list[109]
        assert fq12a+fqd == fq12_res_list[110]
        assert fq12b+fqc == fq12_res_list[111]
        assert fq12b+fqd == fq12_res_list[112]
        assert fq12c+fqd == fq12_res_list[113]

        assert fq12a*fqb == fq12_res_list[114]
        assert fq12a*fqc == fq12_res_list[115]
        assert fq12a*fqd == fq12_res_list[116]
        assert fq12b*fqc == fq12_res_list[117]
        assert fq12b*fqd == fq12_res_list[118]
        assert fq12c*fqd == fq12_res_list[119]

        assert fq12a-fqb == fq12_res_list[120]
        assert fq12a-fqc == fq12_res_list[121]
        assert fq12a-fqd == fq12_res_list[122]
        assert fq12b-fqc == fq12_res_list[123]
        assert fq12b-fqd == fq12_res_list[124]
        assert fq12c-fqd == fq12_res_list[125]

        # int op fq12
        assert a+fq12b == fq12_res_list[90]
        assert a+fq12c == fq12_res_list[91]
        assert a+fq12d == fq12_res_list[92]
        assert b+fq12c == fq12_res_list[93]
        assert b+fq12d == fq12_res_list[94]
        assert c+fq12d == fq12_res_list[95]

        assert a*fq12b == fq12_res_list[96]
        assert a*fq12c == fq12_res_list[97]
        assert a*fq12d == fq12_res_list[98]
        assert b*fq12c == fq12_res_list[99]
        assert b*fq12d == fq12_res_list[100]
        assert c*fq12d == fq12_res_list[101]

        assert a-fq12b == fq12_res_list[102]
        assert a-fq12c == fq12_res_list[103]
        assert a-fq12d == fq12_res_list[104]
        assert b-fq12c == fq12_res_list[105]
        assert b-fq12d == fq12_res_list[106]
        assert c-fq12d == fq12_res_list[107]

        # fq12 op int
        assert fq12a+b == fq12_res_list[108]
        assert fq12a+c == fq12_res_list[109]
        assert fq12a+d == fq12_res_list[110]
        assert fq12b+c == fq12_res_list[111]
        assert fq12b+d == fq12_res_list[112]
        assert fq12c+d == fq12_res_list[113]

        assert fq12a*b == fq12_res_list[114]
        assert fq12a*c == fq12_res_list[115]
        assert fq12a*d == fq12_res_list[116]
        assert fq12b*c == fq12_res_list[117]
        assert fq12b*d == fq12_res_list[118]
        assert fq12c*d == fq12_res_list[119]

        assert fq12a-b == fq12_res_list[120]
        assert fq12a-c == fq12_res_list[121]
        assert fq12a-d == fq12_res_list[122]
        assert fq12b-c == fq12_res_list[123]
        assert fq12b-d == fq12_res_list[124]
        assert fq12c-d == fq12_res_list[125]

        # op fq12
        assert -fq12a == fq12_res_list[126]
        assert -(-fq12a) == fq12a
        assert -fq12b == fq12_res_list[127]
        assert -(-fq12b) == fq12b
        assert ~fq12a == fq12_res_list[128]
        assert ~(~fq12a) == fq12a
        assert ~fq12b == fq12_res_list[129]
        assert ~(~fq6b) == fq6b


"""
Copyright 2018 Chia Network Inc

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
