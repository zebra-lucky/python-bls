from typing import List

from .ec import AffinePoint, default_ec, generator_Fq
from .fields import Fq, Fq2
from .signature import Signature


class Threshold:
    """
    Utility functions for threshold signatures.

    The end result of initializing a T of N Joint-Feldman scheme is
    that each of N players has a secret share (private key) they can
    use to sign messages; and there is a master secret and public key.

    The signatures of any T players (along with the indices of those
    players) can be combined to form a valid signature for that master
    key pair.

    To initialize a T of N threshold key under a Joint-Feldman scheme:

    1. Each player calls PrivateKey.new_threshold(T, N)
       to create a secret key, commitment to a polynomial, and
       secret fragments.
       They send everyone the commitment, and send player j
       secret_fragments[j].

    2. Each player calls Threshold.verify_secret_fragment
       on all secret fragments they receive.  If any verify False,
       they complain to abort the scheme.  (Note that repeatedly
       aborting or 'speaking' last, can bias the master public key.)

    3. Each player computes the shared, master public key:
       master_pubkey = BLS.aggregate_pub_keys(
           [PublicKey.from_g1(cpoly[0].to_jacobian())
            for cpoly in commitments],
           False)

       They also create their secret share from all secret
       fragments received (now verified):
       secret_share = BLS.aggregate_priv_keys(
           [PrivateKey(frag) for frag in share_fragments],
           None,
           False)

    4. Player P may create a signature share with respect to T players:
       sig_share = secret_share.sign_threshold(message, P, players),
       where 'players' is a list of T different player indices
       participating.

       These signature shares can be combined to sign the message:
       signature = BLS.aggregate_sigs_simple(sig_shares).
    """


    @staticmethod
    def lagrange_coeffs_at_zero(X: List[int], ec=default_ec) -> List[Fq]:
        """
        We have k+1 integers X[i], all less than ec.n and non-zero.
        The points (X[i], P(X[i])) interpolate into P(X), a degree k polynomial.
        Returns coefficients L_i such that P(0) = sum L_i * P(X[i]).
        """
        N = len(X)

        # Check all x values are different, less than ec.n, and non-zero.
        assert len(set(X)) == N and all(0 != x < ec.n for x in X)

        def weight(j):
            ans = Fq(ec.n, 1)
            for i in range(N):
                if i != j:
                    ans *= Fq(ec.n, X[j] - X[i])
            return ~ans

        # Using the second barycentric form,
        # P(0) = (sum_j (y_j * w_j / x_j)) / (sum_j w_j/x_j)
        # If desired, the weights can be precomputed.

        ans = []
        denominator = Fq(ec.n, 0)
        for j in range(N):
            shift = weight(j) * ~Fq(ec.n, -X[j])
            ans.append(shift)
            denominator += shift
        denominator = ~denominator
        for i in range(len(ans)):
            ans[i] *= denominator
        return ans


    @staticmethod
    def interpolate_at_zero(X: List[int], Y: List[Fq], ec=default_ec) -> Fq:
        """
        The k+1 points (X[i], Y[i]) interpolate into P(X),
        a degree k polynomial.
        Returns P(0).
        """
        ans = Fq(ec.n, 0)
        for lamb, y in zip(Threshold.lagrange_coeffs_at_zero(X, ec), Y):
            ans += lamb * y
        return ans


    @staticmethod
    def verify_secret_fragment(T: int, secret_fragment: Fq,
                               player: int, commitment: List[AffinePoint],
                               ec=default_ec) -> bool:
        """
        You are player, and have received a secret share fragment,
        claimed to be shares[i] = P(player) from a polynomial P
        with the given commitment.

        Return True if the share given to you is correct wrt that commitment.
        """
        assert len(commitment) == T
        assert secret_fragment != 0
        assert player != 0

        g1 = generator_Fq(ec)
        lhs = g1 * secret_fragment
        rhs = commitment[0]
        for k in range(1, len(commitment)):
            rhs += commitment[k] * pow(player, k, ec.n)

        return lhs == rhs

    @staticmethod
    def aggregate_unit_sigs(signatures: List[Signature], players: List[int],
            T: int, ec=default_ec) -> Signature:

        lambs = Threshold.lagrange_coeffs_at_zero(players, ec)
        agg = (AffinePoint(Fq2.zero(ec.q), Fq2.zero(ec.q), True, ec)
               .to_jacobian())
        for i, sig in enumerate(signatures):
            agg += sig.value * lambs[i]
        return Signature.from_g2(agg)

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
