#!/usr/bin/env python3
"""Acceptance checks for replaying the frozen certificate.

These tests exercise rejection of altered proof inputs and stored claims;
they do not duplicate the certificate's mathematical formulas.
Run from any directory with: python path/to/test_certificate.py
"""
import copy
import json
from pathlib import Path
import unittest

import certify_core


class CertificateReplayAcceptanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = Path(__file__).resolve().parent / "results" / "certificate_k3_D4.json"
        cls.reference = json.loads(fixture.read_text(encoding="utf-8"))

    def changed_certificate(self):
        return copy.deepcopy(self.reference)

    def assert_rejected(self, document):
        with self.assertRaises(certify_core.CertificateFailure):
            certify_core.verify_certificate(document)

    def test_unmodified_frozen_certificate_is_accepted(self):
        result = certify_core.verify_certificate(self.reference)
        self.assertIs(result["verified"], True)
        self.assertEqual((result["k"], result["D"]), (3, 4))

    def test_root_enclosure_must_exhibit_the_required_sign(self):
        document = self.changed_certificate()
        document["root"]["balance_at_lower"]["lower"] = "-1"
        self.assert_rejected(document)

    def test_exterior_cover_cannot_omit_an_interval(self):
        document = self.changed_certificate()
        del document["exterior"][0]["cells"][0]
        self.assert_rejected(document)

    def test_changing_a_chernoff_tilt_requires_new_valid_bounds(self):
        document = self.changed_certificate()
        document["exterior"][0]["cells"][0]["tilt_A"] = "0"
        self.assert_rejected(document)

    def test_forged_variance_enclosure_is_rejected(self):
        document = self.changed_certificate()
        document["profiles"]["V_A"] = {"lower": "0", "upper": "0"}
        self.assert_rejected(document)

    def test_degree_profile_cannot_be_truncated(self):
        document = self.changed_certificate()
        document["profiles"]["A"].pop()
        self.assert_rejected(document)

    def test_published_numerical_target_cannot_be_relaxed(self):
        document = self.changed_certificate()
        document["published_upper_targets"]["exterior_U"] = "0"
        self.assert_rejected(document)

    def test_saved_curvature_bound_must_remain_strictly_negative(self):
        document = self.changed_certificate()
        document["central"]["cells"][0]["gamma_second_derivative"]["upper"] = "0"
        self.assert_rejected(document)

    def test_recorded_interval_runtime_must_match(self):
        document = self.changed_certificate()
        document["mpmath_version"] = "unverified-version"
        self.assert_rejected(document)

    def test_summary_cannot_claim_a_stronger_unverified_bound(self):
        document = self.changed_certificate()
        document["summary"]["exterior_U_upper"] = "-10"
        self.assert_rejected(document)


if __name__ == "__main__":
    unittest.main(verbosity=2)
